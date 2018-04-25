################################################################################
# Author: W. Jenny Shi    
# 
# About: This script incorporates subsampling and possibly random partition of  
#   data features to SmCCA and SsCCA. Here we assume two data types (e.g. miRNA
#   and genes) and one quantitative phenotype measured for the same subjects.
#
################################################################################

library(pbapply)
library(igraph)


getRobustPseudoWeights <- function(X1, X2, Traits, Lambda1, Lambda2, 
                                   s1 = 0.9, s2 = 1, NoTrait = FALSE, 
                                   FilterByTrait = FALSE, Bipartite = FALSE, 
                                   SubsamplingNum = 100, PartitionNum = 100, 
                                   CCcoef = NULL, trace = FALSE){
  # Compute aggregated canonical weights.
  #
  # X1: An n by p1 mRNA expression matrix. 
  # X2: An n by p2 miRNA expression matrix.
  # Traits: An n by k trait data for the same samples (k >=1).
  # Lambda1, Lambda2: LASSO pentalty parameters, need to be between 0 and 1. 
  # s1: Percent of mRNA features to be included.
  # s2: Percent of miRNA features to be included.
  # NoTrait: Logical. Whether trait information is provided.
  # FilterByTrait: Logical. Whether only the features with highest correlation 
  #   to traits will be assigned nonzero weights. The top 80% features are reserved. 
  # Bipartite: Logical. Whether to include random partition.
  # SubsamplingNum: Number of subsampling features.
  # PartitionNum: Number of random partition.
  # CCcoef: Optional coefficients for the pairwise canonical correlations (CC).  
  #   If CCcoef = NULL (default), then the objective function is the total sum  
  #   of all pairwise CC. It can also be a coefficient vector that follows the 
  #   column order of combn(K, 2).
  # trace: Logical. Whether to display CCA algorithm trace.
  
  if(min(s1, s2) == 0){
    stop("Subsampling proprotion needs to be greater than zero.")
  }else{
    if((abs(s1 - 0.5) > 0.5) | (abs(s2 - 0.5) > 0.5)){
      stop("Subsampling proportions can not exceed one.")}
  }
  
  if(abs(Lambda1 - 0.5) > 0.5 | Lambda1 == 0){
    stop("Invalid penalty parameter. Lambda1 needs to be between zero and one.")}
  if(abs(Lambda2 - 0.5) > 0.5 | Lambda2 == 0){
    stop("Invalid penalty parameter. Lambda2 needs to be between zero and one.")}
  
  p1 <- ncol(X1); p2 <- ncol(X2); p <- p1 + p2
  p1.sub <- ceiling(s1 * p1);   p2.sub <- ceiling(s2 * p2)
  X <- cbind(X1, X2)
  
  beta <- pbsapply(1:SubsamplingNum, function(x){
    # Subsample features	
    samp1 <- sort(sample(1:p1, p1.sub, replace = FALSE))
    samp2 <- sort(sample(1:p2, p2.sub, replace = FALSE))
    
    if(Bipartite){ # If only subsample features.
      x1.par <- scale(X1[ , samp1], center = TRUE, scale = TRUE)
      x2.par <- scale(X2[ , samp2], center = TRUE, scale = TRUE)
      
      out <- getCCAout(x1.par, x2.par, Traits, Lambda1, Lambda2, 
                       NoTrait = NoTrait, FilterByTrait = FilterByTrait, 
                       trace = trace, CCcoef = CCcoef)
      
      w <- rep(0, p)
      w[samp1] <- out$u
      w[samp2 + p1] <- out$v
      coeff.avg <- w
      
    }else{ # If perform both subsampling and random partition.
      samp <- c(samp1, samp2 + p1) 
      subSamp <- scale(X[ , samp], center = TRUE, scale = TRUE)
      p.sub <- ncol(subSamp); p.par <- p.sub %/% 2
      
      coeff <- pbsapply(1:PartitionNum, function(y){
        gp1 <- sample(1:p.sub, p.par, replace = FALSE)
        x1.par <- subSamp[ , gp1]
        gp2 <- sample((1:p.sub)[-gp1], p.sub-p.par, replace = FALSE)
        x2.par <- subSamp[ , gp2]
        
        out <- getCCAout(x1.par, x2.par, Traits, Lambda1, Lambda2,
                         NoTrait = NoTrait, FilterByTrait = FilterByTrait, 
                         trace = trace, CCcoef = CCcoef)
        
        w <- rep(0, p)
        w[samp[gp1]] <- out$u
        w[samp[gp2]] <- out$v
        return(w)
      })
      
      coeff.avg <- apply(abs(coeff), 1, mean)
      
    }
    
    return(coeff.avg)
  })
  
  return(beta)
}


computeEdgeWeightMatrix <- function(W){
  # Compute edge weight matrix. It first take the outer product of absolute of 
  #   SCCA weights, then average over the colnums of W and normalized by the 
  #   maximum value of W. 
  # 
  # W: a weight vector or matrix from SCCA. If W is a matrix, then each column 
  #   corresponds to one weight vector and the resulting edge weights will be 
  #   averaged. 
  
  if(is.null(dim(W))){
    Abar <- abs(W) %o% abs(W)
  }else{
    b <- nrow(W)
    Abar <- matrix(0, nrow = b, ncol = b)
    for(ind in 1:ncol(W)){
      w <- abs(W[ , ind])
      A <- w %o% w
      Abar <- Abar + A
    }
  }
  
  diag(Abar) <- 0
  Abar <- Abar/max(Abar)
  return(Abar)  
}


GetFuncGrp <- function(sccaFile, SaveResultFile, AbarLabel = NULL,
                       Cutoff = 1, SavePlotFile = NULL,
                       PlotCex = 0.6, SavePlotWidth = 16){
  # Report functional group obtained using hierarchical tree cut with the 
  #   complete linkage. 
  #
  # sccaFile: SCCA result file.
  # SaveResultFile: File name for the result to be saved as.
  # AbarLabel: Optional feature labels.
  # Cutoff: Tree cut height. Default at 1. 
  # SavePlotFile: File name for the hierarchical tree plot to be saved as.


  load(sccaFile)
  if(is.null(AbarLabel)) AbarLabel <- 1:nrow(Abar)
  
  rownames(Abar) <- colnames(Abar) <- AbarLabel
  hc <- hclust(as.dist(1 - Abar))
  
  if(is.null(SavePlotFile)){
    plot(hc, main = paste0("pen=(", l1, ",", l2,")"), cex = PlotCex)
  }else{
    pdf(width = SavePlotWidth, SavePlotFile)
    plot(hc, main = paste0("pen=(", l1, ",", l2,")"), cex = PlotCex)
    dev.off() 
  }
  
  cut.merge <- hc$merge[hc$height < Cutoff, ]
  lower.leaves <- sort(-cut.merge[cut.merge<0])

  save(hc, l1, l2, Cutoff, lower.leaves, file = SaveResultFile)
  return(lower.leaves)
}


Leaves2Module <- function(Leaves, GrpID){
  # Create module index list.
  # 
  # Leaves: indices for the nodes below the threshold in the hierarchical tree.
  # GrpID: group membership for all nodes.
  
  id <- GrpID[Leaves]
  M <- lapply(1:length(unique(id)), function(x){
    M.x <- Leaves[which(id == unique(id)[x])]
    return(M.x)
  })
  
  return(M)
}



PlotReducedNetwork <- function(EdgeMatrix, CorrMatrix, mRNAmiRNAmodules, 
                               ModuleIdx, p1, EdgeCut, FeatureLabel = NULL, 
                               AddCorrSign = TRUE, SaveFile = NULL, 
                               ShowType1Label = TRUE, ShowType2Label = TRUE, 
                               PlotTitle = "", NetLayout = "lgl",
                               ShowNodes = TRUE, 
                               VertexLabelCex = 1, VertexSize = 1){
  # Plot a trimmed module.
  # 
  # EdgeMatrix: The similarity matrix based on (average pseudo) canonical 
  #   weights. All entries are non-negative.
  # CorrMatrix: A pairwise correlation matrix that provides sign information 
  #   for the network.
  # mRNAmiRNAmodules: A list of full modules. It contains two data types.
  # ModuleIdx: The index of the full module to be plotted. 
  # AddCorrSign: Add signs to the network edges based on CorrMatrix.
  # p1: Number of feature for the first data type.
  # FeatureLabel: A long vector indicating feature names. 
  # SaveFile: File name for the output figure. 
  # PlotTitle: Title for the figure.
  # NetLayout: Graphical layout for the network. Possible options are "circle",
  #   "sphere" for 3D sphere, "fr" for Fruchterman-Reinhold, and "lgl" for the 
  #   LGL algorithm. Refer to igraph manual for more details on the layout 
  #   options.
  # ShowType1Label, ShowType2Label: Logical. Whether to show labels for either 
  #   data type.
  # ShowNodes: Logical. Whether to show nodes.
  # VertexLavelCex: Scaling factor for vertex label.
  # VertexSize: Size of the vertex.
  
  EdgeMatrix <- as.matrix(EdgeMatrix)
  p <- ncol(EdgeMatrix)
  if(is.null(FeatureLabel)){FeatureLabel <- 1:p}
  colnames(EdgeMatrix) <- rownames(EdgeMatrix) <- FeatureLabel[1:p]
  colnames(CorrMatrix) <- rownames(CorrMatrix) <- FeatureLabel[1:p]
    
  grp <- mRNAmiRNAmodules[[ModuleIdx]]
  grp.memb <- colnames(EdgeMatrix)[grp]
  M.node <- grp.memb

  # Trim the module by EdgeCut.
  M <- as.matrix(EdgeMatrix[M.node, M.node])
  if(AddCorrSign){M <- M * sign(CorrMatrix[M.node, M.node])}
  M[which(abs(M) < EdgeCut)] <- 0
  newM.node <- M.node[which(apply(abs(M), 1, max) > 0)]
  
  if(length(newM.node) == 0){
    print("No edge passes threshold.")
  }else{
    M <- M[newM.node, newM.node]
    allidx <- matrix(1:p, ncol = 1)
    rownames(allidx) <- rownames(EdgeMatrix)
    
    NodeInfo <- data.frame(id = newM.node, idx = allidx[newM.node, ])
    net <- graph_from_adjacency_matrix(M, weighted = TRUE, 
                                       diag = FALSE, mode = "undirected")
    
    # Define colors and shapes for vertex and label.
    k <- length(newM.node)
    type1 <- which(NodeInfo$idx <= p1)
    type2 <- which(NodeInfo$idx > p1)
    vcol <- rep("purple", k); vcol[type2] <- "dark orange"
    vshape <- rep("square", k); vshape[type2] <- "circle"
    lcol <- vcol
    # If not show nodes.
    if(!ShowNodes){vshape <- "none"}
    # If no label, assign a place holder.
    if(!ShowType1Label){newM.node[type1] <- " "} 
    if(!ShowType2Label){newM.node[type2] <- " "}
    
    
    # Define edge colors.
    ecol <- rep("gray80", ecount(net))
    ew <- abs(edge.attributes(net)$weight) * 5
    ecol[which(edge.attributes(net)$weight < 0)] <- "red"
    
    # Define network layout.
    if(NetLayout == "circle"){
      l <- layout_in_circle(net)
    }else if(NetLayout == "sphere"){
      l <- layout_on_sphere(net)
    }else if(NetLayout == "fr"){
      l <- layout_with_fr(net)
    }else if(NetLayout == "lgl"){
      l <- layout_with_lgl(net)
    }else{
      stop("Unrecognized NetLayout input. Acceptable options are 'circle',
           'sphere', 'fr', and 'lgl'.")
    }
    
    
    if(!is.null(SaveFile)){
      pdf(SaveFile)
      par(bg = "white")
      plot(net, vertex.color = vcol, vertex.shape = vshape, 
           vertex.label.cex = VertexLabelCex, layout = l, 
           vertex.size = VertexSize, vertex.label = newM.node,
           edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
           vertex.label.font = 2)
      title(PlotTitle)
      dev.off()
    }else{
      plot(net, vertex.color = vcol, vertex.shape = vshape, 
           vertex.label.cex = VertexLabelCex, layout = l, 
           vertex.size = VertexSize, vertex.label = newM.node,
           edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
           vertex.label.font = 2, main = PlotTitle)
    }
    
  }
  
}






#######################################################
# Internal functions called by getRobustPseudoWeights #
#######################################################

getCCAout <- function(X1, X2, Traits, Lambda1, Lambda2, CCcoef = NULL, 
                      NoTrait = FALSE, FilterByTrait = FALSE, trace = FALSE){
  # Compute CCA weights.
  # 
  # X1: An n by p1 mRNA expression matrix. 
  # X2: An n by p2 miRNA expression matrix.
  # Traits: An n by k trait data for the same samples (k >=1).
  # Lambda1, Lambda2: LASSO pentalty parameters, need to be between 0 and 1. 
  # CCcoef: A 3 by 1 vector indicating weights for each pairwise canonical 
  #   correlation.
  # NoTrait: Logical. Whether trait information is provided.
  # FilterByTrait: Logical. Whether only the features with highest correlation 
  #   to traits will be assigned nonzero weights. The top 80% features are reserved. 
  # trace: Logical. Whether to display CCA algorithm trace.
  
  if(abs(Lambda1 - 0.5) >= 0.5){
    stop("Invalid penalty parameter. Lambda1 needs to be between zero and one.")}
  if(abs(Lambda2 - 0.5) >= 0.5){
    stop("Invalid penalty parameter. Lambda2 needs to be between zero and one.")}
  
  k <- ncol(Traits)
  if(NoTrait){Traits <- NULL; FilterByTrait <- TRUE; k <- 1}
  
  if(FilterByTrait){
    if(k > 1){
      stop("'FilterByTrait == TRUE' only allows one trait at a time.")
    }else{
      out <- CCA(X1, X2, outcome = "quantitative", y = Traits, 
                 typex = "standard", typez = "standard", penaltyx = Lambda1, 
                 penaltyz = Lambda2, trace = trace)
    }
  }else{
    xlist <- list(x1 = X1, x2 = X2, y = Traits)
    L1 <- max(1, sqrt(ncol(X1)) * Lambda1)
    L2 <- max(1, sqrt(ncol(X2)) * Lambda2)
    out <- myMultiCCA(xlist, penalty = c(L1, L2, sqrt(ncol(Traits))), 
                      CCcoef = CCcoef, trace = trace)
    out$u <- out$ws[[1]]; out$v <- out$ws[[2]]
  }
  
  return(out)
}
