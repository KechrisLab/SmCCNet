################################################################################
# Author: W. Jenny Shi
#
# About: This script incorporates subsampling and possibly random partition of
#   data features to SmCCA and SsCCA. Here we assume two data types (e.g. miRNA
#   and genes) and one quantitative phenotype measured for the same subjects.
#
################################################################################


######## Package that we need ########
# You only need to install once #
# install.packages("Matrix", repos="http://cran.r-project.org")
# install.packages("pbapply", repos="http://cran.r-project.org")
# install.packages("igraph", repos="http://cran.r-project.org")

requireNamespace("Matrix", quietly = TRUE)
requireNamespace("pbapply", quietly = TRUE)
requireNamespace("igraph", quietly = TRUE)


# library(Matrix)
# library(pbapply)
# library(igraph)


################################################################################
### Apply sparse multiple canonical correlation analysis to omics feature
### subsamples

#' Compute canonical weights based on sparse multiple canonical correlations
#' (SmCCA), sparse supervised canonical correlations (SsCCA), or sparse
#' canonical correlations (SCCA).
#'
#' SmCCA and SsCCA take into account of a phenotype/trait. SmCCA maximizes the
#' total (weighted or unweighted) pairwise canonical weights between two omics
#' data type and the trait. It requires the trait to be quantitative. SsCCA
#' prioritizes omics features based on the trait, and assigns non-zero canonical
#' weights to features that are more correlated to the trait. SCCA does not use
#' any trait information for computing the canonical weights. All of these three
#' methods are included in the function, along with an omic feature subsampling
#' scheme.
#'
#'
#' @param X1 An n\eqn{\times}p1 data matrix (e.g. mRNA) with p1 features and n
#'   subjects.
#' @param X2 An n\eqn{\times}p2 data matrix (e.g. miRNA) with p2 features and n
#'   subjects.
#' @param Trait An n\eqn{\times}1 trait data matrix for the same n subjects.
#' @param Lambda1 LASSO pentalty parameter for X1, need to be between 0 and 1.
#' @param Lambda2 LASSO pentalty parameter for X2, need to be between 0 and 1.
#' @param s1 Proportion of mRNA features to be included.
#' @param s2 Proportion of miRNA features to be included.
#' @param NoTrait Logical. Whether trait information is provided.
#' @param FilterByTrait Logical. Whether only the top 80% features with highest
#'   correlation to the trait will be assigned nonzero weights.
#' @param Bipartite Logical. Whether to include random partition.
#' @param SubsamplingNum Number of feature subsamples. Larger number leads to
#'   more accurate results, but at a higher cost. We recommend to subsample at
#'   least 1000 times.
#' @param PartitionNum: Number of random partitions for each . This is only used if
#'   Bipartite = FALSE.
#' @param CCcoef Optional coefficients for the pairwise canonical correlations
#'   (CC). If CCcoef = NULL (default), then the objective function is the total
#'   sum of all pairwise CC. It can also be a coefficient vector that follows
#'   the column order of combn(K, 2).
#' @param trace Logical. Whether to display CCA algorithm trace.
#' @return A canonical weight matrix with p1+p2 rows. Each column is the
#'   canonical weights based on subsampled X1 and X2 features.The column number
#'   equals to SubsamplineNum.
#'
#' @examples
#' \donttest{
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#'
#' # SmCCA
#' W1 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = FALSE,
#'   Bipartite = FALSE, SubsamplingNum = 100, PartitionNum = 5,
#'   CCcoef = NULL, trace = FALSE)
#'
#' # SsCCA
#' W2 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = .05, Lambda2 = 0.5
#'   s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = TRUE,
#'   Bipartite = FALSE, SubsamplingNum = 100, PartitionNum = 5,
#'   CCcoef = NULL, trace = FALSE)
#'
#' # SCCA
#' W3 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = TRUE, Bipartite = FALSE,
#'   SubsamplingNum = 100, PartitionNum = 5, CCcoef = NULL, trace = FALSE)
#'   }
#'
#' @export
getRobustPseudoWeights <- function(X1, X2, Trait, Lambda1, Lambda2,
                                   s1 = 0.9, s2 = 1, NoTrait = FALSE,
                                   FilterByTrait = FALSE, Bipartite = FALSE,
                                   SubsamplingNum = 100, PartitionNum = 100,
                                   CCcoef = NULL, trace = FALSE){

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

      out <- getCCAout(x1.par, x2.par, Trait, Lambda1, Lambda2,
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

        out <- getCCAout(x1.par, x2.par, Trait, Lambda1, Lambda2,
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


################################################################################
### Aggregate pseudo canonical weights and compute the network similarity matrix
### for all omics features.

#' Compute the similarity matrix based on the outer products of absolute canonical
#' weights.
#'
#'
#' @param Ws A canonical weight vector or matrix. If Ws is a matrix, then each
#'   column corresponds to one weight vector.
#' @param FeatureLabel A 1\eqn{\times}p vector indicating feature names. If
#'   FeatureLabel = NULL (default), the feature names will be indices 1
#'   through p, where p is the total number of omics features.
#' @return A p\eqn{\times}p symmetric non-negative matrix.
#'
#' @examples
#' w <- matrix(rnorm(6), nrow = 3)
#' Ws <- apply(w, 2, function(x)return(x/sqrt(sum(x^2))))
#' abar <- getAbar(Ws)
#'
#' @export
getAbar <- function(Ws, FeatureLabel = NULL){

    if(is.null(dim(Ws))){
        Abar <- Matrix::Matrix(abs(Ws) %o% abs(Ws), sparse = TRUE)
    }else{
        b <- nrow(Ws)
        Abar <- matrix(0, nrow = b, ncol = b)
        for(ind in 1:ncol(Ws)){
            w <- abs(Ws[ , ind])
            A <- Matrix::Matrix(w %o% w, sparse = TRUE)
            Abar <- Abar + A
        }
    }

    diag(Abar) <- 0
    Abar <- Abar/max(Abar)

    if(is.null(colnames(Abar))){
        if(is.null(FeatureLabel)){FeatureLabel <- 1:nrow(Abar)}
        colnames(Abar) <- rownames(Abar) <- FeatureLabel
    }

    return(Abar)
}


################################################################################
### Extract multi-omics modules.

#' Apply a hierarchical tree cutting to the similarity matrix and extract
#' modules that contain both omics data types.
#'
#' @param Abar A similary matrix for both omics data types.
#' @param P1 Total number of features for the first omics data type.
#' @param CutHeight Height threshold for the hierarchical tree. Default is
#'   1-.1^10.
#' @param PlotTree Logical. Whether to create the hierarchical tree plot.
#' @return A list of multi-omics modules.
#'
#' @examples
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' abar <- getAbar(w)
#' modules <- getMultiOmicsModules(abar, P1 = 2, CutHeight = 0.5)
#'
#' @export
getMultiOmicsModules <- function(Abar, P1, CutHeight = 1-.1^10, PlotTree = TRUE){

    hc <- hclust(as.dist(1 - Abar))
    if(PlotTree){plot(hc)}
    cut.merge <- hc$merge[hc$height < CutHeight, ]
    lower.leaves <- sort(-cut.merge[cut.merge<0])

    grpID <- cutree(hc, h = CutHeight)
    id <- grpID[lower.leaves]
    M <- lapply(1:length(unique(id)), function(x){
        M.x <- lower.leaves[which(id == unique(id)[x])]
        return(M.x)
    })

    multiOmicsModule <- lapply(M, function(s){
        s.min <- min(s)
        s.max <- max(s)
        if(s.min <= P1 & s.max > P1)return(s)
    })

    if(length(multiOmicsModule) > 1){
        nullSet <- which(sapply(multiOmicsModule, is.null))
        if(length(nullSet) > 0){
            multiOmicsModule <- multiOmicsModule[-nullSet]
        }
    }

    return(multiOmicsModule)
}


################################################################################
### Visualize multi-omics subnetworks.

#' Plot multi-omics modules based on similarity matrix derived from pseudo
#' canonical weights and pairwise feature correlations.
#'
#' @param Abar A p\eqn{\times}p similary matrix for both omics data types
#'   based on pseudo canonical weights. All entries are non-negative.
#' @param CorrMatrix A p\eqn{\times}p correlation matrix that provides sign
#'   information for the network.
#' @param multiOmicsModule A list of multi-omics modules.
#' @param ModuleIdx Index for the module to be plotted. It can not exceed the
#'   length of multiOmicsModule.
#' @param P1 Total number of features for the first omics data type.
#' @param EdgeCut A numerical value between 0 and 1, indicating an edge
#'   threshold for the network. Any features (network nodes) without any edge
#'   strength that passes the threshold are excluded from the figure.
#' @param FeatureLabel A 1\eqn{\times}p vector indicating feature names. If
#'   FeatureLabel = NULL (default), the feature names will be indices 1
#'   through p, where p is the total number of omics features.
#' @param AddCorrSign Logical. Whether to add a positive or negative sign to
#'   each network edge based on pairwise feature correlations.
#' @param SaveFile A pdf file name for the figure output. If SaveFile = NULL
#'   (default), the figure will not be saved.
#' @param ShowType1Label Logical. Whether to label the network nodes for the
#'   first omics data type.
#' @param ShowType2Label Logical. Whether to label the network nodes for the
#'   second omics data type.
#' @param PlotTitle A title for the figure. Default is without any title.
#' @param NetLayout Graphical layout for the network. Possible options are
#'   "circle", "sphere" for 3D sphere, "fr" for Fruchterman-Reinhold, and "lgl"
#'   for the LGL algorithm. Refer to igraph manual for more details on the
#'   layout options.
#' @param ShowNodes Logical. Whether to show network nodes.
#' @param VertexLabelCex Scaling factor for the vertex labels.
#' @param VertexSize Size of the vertices.
#' @return A multi-omics network figure.
#'
#' @examples
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' abar <- getAbar(w)
#' modules <- getMultiOmicsModules(abar, P1 = 2, CutHeight = 0.5)
#' x <- cbind(X1[ ,1:2], X2[ , 1:3])
#' corr <- cor(x)
#'
#' plotMultiOmicsNetwork(abar, corr, modules, ModuleIdx = 1, P1 = 2,
#'   EdgeCut = 0)
#'
#' @export
plotMultiOmicsNetwork <- function(Abar, CorrMatrix, multiOmicsModule,
                               ModuleIdx, P1, EdgeCut, FeatureLabel = NULL,
                               AddCorrSign = TRUE, SaveFile = NULL,
                               ShowType1Label = TRUE, ShowType2Label = TRUE,
                               PlotTitle = "", NetLayout = "lgl",
                               ShowNodes = TRUE,
                               VertexLabelCex = 1, VertexSize = 1){

    p <- ncol(Abar)
    if(is.null(FeatureLabel)){FeatureLabel <- 1:p}
    colnames(Abar) <- rownames(Abar) <- FeatureLabel[1:p]
    colnames(CorrMatrix) <- rownames(CorrMatrix) <- FeatureLabel[1:p]

    grp <- multiOmicsModule[[ModuleIdx]]
    grp.memb <- colnames(Abar)[grp]
    M.node <- grp.memb

    # Trim the module by EdgeCut.
    M <- as.matrix(Abar[M.node, M.node])
    if(AddCorrSign){M <- M * sign(CorrMatrix[M.node, M.node])}
    M[which(abs(M) < EdgeCut)] <- 0
    newM.node <- M.node[which(apply(abs(M), 1, max) > 0)]

    if(length(newM.node) == 0){
        print("No edge passes threshold.")
    }else{
        M <- M[newM.node, newM.node]
        allidx <- matrix(1:p, ncol = 1)
        rownames(allidx) <- rownames(Abar)

        NodeInfo <- data.frame(id = newM.node, idx = allidx[newM.node, ])
        net <- graph_from_adjacency_matrix(M, weighted = TRUE,
                                           diag = FALSE, mode = "undirected")

        # Define colors and shapes for vertex and label.
        k <- length(newM.node)
        type1 <- which(NodeInfo$idx <= P1)
        type2 <- which(NodeInfo$idx > P1)
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
            grDevices::pdf(SaveFile)
            graphics::par(bg = "white")
            graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
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

# (INTERNAL)
getCCAout <- function(X1, X2, Trait, Lambda1, Lambda2, CCcoef = NULL,
                      NoTrait = FALSE, FilterByTrait = FALSE, trace = FALSE){
  # Compute CCA weights.
  #
  # X1: An n by p1 mRNA expression matrix.
  # X2: An n by p2 miRNA expression matrix.
  # Trait: An n by k trait data for the same samples (k >=1).
  # Lambda1, Lambda2: LASSO pentalty parameters, need to be between 0 and 1.
  # CCcoef: A 3 by 1 vector indicating weights for each pairwise canonical
  #   correlation.
  # NoTrait: Logical. Whether trait information is provided.
  # FilterByTrait: Logical. Whether only the features with highest correlation
  #   to Trait will be assigned nonzero weights. The top 80% features are reserved.
  # trace: Logical. Whether to display CCA algorithm trace.

  if(abs(Lambda1 - 0.5) > 0.5){
      stop("Invalid penalty parameter. Lambda1 needs to be between zero and
           one.")}
  if(abs(Lambda2 - 0.5) > 0.5){
      stop("Invalid penalty parameter. Lambda2 needs to be between zero and
           one.")}
  if(min(Lambda1, Lambda2) == 0){
      stop("Invalid penalty parameter. Both Lambda1 and Lambda2 has to be
           greater than 0.")
  }

  k <- ncol(Trait)
  if(NoTrait){Trait <- NULL; FilterByTrait <- TRUE; k <- 1}

  if(FilterByTrait){
    if(k > 1){
      stop("'FilterByTrait == TRUE' only allows one trait at a time.")
    }else{
      out <- CCA(X1, X2, outcome = "quantitative", y = Trait,
                 typex = "standard", typez = "standard", penaltyx = Lambda1,
                 penaltyz = Lambda2, trace = trace)
    }
  }else{
    xlist <- list(x1 = X1, x2 = X2, y = Trait)
    L1 <- max(1, sqrt(ncol(X1)) * Lambda1)
    L2 <- max(1, sqrt(ncol(X2)) * Lambda2)
    out <- myMultiCCA(xlist, penalty = c(L1, L2, sqrt(ncol(Trait))),
                      CCcoef = CCcoef, trace = trace)
    out$u <- out$ws[[1]]; out$v <- out$ws[[2]]
  }

  return(out)
}
