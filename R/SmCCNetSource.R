################################################################################
# Author: W. Jenny Shi
#
# About: This script incorporates a data feature subsampling scheme to SmCCA,  
#   SsCCA, and SCCA. Here we assume two data types (e.g. mRNA and miRNA 
#   expression levels) and one quantitative phenotype measured for the same 
#   subjects.
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



################################################################################
### Apply sparse multiple canonical correlation analysis to omics feature
### subsamples

#' Calculate the canonical correlation weights based on sparse multiple
#' canonical correlation analysis (SmCCA), sparse supervised canonical 
#' correlation analysis (SsCCA), or sparse canonical correlation analysis (SCCA).
#' 
#' Integrate two omics data type (and a quantitative phenotype), and calculate
#' the absolute canonical correlation weights for the omics features using SmCCA
#' SsCCA, or SCCA. SmCCA and SsCCA take into account a phenotype/trait. SmCCA 
#' maximizes the total (weighted or unweighted) pairwise canonical correlation 
#' weights between two omics data types and the trait. It requires the trait to 
#' be quantitative. SsCCA prioritizes omics features based on the trait, and 
#' assigns non-zero canonical weights to features that are more correlated to 
#' the trait. SCCA does not use any trait information for computing the 
#' canonical correlation weights. All of these three methods are included in 
#' this function, along with an omics feature subsampling scheme. 
#' 
#' To choose SmCCA, set \code{NoTrait = FALSE, FilterByTrait = FALSE}.  
#' To choose SsCCA, set \code{NoTrait = FALSE, FilterByTrait = TRUE}.
#' To choose SCCA, set \code{Trait = NULL, NoTrait = TRUE}.
#'
#' @param X1 An \eqn{n\times p_1} data matrix (e.g. mRNA) with \eqn{p_1} 
#'   features and \eqn{n} subjects.
#' @param X2 An \eqn{n\times p_2} data matrix (e.g. miRNA) with \eqn{p_2} 
#'   features and \eqn{n} subjects.
#' @param Trait An \eqn{n\times 1} trait data matrix for the same n subjects.
#' @param Lambda1 LASSO penalty parameter for \code{X1}. \code{Lambda1} needs
#'   to be between 0 and 1.
#' @param Lambda2 LASSO penalty parameter for \code{X2}. \code{Lambda2} needs 
#'   to be between 0 and 1.
#' @param s1 Proportion of mRNA features to be included, default at \code{s1 = 0.7}.
#'   \code{s1} needs to be between 0 and 1.
#' @param s2 Proportion of miRNA features to be included, default at \code{s1 = 0.9}.
#'   \code{s2} needs to be between 0 and 1.
#' @param NoTrait Logical, default is \code{FALSE}. Whether trait information is
#'   provided.
#' @param FilterByTrait Logical, default is \code{FALSE}. Whether only the top 
#'   (\eqn{80\%}) features with highest correlation to the trait will be assigned 
#'   nonzero weights. The choice of \eqn{80\%} is based on the PMA package.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#'   number leads to more accurate results, but at a higher cost. 
#' @param CCcoef Optional coefficients for the SmCCA pairwise canonical 
#'   correlations. If \code{CCcoef = NULL} (default), then the objective 
#'   function is the total sum of all pairwise canonical correlations. It can 
#'   also be a coefficient vector that follows the column order of 
#'   \code{combn(K, 2)}.
#' @param trace Logical. Whether to display the CCA algorithm trace.
#' @return A canonical correlation weight matrix with \eqn{p_1+p_2} rows. Each 
#'   column is the canonical correlation weights based on subsampled \code{X1}
#'   and \code{X2} features. The number of columns is \code{SubsamplingNum}.
#'
#' @examples
#' 
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#'
#' # Unweighted SmCCA
#' W1 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = FALSE,
#'   SubsamplingNum = 5, CCcoef = NULL, trace = FALSE)
#'   
#' # Weighted SmCCA
#' W2 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = FALSE,
#'   SubsamplingNum = 5, CCcoef = c(1, 5, 5), trace = FALSE)
#'   
#' # SsCCA
#' W3 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = .05, Lambda2 = 0.5,
#'   s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = TRUE,
#'   SubsamplingNum = 5, CCcoef = NULL, trace = FALSE)
#'
#' # SCCA
#' W4 <- getRobustPseudoWeights(X1, X2, Trait = NULL, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = TRUE, 
#'   SubsamplingNum = 5, CCcoef = NULL, trace = FALSE)
#'   
#'
#' @export
getRobustPseudoWeights <- function(X1, X2, Trait, Lambda1, Lambda2,
                                   s1 = 0.7, s2 = 0.9, NoTrait = FALSE, 
                                   FilterByTrait = FALSE, SubsamplingNum = 1000,
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

  beta <- pbapply::pbsapply(1:SubsamplingNum, function(x){
    # Subsample features
    samp1 <- sort(sample(1:p1, p1.sub, replace = FALSE))
    samp2 <- sort(sample(1:p2, p2.sub, replace = FALSE))

    x1.par <- scale(X1[ , samp1], center = TRUE, scale = TRUE)
    x2.par <- scale(X2[ , samp2], center = TRUE, scale = TRUE)

    out <- getCCAout(x1.par, x2.par, Trait, Lambda1, Lambda2,
                       NoTrait = NoTrait, FilterByTrait = FilterByTrait,
                       trace = trace, CCcoef = CCcoef)

    w <- rep(0, p)
    w[samp1] <- out$u
    w[samp2 + p1] <- out$v
    coeff.avg <- w

    return(coeff.avg)
  })

  return(beta)
}


################################################################################
### Aggregate pseudo canonical weights and compute the network similarity matrix
### for all omics features.

#' Compute the similarity matrix based on one or more canonical correlation 
#' weight vectors.
#' 
#' Compute the similarity matrix based on the outer products of absolute 
#' canonical correlation weights.
#'
#'
#' @param Ws A canonical correlation weight vector or matrix. If \code{Ws} is a
#'   matrix, then each column corresponds to one weight vector.
#' @param P1 Total number of features for the first omics data type. 
#' @param FeatureLabel If \code{FeatureLabel = NULL} (default), the feature 
#'   names will be \eqn{\{Gene_1, \cdots, Gene_{p_1}, Mir_1, \cdots, Mir_{p-p_1}\}}, 
#'   where \eqn{p_1 = }\code{P1}, and \eqn{p} is the total number of omics features.
#' @return A \eqn{p\times p} symmetric non-negative matrix.
#'
#' @examples
#' w <- matrix(rnorm(6), nrow = 3)
#' Ws <- apply(w, 2, function(x)return(x/sqrt(sum(x^2))))
#' abar <- getAbar(Ws, P1 = 2, FeatureLabel = NULL)
#'
#' @export
getAbar <- function(Ws, P1 = NULL, FeatureLabel = NULL){
    

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
        if(is.null(FeatureLabel)){
            if(is.null(P1)){
                stop("Need to provide FeatureLabel or the number of features 
                    for the first data type P1.")
            }else{
                p <- ncol(Abar)
                FeatureLabel <- c(paste0("Gene_", 1:P1), paste0("Mir_", 1:(p-P1)))
            }
        }
        colnames(Abar) <- rownames(Abar) <- FeatureLabel
    }

    return(Abar)
}


################################################################################
### Extract multi-omics modules.

#' Extract multi-omics modules based on the similarity matrix.
#' 
#' Apply hierarchical tree cutting to the similarity matrix and extract
#' modules that contain both omics data types.
#'
#' @param Abar A similary matrix for all features (both omics data types).
#' @param P1 Total number of features for the first omics data type.
#' @param CutHeight Height threshold for the hierarchical tree cutting. Default 
#'   is \eqn{1-0.1^{10}}.
#' @param PlotTree Logical. Whether to create a hierarchical tree plot.
#' @return A list of multi-omics modules.
#'
#' @examples
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' abar <- getAbar(w, P1 = 2, FeatureLabel = NULL)
#' modules <- getMultiOmicsModules(abar, P1 = 2, CutHeight = 0.5)
#'
#' @export
getMultiOmicsModules <- function(Abar, P1, CutHeight = 1-.1^10, PlotTree = TRUE){

    hc <- stats::hclust(stats::as.dist(1 - Abar))
    if(PlotTree){graphics::plot(hc)}
    cut.merge <- hc$merge[hc$height < CutHeight, ]
    lower.leaves <- sort(-cut.merge[cut.merge<0])

    grpID <- stats::cutree(hc, h = CutHeight)
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

#' Plot multi-omics module networks.
#' 
#' Plot multi-omics modules based on similarity matrix derived from pseudo
#' canonical weights and pairwise feature correlations.
#'
#' @param Abar A \eqn{p\times p} similary matrix for both omics data types
#'   based on pseudo canonical correlation weights. \eqn{p} is the number of 
#'   total features for the two omics data types. All entries are non-negative.
#' @param CorrMatrix A \eqn{p\times p} correlation matrix that provides sign
#'   information for the network.
#' @param multiOmicsModule A list of multi-omics modules.
#' @param ModuleIdx Index for the module to be plotted. It can not exceed the
#'   length of \code{multiOmicsModule}.
#' @param P1 Total number of features for the first omics data type.
#' @param EdgeCut A numerical value between 0 and 1, indicating an edge
#'   threshold for the network. Any features (network nodes) without any edge
#'   strength that passes the threshold are excluded from the figure. If 
#'   \code{EdgeCut = 0} (default), then the full module network will be created. 
#' @param FeatureLabel A \eqn{1\times p} vector indicating feature names. If
#'   \code{FeatureLabel = NULL} (default), the feature names will be 
#'   \eqn{\{Gene_1, \cdots, Gene_{p_1}, Mir_1, \cdots, Mir_{p-p_1}\}}, where 
#'   \eqn{p_1 = }\code{P1}.
#' @param AddCorrSign Logical. Whether to add a positive or negative sign to
#'   each network edge based on pairwise feature correlations.
#' @param SaveFile A pdf file name for the figure output. 
#'   If \code{SaveFile = NULL} (default), the figure will not be saved.
#' @param ShowType1Label Logical. Whether to label the network nodes for the
#'   first omics data type.
#' @param ShowType2Label Logical. Whether to label the network nodes for the
#'   second omics data type.
#' @param PlotTitle A title for the figure. Default is without any title.
#' @param NetLayout Graphical layout for the network. Possible options are
#'   \code{circle} for circle layout, \code{sphere} for 3D sphere, \code{fr} for
#'   Fruchterman-Reinhold, and \code{lgl} for the LGL algorithm. Refer to igraph 
#'   manual for more details on the layout options.
#' @param ShowNodes Logical. Whether to show network nodes.
#' @param VertexLabelCex Scaling factor for the vertex labels.
#' @param VertexSize Size of the vertices.
#' @return A multi-omics network figure.
#'
#' @examples
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' abar <- getAbar(w, P1 = 2, FeatureLabel = NULL)
#' modules <- getMultiOmicsModules(abar, P1 = 2, CutHeight = 0.5)
#' x <- cbind(X1[ ,1:2], X2[ , 1:3])
#' corr <- cor(x)
#'
#' plotMultiOmicsNetwork(abar, corr, modules, ModuleIdx = 1, P1 = 2)
#'
#' @export
plotMultiOmicsNetwork <- function(Abar, CorrMatrix, multiOmicsModule,
                               ModuleIdx, P1, EdgeCut = 0, FeatureLabel = NULL,
                               AddCorrSign = TRUE, SaveFile = NULL,
                               ShowType1Label = TRUE, ShowType2Label = TRUE,
                               PlotTitle = "", NetLayout = "lgl",
                               ShowNodes = TRUE,
                               VertexLabelCex = 1, VertexSize = 1){

    p <- ncol(Abar)
    if(is.null(FeatureLabel)){
        FeatureLabel <- c(paste0("Gene_", 1:P1), paste0("Mir_", 1:(p-P1)))
        }
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
        net <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
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
        ecol <- rep("gray80", igraph::ecount(net))
        ew <- abs(igraph::edge.attributes(net)$weight) * 5
        ecol[which(igraph::edge.attributes(net)$weight < 0)] <- "red"

        # Define network layout.
        if(NetLayout == "circle"){
            l <- igraph::layout_in_circle(net)
        }else if(NetLayout == "sphere"){
            l <- igraph::layout_on_sphere(net)
        }else if(NetLayout == "fr"){
            l <- igraph::layout_with_fr(net)
        }else if(NetLayout == "lgl"){
            l <- igraph::layout_with_lgl(net)
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
            graphics::title(PlotTitle)
            grDevices::dev.off()
        }else{
            graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
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
  if(NoTrait | is.null(k)){
      out <- PMA::CCA(X1, X2, typex = "standard", typez = "standard", 
                      penaltyx = Lambda1, penaltyz = Lambda2, K = 1, 
                      trace = trace)
  }else{
      if(FilterByTrait){
          if(k > 1){
              stop("'FilterByTrait == TRUE' only allows one trait at a time.")
          }else{
              out <- PMA::CCA(X1, X2, outcome = "quantitative", y = Trait,
                              typex = "standard", typez = "standard", 
                              penaltyx = Lambda1, penaltyz = Lambda2, K = 1, 
                              trace = trace)
          }
      }else{
          xlist <- list(x1 = X1, x2 = X2, y = scale(Trait))
          L1 <- max(1, sqrt(ncol(X1)) * Lambda1)
          L2 <- max(1, sqrt(ncol(X2)) * Lambda2)
          out <- myMultiCCA(xlist, penalty = c(L1, L2, sqrt(ncol(Trait))),
                            CCcoef = CCcoef, trace = trace)
          out$u <- out$ws[[1]]; out$v <- out$ws[[2]]
      }
  }

  return(out)
}
