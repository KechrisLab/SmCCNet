#' Calculate similarity matrix based on canonical weights.
#' 
#' Compute the similarity matrix based on the outer products of absolute
#' canonical correlation weights, can be used for both single and multi-omics setting.
#' 
#' 
#' @param Ws A canonical correlation weight vector or matrix. If \code{Ws} is a
#' matrix, then each column corresponds to one weight vector.
#' @param FeatureLabel A vector of feature labels for each feature in the adjacency matrix
#' @return A \eqn{p\times p} symmetric non-negative matrix.
#' @examples
#' 
#' w <- matrix(rnorm(6), nrow = 3)
#' Ws <- apply(w, 2, function(x)return(x/sqrt(sum(x^2))))
#' abar <- getAbar(Ws,  FeatureLabel = c('omics1', 'omics2', 'omics3'))
#' 
#' @export
getAbar <- function(Ws, FeatureLabel = NULL){
  # require(Matrix)
  if(is.null(dim(Ws))){
    Abar <- Matrix::Matrix(abs(Ws) %o% abs(Ws), sparse = TRUE)
  }else{
    b <- nrow(Ws)
    Abar <- matrix(0, nrow = b, ncol = b)
    for(ind in seq_len(ncol(Ws))){
      # w <- methods::as(abs(Ws[ , ind]), "sparseVector")
      w <- abs(Ws[ , ind])
      A <- Matrix::tcrossprod(w, w)
      Abar <- Abar + A
    }
  }
  # set the diagonal of the adjacency matrix to 0
  diag(Abar) <- 0
  Abar <- Abar/max(Abar)
  if(is.null(FeatureLabel)){
      stop("Need to provide FeatureLabel.")
  }
  # set feature names in adjacency matrix
  colnames(Abar) <- rownames(Abar) <- FeatureLabel
  # return the adjacency matrix
  return(Abar)
}






#' preprocess a omics dataset before running omics SmCCNet
#' 
#' Data preprocess pipeline to: (1) filter by coefficient of variation (cv), (2) center or scale data
#' and (3) adjust for clinical covariates.
#' 
#' 
#' @param X dataframe with the size of \eqn{n} by \eqn{p}, where \eqn{n} is the sample size and \eqn{p} is the feature size.
#' @param covariates dataframe with covariates to be adjusted for.
#' @param is_cv Whether to use coefficient of variation filter (small cv filter out). 
#' @param cv_quantile CV filtering quantile.
#' @param center Whether to center the dataset X.
#' @param scale Whether to scale the dataset X.
#' @return Processed omics data with the size of nxp. 
#' @examples
#' 
#' X1 <- as.data.frame(matrix(rnorm(600, 0, 1), nrow = 60))
#' covar <- as.data.frame(matrix(rnorm(120, 0, 1), nrow = 60))
#' processed_data <- dataPreprocess(X = X1, covariates = covar, is_cv = TRUE, 
#' cv_quantile = 0.5, center = TRUE, scale = TRUE)
#' 
#' @export
#' 
dataPreprocess <- function(X, covariates = NULL, 
                            is_cv = FALSE, cv_quantile = 0, center = TRUE, scale = TRUE){
  # placeholder for X
  data <- X
 
  
  if (is_cv == TRUE)
  {
    cv_value = apply(data, 2, EnvStats::cv)
    cv_value = abs(cv_value)
    if (is.null(cv_quantile))
      stop("cv filtering quantile must be provided!")
    # filter data column based on cv quantile
    data <- data[, which(cv_value > stats::quantile(cv_value, cv_quantile))]
  }
  
  # center and scale data
  if (center == TRUE)
    data <- as.data.frame(scale(data, center = TRUE, scale = FALSE))
  if (scale == TRUE)
    data <- as.data.frame(scale(data, center = FALSE, scale = TRUE))
  
  if(!is.null(covariates))
  {
    if(!is.data.frame(covariates))
      stop('Covariate data should be in dataframe format!')
    if(ncol(covariates) == 0)
      stop('Covariate dataframe must have at least 1 column!')
    # regress out covariates
    
    
    #map_data <- data %>%  
     # purrr::map(function(x){
      #  formula <- stats::as.formula(paste("x ~ ", paste(paste0('covariates$',colnames(covariates)), collapse= "+")))
       # return(stats::lm(formula, data = data))
      #}) %>%
      #purrr::map(stats::resid)
    
    map_data <- purrr::map(data, function(x){
        formula <- stats::as.formula(paste("x ~ ", paste(paste0('covariates$',colnames(covariates)), collapse= "+")))
        return(stats::lm(formula, data = data))
      })
    map_data <- (purrr::map(map_data, stats::resid))
    
    
    
    adjusted_data <- purrr::map_df(map_data, ~as.data.frame(t(.)))
    adjusted_data <- t(adjusted_data)
    colnames(adjusted_data) <- colnames(data)
    data <- adjusted_data
    
  }
  return(data)
  
  
}






#' Internal functions called by getRobustPseudoWeights_single.
#' 
#' 
#' @param X1 data.
#' @param Trait phenotype.
#' @param Lambda1 penalty term
#' @param trace Whether to display CCA algorithm trace.
#' 
#' @keywords internal


getCCAout_single <- function(X1, Trait, Lambda1, trace = FALSE){
  # Compute CCA weights.
  #
  # Protein: An n by p1 mRNA expression matrix.
  # Trait: An n by k trait data for the same samples (k >=1).
  # Lambda1: LASSO pentalty parameters, need to be between 0 and 1.
  
  
  if(abs(Lambda1 - 0.5) > 0.5){
    stop("Invalid penalty parameter. Lambda1 needs to be between zero and
           one.")}
  
  
  k <- ncol(Trait)    
  
  
  xlist <- list(x1 = X1, y = scale(Trait))
  L1 <- max(1, sqrt(ncol(X1)) * Lambda1)
  out <- myMultiCCA(xlist, penalty = c(L1, sqrt(ncol(Trait))), trace = trace)
  out$u <- out$ws[[1]]
  
  return(out)
}







#' Single-omics SmCCA with Binary Phenotype
#' 
#' Compute aggregated (SmCCA) canonical weights for single omics data with quantitative phenotype (subampling enabled). 
#' 
#' 
#' @param X1 An \eqn{n\times p_1} data matrix (e.g. mRNA) with \eqn{p_1}
#' features and \eqn{n} subjects.
#' @param Trait An \eqn{n\times 1} trait (phenotype) data matrix for the same \eqn{n} subjects.
#' @param Lambda1 LASSO penalty parameter for \code{X1}. \code{Lambda1} needs
#' to be between 0 and 1.
#' @param s1 Proportion of mRNA features to be included, default at \code{s1 =
#' 0.7}. \code{s1} needs to be between 0 and 1, default is set to 0.7.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#' number leads to more accurate results, but at a higher computational cost.
#' @param K Number of hidden components for PLSDA, default is set to 3.
#' @return A partial least squared weight matrix with \eqn{p_1} rows. Each
#' column is the canonical correlation weights based on subsampled \code{X1}
#' features. The number of columns is \code{SubsamplingNum}.
#' @examples
#' 
#' 
#' X <- matrix(rnorm(600,0,1), nrow = 60)
#' Y <- rbinom(60,1,0.5)
#' Ws <- getRobustWeightsSingleBinary(X1 = X, Trait = as.matrix(Y), Lambda1 = 0.8, 
#' 0.7, SubsamplingNum = 10)
#' 
#' @export
#' 

getRobustWeightsSingleBinary <- function(X1, Trait, Lambda1, 
                                         s1 = 0.7, SubsamplingNum = 1000, 
                                         K = 3
){
  
  # This is a source function when binary outcome variable is used, it takes in the covariates data
  # phenotype data, penalty as well as the subsampling rate and the number of subsamples
  
  p1 <- ncol(X1); 
  p1.sub <- ceiling(s1 * p1)
  # Xdata <- X1
  
  beta <- sapply(seq_len(SubsamplingNum), function(x){
    # Subsample features
    samp1 <- sort(sample(seq_len(p1), p1.sub, replace = FALSE))
    
    x1.par <- scale(X1[ , samp1], center = TRUE, scale = TRUE)
    
    out <- spls::splsda(x = x1.par, y = Trait, K = K, eta = Lambda1, kappa=0.5,
                        classifier = 'logistic', scale.x=FALSE)
    u <- rep(0, p1.sub)
    w <- rep(0, p1)
    
    # fit logistic regression model with latent factors
    logit.data <- as.data.frame(cbind(Trait, out[["T"]]))
    colnames(logit.data)[1] <-'Y'
    glm.result <- stats::glm(Y~.-1, family = 'binomial', data = logit.data)
    # average latent weight based on logistic regression coefficients
    u[out[["spls.fit"]][["A"]]] <- abs(out[["W"]]) %*% abs(glm.result$coefficients)
    # normalize
    u[out[["spls.fit"]][["A"]]] <- u[out[["spls.fit"]][["A"]]]/pracma::Norm(u[out[["spls.fit"]][["A"]]])
    w[samp1] <- u
    coeff.avg <- w
    
    return(coeff.avg)
  })
  
  return(beta)
}




#' Single-omics SmCCA with Quantitative Phenotype
#' 
#' Compute aggregated (SmCCA) canonical weights for single omics data with quantitative phenotype (subampling enabled). 
#' 
#' 
#' @param X1 An \eqn{n\times p_1} data matrix (e.g. mRNA) with \eqn{p_1}
#' features and \eqn{n} subjects.
#' @param Trait An \eqn{n\times 1} trait (phenotype) data matrix for the same \eqn{n} subjects.
#' @param Lambda1 LASSO penalty parameter for \code{X1}. \code{Lambda1} needs
#' to be between 0 and 1.
#' @param s1 Proportion of features in \code{X1} to be included, default at \code{s1 =
#' 0.7}. \code{s1} needs to be between 0 and 1, default is set to 0.7.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#' number leads to more accurate results, but at a higher computational cost.
#' @param trace Whether to display the CCA algorithm trace, default is set to FALSE.
#' @return A canonical correlation weight matrix with \eqn{p_1} rows. Each
#' column is the canonical correlation weights based on subsampled \code{X1}
#' features. The number of columns is \code{SubsamplingNum}.
#' @examples
#' 
#' 
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#' 
#' # Single Omics SmCCA
#' W1 <- getRobustWeightsSingle(X1, Trait = Y, Lambda1 = 0.05,
#'   s1 = 0.7, 
#'   SubsamplingNum = 5, trace = FALSE)
#'   
#'   
#' @export
#' 
getRobustWeightsSingle <- function(X1, Trait, Lambda1, 
                                     s1 = 0.7, SubsamplingNum = 1000, trace = FALSE){
  
  
  
  p1 <- ncol(X1); p <- p1 
  p1.sub <- ceiling(s1 * p1)
  X <- X1
  
  beta <- pbapply::pbsapply(seq_len(SubsamplingNum), function(x){
    # Subsample features
    samp1 <- sort(sample(seq_len(p1), p1.sub, replace = FALSE))
    
    
    x1.par <- scale(X1[ , samp1], center = TRUE, scale = TRUE)
    
    
    out <- getCCAout_single(x1.par, Trait, Lambda1,
                            trace = trace)
    w <- rep(0, p)
    w[samp1] <- out$u
    coeff.avg <- w
    
    return(coeff.avg)
  })
  
  return(beta)
}










#' Extract Omics Modules based on Similarity Matrix.
#' 
#' Apply hierarchical tree cutting to the similarity matrix and extract multi/single-omics network modules.
#' 
#' 
#' 
#' @param Abar A similary matrix for all features (all omics data types).
#' @param CutHeight Height threshold for the hierarchical tree cutting. Default
#' is \eqn{1-0.1^{10}}.
#' @param PlotTree Logical. Whether to create a hierarchical tree plot, default is set to \code{TRUE}.
#' @return A list of multi/single-omics modules.
#' @examples
#' 
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' feature_name <- paste0('feature_', 1:5)
#' abar <- getAbar(w, FeatureLabel = feature_name)
#' modules <- getOmicsModules(abar, CutHeight = 0.5)
#' @export
#' 
getOmicsModules <- function(Abar, CutHeight = 1-.1^10, PlotTree = TRUE){
  
  hc <- stats::hclust(stats::as.dist(1 - Abar))
  if(PlotTree){graphics::plot(hc)}
  cut.merge <- hc$merge[hc$height < CutHeight, ]
  lower.leaves <- sort(-cut.merge[cut.merge<0])
  
  grpID <- stats::cutree(hc, h = CutHeight)
  id <- grpID[lower.leaves]
  M <- lapply(seq_len(length(unique(id))), function(x){
    M.x <- lower.leaves[which(id == unique(id)[x])]
    return(M.x)
  })
  
  multiOmicsModule <- lapply(M, function(s){
    s.min <- min(s)
    s.max <- max(s)
    return(s)
  })
  
  if(length(multiOmicsModule) > 1){
    nullSet <- which(vapply(multiOmicsModule, is.null, logical(1)))
    if(length(nullSet) > 0){
      multiOmicsModule <- multiOmicsModule[-nullSet]
    }
  }
  
  return(multiOmicsModule)
}





#' Prunes Subnetwork and Return Final Pruned Subnetwork Module
#' 
#' Prunes subnetworks with network pruning algorithm (see multi-omics vignette for detail), and save the final pruned subnetwork to the user-defined directory.
#' The final subnetwork is an .Rdata file with a name 'size_m_net_ind.Rdata', where \eqn{m} is the final pruned network size, and ind is the index of the subnetwork module after hierarchical clustering.
#' 
#' 
#' @param Abar Adjacency matrix of subnetwork with size \eqn{m^{*}} by \eqn{m^{*}} after hierarchical clustering.
#' @param CorrMatrix The correlation matrix of features in \code{Abar}, it should be
#' \eqn{m^{*}} by \eqn{m^{*}} as well.
#' @param data The omics data for the subnetwork.
#' @param Pheno The trait (phenotype) data used for network pruning.
#' @param ModuleIdx The index of the network module that summarization score is
#' intended to be stored, this is used for naming the subnetwork file in user-defined directory.
#' @param min_mod_size The minimally possible subnetwork size for the pruned network module,
#' should be an integer from 1 to the largest possible size of the subnetwork, default is set to 10.
#' @param max_mod_size the maximally possible subnetwork size for the pruned network module,
#' should be an integer from 1 to the largest possible size of the subnetwork, and it needs to be greater than the value specified in \code{min_mod_size}.
#' @param type A vector with length equal to total number of features in the adjacency matrix
#' indicating the type of data for each feature. For instance, for a subnetwork with 2 genes and a protein, the \code{type} argument should be set to \code{c('gene', 'gene', 'protein')}, see multi-omics vignette for more information.
#' @param damping damping parameter for the PageRank algorithm, default is set to 0.9, see \code{igraph} package for more detail.
#' @param method Selection between NetSHy' and 'PCA', specifying the network summarization method used for network pruning, default is set to NetSHy.
#' @param saving_dir User-defined directory to store pruned subnetwork.
#' @return A file stored in the user-defined directory, which contains the
#' following: (1) correlation_sub: correlation matrix for the subnetwork.
#' (2) M: adjacency matrix for the subnetwork.
#' (3) omics_corelation_data: individual molecular feature correlation with phenotype.
#' (4) pc_correlation: first 3 PCs correlation with phenotype.
#' (5) pc_loading: principal component loadings.
#' (6) pca_x1_score: principal component score and phenotype data.
#' (7) mod_size: number of molecular features in the subnetwork.
#' (8) sub_type: type of feature for each molecular features.
#' @examples
#' library(SmCCNet)
#' set.seed(123)
#' w <- rnorm(20)
#' w <- w/sqrt(sum(w^2))
#' labels <- paste0('feature_', 1:20)
#' abar <- getAbar(w, FeatureLabel = labels)
#' modules <- getOmicsModules(abar, CutHeight = 0.1)
#' x <- X1[ ,seq_len(20)]
#' corr <- stats::cor(x)
#' # display only example
#' # networkPruning(abar, corr, data = x, Pheno = Y,
#' # ModuleIdx = 1,  min_mod_size = 3, max_mod_size = 10, method = 'NetSHy', saving_dir = 
#' # )
#' 
#' @export
#' 

networkPruning <- function(Abar, CorrMatrix, data,
                                             Pheno, type, ModuleIdx,
                                             min_mod_size = 10, max_mod_size, damping = 0.9,
                                             method = 'NetSHy', saving_dir){
  
  
  
  P1 <- p <- ncol(Abar)
  # Trim module by PPR
  net_ppr <- igraph::graph_from_adjacency_matrix(Abar, weighted = TRUE,
                                                 diag = FALSE, mode = "undirected")
  igraph::set_vertex_attr(net_ppr, "type", index = igraph::V(net_ppr), as.factor(type))
  # All parameter setups are based on the recommendation
  ranking <- igraph::page_rank(net_ppr, directed = FALSE, damping = damping, 
                               options = list(niter = 10^5,eps = 1e-06))
  # Obtain ranked protein names
  rank_names <- names(sort(ranking$vector, decreasing = TRUE))
  rank_value <- sort(ranking$vector, decreasing = TRUE)
  names(rank_value) <- rank_names
  # Choose to include only the top proteins as we desired
  summary_correlation <- c()
  if (max_mod_size > nrow(Abar))
    max_mod_size <- nrow(Abar)
  
  for (i in min_mod_size : max_mod_size)
  { 
    # print(paste0('iteration:', i))
    newM.node <- which(colnames(Abar) %in% rank_names[1:i])
    M <- Abar[newM.node, newM.node]
    
    ###### section for principal component analysis
    X1_PC <- data[,which(colnames(data) %in% rank_names[1:i])]
    #print(X1_PC)
    # calculate individual protein correlation
    protein_correlation <- stats::cor(X1_PC, Pheno)
    protein_correlation_data <- data.frame(name = colnames(X1_PC), correlation = protein_correlation)
    # PCA function
    if (method == 'PCA')
    {
      # run pca
      pca_x1_summary <- stats::prcomp(X1_PC)
      # Obtain summary
      summary_result <- summary(pca_x1_summary)
      # Extract the first three PC scores
      pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary$x[,1], 
                               pc2 = pca_x1_summary$x[,2],
                               pc3 = pca_x1_summary$x[,3],
                               y = Pheno)
    }else if(method == 'NetSHy'){
      pca_x1_summary <- summarizeNetSHy(X1_PC, M, 
                                     npc = 3)
      # Extract the first three PC scores
      pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary[[1]][,1], 
                               pc2 = pca_x1_summary[[1]][,2],
                               pc3 = pca_x1_summary[[1]][,3],
                               y = Pheno)
      
    }else{
      stop('Either method of PCA or NetSHy should be provided.')
    } 
    
    
    pc_correlation <- stats::cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4])
    summary_correlation <- c(summary_correlation, pc_correlation[1])
    #cor_index <- which.max(abs(cor(pca_x1_pc1, Y)))
    #if (i == min_mod_size)
     # score_frame <- data.frame(pca_x1_pc1$pc1)
    #else
      #score_frame = cbind(score_frame, pca_x1_pc1$pc1)
    
    if (i == min_mod_size)
    {
      temp_cor <- abs(stats::cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4]))
      score_frame <- data.frame(pca_x1_pc1[,which.max(temp_cor)])
    } else
    {  
      temp_cor <- abs(stats::cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4]))
      score_frame = cbind(score_frame, pca_x1_pc1[,which.max(temp_cor)])
    }  
  }
  
  corr_pheno <- abs(stats::cor(score_frame, Pheno))
  
  # selecting the optimal network size
  # candidate_size_1 <- min(which(corr_pheno > (0.9 * max(corr_pheno))))
  candidate_size_1 <- min(which.max(corr_pheno))
  cormat <-  abs(round(x = stats::cor(score_frame[,candidate_size_1],score_frame[,candidate_size_1:(max_mod_size - min_mod_size + 1)]), digits = 2))
  candidate_size_2 <- max(which(cormat > 0.8))
  #print(candidate_size_2)
  candidate_size_3 <- max(which(corr_pheno[candidate_size_1: (candidate_size_1 + candidate_size_2 - 1)] > (0.9 * max(corr_pheno))))
  mod_size <- candidate_size_3 + candidate_size_1 - 2 + min_mod_size
  #print(c(candidate_size_1, candidate_size_2, candidate_size_3))
  
  # finally calculate the optimal network size 
  #  mod_size <- selection(cormat[,1], summary_correlation,cor_cut = 0.8, default_size = min_mod_size, network_preference = network_preference)
  # print(mod_size)
  # obtain optimal result
  newM.node <- which(colnames(Abar) %in% rank_names[1:mod_size])
  sub_type <- type[newM.node]
  M <- Abar[newM.node, newM.node]
  
  # print(corr_pheno)
  
  
  ###### section for principal component analysis
  X1_PC <- data[,which(colnames(data) %in% rank_names[1:mod_size])]
  # calculate individual protein correlation
  protein_correlation <- stats::cor(X1_PC, Pheno)
  protein_correlation_test <- rep(0, ncol(X1_PC))
  for (i in 1:ncol(X1_PC))
  {
    protein_correlation_test[i] <- stats::cor.test(X1_PC[,i], Pheno)$p.value
  }
  omics_correlation_data <- data.frame(name = colnames(X1_PC), correlation = protein_correlation, p = protein_correlation_test)
  if (method == 'PCA')
  {
    # run pca
    pca_x1_summary <- stats::prcomp(X1_PC)
    # Obtain summary
    summary_result <- summary(pca_x1_summary)
    # Extract the first three PC scores
    pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary$x[,1], 
                             pc2 = pca_x1_summary$x[,2],
                             pc3 = pca_x1_summary$x[,3],
                             y = Pheno)
    pc_loading <- summary_result[["rotation"]]
  }else if(method == 'NetSHy'){
    pca_x1_summary <- summarizeNetSHy(X1_PC, M, 
                                   npc = 3)
    # Extract the first three PC scores
    pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary[[1]][,1], 
                             pc2 = pca_x1_summary[[1]][,2],
                             pc3 = pca_x1_summary[[1]][,3],
                             y = Pheno)
    pc_loading <- pca_x1_summary[[3]]
    
  } 
  pc_correlation <- stats::cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4])
  summary_correlation <- c(summary_correlation, pc_correlation[1])
  
  correlation_sub <- CorrMatrix[newM.node, newM.node]
  
  
  
  
  
  # Save all the cc results into a same data file
  save(M = M, pca_score = pca_x1_pc1,pc_loading, rank_value = rank_value, pc_correlation, omics_correlation_data,
       mod_size, sub_type, summary_correlation, correlation_sub,
       file = paste0(saving_dir, "/size_", mod_size,"_net_",ModuleIdx,".Rdata"))
  max_cor_index <- which.max(abs(pc_correlation))
  cat(paste0('The final network size is: ', nrow(M), ' with maximum PC correlation w.r.t. phenotype to be: ', round(pc_correlation[max_cor_index], 3)))
  #return(c(nrow(M),pc_correlation))
}








#' Aggregate and Save Cross-validation Result for Single-omics Analysis
#' 
#' Saves cross-validation results in a table with the
#' user-defined directory and outputs penalty term with the highest testing canonical correlation,
#' lowest prediction error, and lowest scaled prediction error.
#' 
#' 
#' @param CVDir A directory where the result is stored.
#' @param SCCAmethod The canonical correlation analysis method that is used in
#' the model, used to name cross-validation table file, default is set to 'SmCCA'.
#' @param K number of folds for cross-validation.
#' @param NumSubsamp Number of subsampling used.
#' @return A vector of length 3 with indices of the penalty term that (1) maximize the testing canonical correlation,
#' (2) minimize the prediction error and (3) minimize the scaled prediction error.
#' 
#' @export

aggregateCVSingle <- function(CVDir, SCCAmethod = 'SmCCA', K = 5, NumSubsamp = 500){
  plotD <- paste0(CVDir, "Figures/")
  saveD <- paste0(CVDir, "Results/")
  dir.create(plotD); dir.create(saveD)
  
  testCC <- predError <- NULL
  for(j in 1:K){
    resultT <- paste0(CVDir, "CV_", j, "/", SCCAmethod,
                      "/SCCA_", NumSubsamp, "_allDeltaCor.csv")
    dCorT <- utils::read.csv(resultT)[ , -1]
    testCC <- cbind(testCC, abs(dCorT[ , 3]))
    predError <- cbind(predError, dCorT[ , 4])
  }
  
  S1 <- rowMeans(testCC)
  S2 <- rowMeans(predError)
  S3 <- abs(S2/abs(S1))
  T12 <- dCorT[ , -2]; T12[ , 2] <- S1; T12[ , 3] <- S2
  utils::write.csv(T12, file = paste0(saveD, SCCAmethod, "CVmeanDeltaCors.csv"))
  
  print(paste0("testCC choice: ", which(S1 == max(S1))))
  print(paste0("CC Pred. Error choice: ", which(S2 == min(S2))))
  choices <- c(which(S1 == max(S1)), which(S2 == min(S2)), which(S3 == min(S3)))
  names(choices) <- c('test.cc', 'pred.error', 'scaled.pred.error')
  return(c(which(S1 == max(S1)), which(S2 == min(S2)), which(S3 == min(S3))))
}





#'@title NetSHy Summarization Score 
#'
#'@description Implement NetSHy network summarization via a hybrid approach (Vu et al.,) to 
#'summarize network by considering the network topology with Laplacian matrix.
#'
#'@param X An \eqn{n\times m} data matrix  with \eqn{m}
#' features and \eqn{n} subjects.
#'@param A Corresponding adjacency matrix of size \eqn{p} by \eqn{p}.
#'@param npc Number of principal components used to summarize the network, default is set to 1.
#'
#'
#'@return A list consists of (1) subject-level network summarization score,
#'(2) principal component importance information: standard deviation, percent of variance explained, and cumulative proportion of variance explained, and (3) principal component feature-level loadings.
#'
#'@examples 
#'# simulate omics data
#'OmicsData <- matrix(rnorm(200,0,1), nrow = 10, ncol = 20)
#'# simulate omics adjacency matrix
#'set.seed(123)
#'w <- rnorm(20)
#'w <- w/sqrt(sum(w^2))
#'featurelabel <- paste0('omics',1:20)
#'abar <- getAbar(w, FeatureLabel = featurelabel)
#'# extract NetSHy summarization score
#'netshy_score <- summarizeNetSHy(OmicsData, abar)
#'
#'@references 
#'Vu, Thao, Elizabeth M. Litkowski, Weixuan Liu, Katherine A. Pratte, Leslie Lange, Russell P. Bowler, Farnoush Banaei-Kashani, and Katerina J. Kechris. "NetSHy: network summarization via a hybrid approach leveraging topological properties." Bioinformatics 39, no. 1 (2023): btac818.
#'
#'@export

summarizeNetSHy = function(X, A, npc=1){
  #X: data matrix (n, p)
  #A: corresponding adjacency matrix
  #pc_id: PC index
  g = igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  # laplacian
  L2 <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE,diag = FALSE) 
  L2 <- as.matrix(igraph::graph.laplacian(L2, normalized = F))
  # TOM
  # TOM_A = WGCNA::TOMsimilarity(as.matrix(A), verbose = FALSE)
  
  alpha = igraph::graph.density(g)
  X= scale(X, center = TRUE, scale = TRUE)
  # weighted approach
  
  temp = (X %*% L2)
  temp = summary(stats::prcomp(temp))
  h_score = temp$x[,1:npc] 
  importance = temp$importance[,1:npc]
  loading = temp$rotation[,1:npc]
  return(list(h_score, importance, loading))
  
}





#' Run Sparse multiple Canonical Correlation Analysis and Obtain Canonical Weights (with Subsampling)
#' 
#' SmCCNet algorithm with multi-omics data and binary phenotype. This is a stepwise approach 
#' (1) use SmCCA to identify relationship between omics (exlude phenotype), (2) within highly connected omics features
#' selected in step 1, identify relationship between these selected omics features and phenotype of interest with 
#' sparse PLS. First, it computes PLSDA by assuming outcome is continuous to extract multiple latent factors, then uses latent factors to fit logistic regression, and weight latent factor by
#' regression parameters. Refer to multi-omics vignette for more detail. 
#' 
#' @param X A list of omics data each with n subjects.
#' @param Y A vector of binary variable, user needs to set the level of this 
#' variable to 0 and 1.
#' @param Between_Discriminate_Ratio A vector with length 2 specifying the relative importance
#' of between-omics relationship and omics-phenotype relationship. For instance a ratio of 1:1 (c(1,1) in the argument)
#' means between-omics relationship and omics-phenotype relationship contribute equally to the canonical weights extraction.
#' @param SubsamplingPercent A vector with length equal to the number of omics data (\code{X}), specifying the 
#' percentage of omics feature being subsampled at each subsampling iteration.
#' @param CCcoef A vector of scaling factors only for between-omics relationship (exclude omics-phenotype). This 
#' coefficient vector follows the column order of \code{combn(T, 2)} when there are \code{T} omics data.
#' @param LambdaBetween A vector of sparsity penalty value for each omics data to run the between-omics SmCCA, each 
#' penalty term should be within the range of 0 and 1.
#' @param LambdaPheno A penalty term when running the sparse PLS with phenotype, penalty term should be within the range of 0 and 1.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#' number leads to more accurate results, but at a higher computational cost, default is set to 1000.
#' @param ncomp_pls Number of latent components for PLS, default set to 3.
#' @param EvalClassifier If \code{TRUE}, the algorithm is at the phase of evaluating classification performance, and the latent factors from SPLSDA will be returned; if FALSE, the algorithm is at the phase of constructing multi-omics network, canonical weight will be returned. 
#' Default is set to \code{FALSE}. 
#' @param testData A list of testing omics data matrix, should have the exact same order as data list X, only used when EvalClassifier is set to \code{TRUE} for performing cross-validation, refer to multi-omics vignette for detail.
#' @return If \code{EvalClassifier} is set to \code{FALSE}, a canonical correlation weight matrix is returned with combined omics data. Each
#' column is the canonical correlation weights based on subsampled X features. The number of columns is \code{SubsamplingNum}. If \code{EvalClassifier} is set to \code{TRUE}, then latent factors from training and testing data will be returned for classifier evaluation. 
#' @examples
#' 
#' 
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#' X1 <- matrix(rnorm(600,0,1), nrow = 60)
#' X2 <- matrix(rnorm(600,0,1), nrow = 60)
#' Y_binary <- rbinom(60,1,0.5)
#' 
#' Ws <- getRobustWeightsMultiBinary(list(X1,X2), Y_binary, 
#'       SubsamplingPercent = c(0.8,0.8), CCcoef = NULL,
#'       LambdaBetween = c(0.5,0.5), LambdaPheno = 0.1, SubsamplingNum = 10)
#'   
#' @export
#' 

getRobustWeightsMultiBinary <- function(X, Y, Between_Discriminate_Ratio = c(1,1),
                                              SubsamplingPercent = NULL, CCcoef = NULL, LambdaBetween, LambdaPheno = NULL,
                                              SubsamplingNum = 1000, ncomp_pls = 3, EvalClassifier = FALSE, testData = NULL)
{
  eta <- LambdaPheno
  # run between-omics SmCCA
  BetweenOmicsWeight <- getRobustWeightsMulti(X, Trait = NULL, NoTrait = TRUE, CCcoef = CCcoef,
                                         Lambda = LambdaBetween, s = SubsamplingPercent, SubsamplingNum = SubsamplingNum)
  # column bind all data
  X_all <- rlist::list.cbind(X)
  # Feature type index
  type_index <- unlist(purrr::map(1:length(X), function(h){
    rep(h, ncol(X[[h]]))
  }))
  if (EvalClassifier == FALSE)
  {  
  # create empty matrix to store the omics-phenotype weight
  OmicsPhenoWeight <- matrix(0, nrow = nrow(BetweenOmicsWeight), ncol = ncol(BetweenOmicsWeight))
  # run omics-phenotype PLS  
  for (iii in 1:ncol(BetweenOmicsWeight))
  {
    # subset selected features
    X_subset <- X_all[ ,which(BetweenOmicsWeight[,iii] != 0)]
    # run omics-phenotype SmCCA based on selected molecular features
      Ws_pheno <- getRobustWeightsSingleBinary(X1 = X_subset, Trait = matrix(Y, ncol = 1), Lambda1 = as.numeric(eta), 
                                               s1 = 1, SubsamplingNum = 1, K = ncomp_pls)
      # store omics-phenotype canonical weight
      OmicsPhenoWeight[which(BetweenOmicsWeight[,iii] != 0),iii] <- as.numeric(Ws_pheno)
      
   
    

    # normalize each data type 
    for (j in 1:length(X))
    {
      OmicsPhenoWeight[which(type_index == j),iii] <- OmicsPhenoWeight[which(type_index == j),iii]/pracma::Norm(OmicsPhenoWeight[which(type_index == j),iii])
    }

    
  }  

  # set part of the between-omics weight to 0
  BetweenOmicsWeight[OmicsPhenoWeight == 0] <- 0
  # remove all zero columns
  if (SubsamplingNum > 1)
  {
    # find zero columns and NAN columns
    zero_cols <- which(apply(OmicsPhenoWeight, 2, function(x) all(x == 0) | any(is.nan(x))))
    # remove all these columns
    BetweenOmicsWeight <- BetweenOmicsWeight[,-zero_cols]
    OmicsPhenoWeight <- OmicsPhenoWeight[,-zero_cols]
  }  
  # aggregate canonical weight (trade-off: between-omics, omics-phenotype)
  CCWeight <- (Between_Discriminate_Ratio[1]/sum(Between_Discriminate_Ratio)) * BetweenOmicsWeight + 
    (Between_Discriminate_Ratio[2]/sum(Between_Discriminate_Ratio)) * OmicsPhenoWeight
  # set row names
  row.names(CCWeight) <- colnames(X_all)
  
  return(CCWeight)
  }else{
    if (SubsamplingNum != 1)
      stop('Supsampling number must be 1 when evaluating the classifier.')
    
    # subset selected features
    X_subset <- X_all[ ,which(BetweenOmicsWeight[,1] != 0)]
    # run omics-phenotype SmCCA based on selected molecular features
    has_error <- FALSE
    tryCatch({
      # run omics-phenotype SmCCA based on selected molecular features
      out <- spls::splsda(x = X_subset, y = as.matrix(Y), K = ncomp_pls, eta = LambdaPheno, kappa=0.5,
                          classifier = 'logistic', scale.x=TRUE)
      # define testing data
      X_all_test <- rlist::list.cbind(testData)
      X_subset_test <- X_all_test[ ,which(BetweenOmicsWeight[,1] != 0)]
      out.data <- matrix(0, nrow = ncol(X_subset), ncol = ncomp_pls)
      out.data[out[["spls.fit"]][["A"]], ] <- out[["W"]]
      out.test <- as.matrix(X_subset_test) %*% out.data
      colnames(out[['T']]) <- colnames(out.test) <- paste0('lat', 1:ncomp_pls)
      return(list(out.train = out[['T']], out.test = out.test))
    },
    error = function(e) {
      cat("Caught an error:", e$message, "on iteration", "\n")
      has_error <- TRUE
    })
    # skip current iteration if error occurs
    if (has_error) {
      out.train <- matrix(0, nrow = nrow(X_all), ncol = ncomp_pls)
      out.test <- matrix(0, nrow = nrow(X_all), ncol = ncomp_pls)
      return(list(out.train, out.test))
    }
    
  }
  
  
}  

