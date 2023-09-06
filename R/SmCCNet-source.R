#' Compute the similarity matrix based on one or more canonical correlation
#' weight vectors.
#' 
#' Compute the similarity matrix based on the outer products of absolute
#' canonical correlation weights.
#' 
#' 
#' @param Ws A canonical correlation weight vector or matrix. If \code{Ws} is a
#' matrix, then each column corresponds to one weight vector.
#' @param P1 Total number of features for the first omics data type.
#' @param FeatureLabel If \code{FeatureLabel = NULL} (default), the feature
#' names will be \eqn{\{TypeI_1, \cdots, TypeI_{p_1}, TypeII_1, \cdots,
#' TypeII_{p-p_1}\}}, where \eqn{p_1 = }\code{P1}, and \eqn{p} is the total
#' number of omics features.
#' @return A \eqn{p\times p} symmetric non-negative matrix.
#' @examples
#' 
#' w <- matrix(rnorm(6), nrow = 3)
#' Ws <- apply(w, 2, function(x)return(x/sqrt(sum(x^2))))
#' abar <- getAbar(Ws, P1 = 2, FeatureLabel = NULL)
#' 
#' @export
getAbar <- function(Ws, P1 = NULL, FeatureLabel = NULL){
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
  
  diag(Abar) <- 0
  Abar <- Abar/max(Abar)
  
  if(is.null(colnames(Abar))){
    if(is.null(FeatureLabel)){
      if(is.null(P1)){
        stop("Need to provide FeatureLabel or the number of features 
                    for the first data type P1.")
      }else{
        p <- ncol(Abar)
        FeatureLabel <- c(paste0("TypeI_", seq_len(P1)), 
                          paste0("TypeII_", seq_len(p-P1)))
      }
    }
    colnames(Abar) <- rownames(Abar) <- FeatureLabel
  }
  
  return(Abar)
}






#' preprocess single omics dataset before running single omics SmCCNet
#' 
#' Data preprocess pipeline that: (1) cv filtering, (2) center or scale data
#' and (3) adjust for clinical covariates.
#' 
#' 
#' @param X dataframe with the size of nxp.
#' @param covariates dataframe with covariates to be adjusted for.
#' @param is_cv Whether to use coefficient of variation (CV) filter filter
#' data.
#' @param cv_quantile CV filtering quantile.
#' @param center Whether to center the dataset X1.
#' @param scale Whether to scale the dataset X1.
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







#' One omics SmCCNet when the outcome variable is categorical
#' 
#' One omics SmCCNet when there is only one omics, and the outcome variable is
#' categorical, implement SPLS instead of the regular canonical correlation
#' analysis.
#' 
#' 
#' @param X1 An \eqn{n\times p_1} data matrix (e.g. mRNA) with \eqn{p_1}
#' features and \eqn{n} subjects.
#' @param Trait An \eqn{n\times 1} trait data matrix for the same n subjects.
#' @param Lambda1 LASSO penalty parameter for \code{X1}. \code{Lambda1} needs
#' to be between 0 and 1.
#' @param s1 Proportion of mRNA features to be included, default at \code{s1 =
#' 0.7}. \code{s1} needs to be between 0 and 1.
#' @param SubsamplingNum The numebr of subsampling used.
#' @param K Number of hidden components for PLSDA
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




#' SmCCA Algorithm with Single Omics Data
#' 
#' Calculate the canonical correlation weights based on sparse canonical
#' correlation analysis (SCCA). Integrate one omics data type (and a
#' quantitative phenotype), and calculate the absolute canonical correlation
#' weights for the omics features using SCCA take into account a
#' phenotype/trait. SCCA maximizes the pairwise canonical correlation weights
#' between omics data and the trait. It requires the trait to be quantitative.
#' The algorithm is computed along with an omics feature subsampling scheme.
#' 
#' To choose SmCCA, set \code{NoTrait = FALSE, FilterByTrait = FALSE}.  To
#' choose SsCCA, set \code{NoTrait = FALSE, FilterByTrait = TRUE}. To choose
#' SCCA, set \code{Trait = NULL, NoTrait = TRUE}.
#' 
#' @param X1 An \eqn{n\times p_1} data matrix (e.g. mRNA) with \eqn{p_1}
#' features and \eqn{n} subjects.
#' @param Trait An \eqn{n\times 1} trait data matrix for the same n subjects.
#' @param Lambda1 LASSO penalty parameter for \code{X1}. \code{Lambda1} needs
#' to be between 0 and 1.
#' @param s1 Proportion of mRNA features to be included, default at \code{s1 =
#' 0.7}. \code{s1} needs to be between 0 and 1.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#' number leads to more accurate results, but at a higher cost.
#' @param trace Logical. Whether to display the CCA algorithm trace.
#' @return A canonical correlation weight matrix with \eqn{p_1+p_2} rows. Each
#' column is the canonical correlation weights based on subsampled \code{X1}
#' and \code{X2} features. The number of columns is \code{SubsamplingNum}.
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










#' Apply hierarchical tree cutting to the similarity matrix and extract multi/single-omics network modules.
#' 
#' Extract single-omics modules based on the similarity matrix.
#' 
#' 
#' @param Abar A similary matrix for all features (all omics data types).
#' @param CutHeight Height threshold for the hierarchical tree cutting. Default
#' is \eqn{1-0.1^{10}}.
#' @param PlotTree Logical. Whether to create a hierarchical tree plot.
#' @return A list of multi/single-omics modules.
#' @examples
#' 
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' feature_name <- paste0('feature_', 1:5)
#' abar <- getAbar(w, P1 = 5, FeatureLabel = feature_name)
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





#' Extract Network Summarization Result from Sub-network
#' 
#' Extract summarization scores (the first 3 prinicipal components) for
#' specified network module with given network size. The proteins will be
#' ranked based on PageRank algorithm, then the top k proteins (where k is the
#' specified subnetwork size) will be included into the final subnetwork to
#' generate the summarization score. For the PC score, the correlation with
#' respect to the phenotype of interest will be calculated and stored. In
#' addition, the correlation between individual proteins and phenotype of
#' interest will also be recorded. The final subnetwork adjacency matrix will
#' be stored into the user-specified working directory of interest.
#' 
#' 
#' @param Abar Adjacency matrix of size pxp extracted from the SmCCA step
#' @param CorrMatrix The correlation matrix calculaed based on X1, it should be
#' pxp as well
#' @param data the original protein data.
#' @param Pheno the original phenotype data
#' @param ModuleIdx the index of the network module that summarization score is
#' intended to be stored
#' @param min_mod_size the minimally possible subnetwork size for the given network module,
#' should be an integer from 1 to the largest possible size of the protein
#' network
#' @param max_mod_size the minimally possible subnetwork size for the given network module,
#' should be an integer from 1 to the largest possible size of the protein
#' network, and it needs to be greater than the specified minimally possible network size.
#' @param type A vector with length equal to total number of features in the adjacency matrix
#' indicating the type of data for each feature, for instance, it could be genes, or proteins.
#' @param damping damping parameter for the pagerank algorithm
#' @param method Either NetSHy'or 'PCA' indicating which summarization method to use
#' @param saving_dir Directory where user prefers to store the result
#' @return a file stored in the local designated directory, which contains the
#' following: (1) M: subnetwork adjacency matrix. (2) pca_score: a dataframe
#' contains regular PC scores for the first 3 PCs as well as the phenotype. (3)
#' rank_value: PageRank score for each individual feature in the subnetwork.\
#' (4) pca_hybrid: a list that contains: PC score, PC loadings and PC varianace
#' explained for hybrid PC method. (5) pca_hybrid: a list that contains: PC
#' score, PC loadings and PC varianace explained for hybrid PC method with
#' alpha = 0. (6) pc_correlation: Regular PC score's correlation with respect
#' to phenotype. (7) correlation_hybrid: Hybrid PC score's correlation with
#' respect to phenotype. (8) correlation_hybrid_zero: Hybrid PC score's
#' correlation with respect to phenotype with alpha = 0.  (9)
#' omics_correlation_data: Individual omics feature correlation with respect to
#' phenotype (10) pc_loading: Regular PC loadings.
#' @examples
#' library(SmCCNet)
#' set.seed(123)
#' w <- rnorm(20)
#' w <- w/sqrt(sum(w^2))
#' labels <- paste0('feature_', 1:20)
#' abar <- getAbar(w, P1 = 20, FeatureLabel = labels)
#' modules <- getOmicsModules(abar, CutHeight = 0.1)
#' x <- X1[ ,seq_len(20)]
#' corr <- stats::cor(x)
#' # display only example
#' # networkPruning(abar, corr, modules, data = x, Pheno = Y,folder = 'My_Example',
#' # ModuleIdx = 1, pheno_name = 'My_Example', mod_size = 3
#' # )
#' 
#' @export
#' 

networkPruning <- function(Abar, CorrMatrix, data,
                                             Pheno, type, ModuleIdx,
                                             min_mod_size, max_mod_size, damping = 0.9,
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
                                     npc = 3, is_alpha = FALSE)
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
                                   npc = 3, is_alpha = FALSE)
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








#' Saving cross-validation result as the cross-validation table into the
#' working directory and provide recommendation on the penalty term selection.
#' 
#' Saving cross-validation result as the cross-validation table into the
#' working directory
#' 
#' 
#' @param CVDir A directory where the result is stored.
#' @param SCCAmethod The canonical correlation analysis method that is used in
#' the model.
#' @param K number of folds for cross-validation.
#' @param NumSubsamp The numebr of subsampling used.
#' 
#' @export

aggregateCVSingle <- function(CVDir, SCCAmethod, K = 5, NumSubsamp = 500){
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
#'@description Implement NetSHy network summarization via a hybrid approach to 
#'summarize network by considering the network topology with laplacian matrix (and TOM matrix).
#'
#'@param X Data matrix of size N by P, where N is the number of subjects and P is the number 
#'of features.
#'@param A Corresponding adjacency matrix of size P by P.
#'@param is_alpha Whether TOM matrix is considered as part of the network topology, 
#'default is set to TRUE
#'@param npc Number of principal components used to summarize the network.
#'
#'
#'@return a list consists of [[1]] Subject-level network summarization score,
#'[[2]] Principal component importance information such as percent of variance explained,
#'(3) Principal component feature-level loadings.
#'
#'@examples 
#'# simulate omics data
#'OmicsData <- matrix(rnorm(200,0,1), nrow = 10, ncol = 20)
#'# simulate omics adjacency matrix
#'set.seed(123)
#'w <- rnorm(20)
#'w <- w/sqrt(sum(w^2))
#'abar <- getAbar(w, P1 = 2, FeatureLabel = NULL)
#'# extract NetSHy summarization score
#'netshy_score <- summarizeNetSHy(OmicsData, abar)
#'
#'@references 
#'Vu, T., Litkowski, E. M., Liu, W., Pratte, K. A., Lange, L., Bowler, R. P., ... & Kechris, K. J. (2023). NetSHy: network summarization via a hybrid approach leveraging topological properties. Bioinformatics, 39(1), btac818.
#'
#'@export

summarizeNetSHy = function(X, A, is_alpha = TRUE, npc=1){
  #X: data matrix (n, p)
  #A: corresponding adjacency matrix
  #pc_id: PC index
  g = igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  # laplacian
  L2 <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE,diag = FALSE) 
  L2 <- as.matrix(igraph::graph.laplacian(L2, normalized = F))
  # TOM
  TOM_A = WGCNA::TOMsimilarity(as.matrix(A), verbose = FALSE)
  
  alpha = igraph::graph.density(g)
  X= scale(X, center = TRUE, scale = TRUE)
  # weighted approach
  if (is_alpha == TRUE)
    temp  = (1-alpha)*(X %*% L2) + alpha*(X %*% TOM_A)
  else
    temp = (X %*% L2)
  temp = summary(stats::prcomp(temp))
  h_score = temp$x[,1:npc] 
  importance = temp$importance[,1:npc]
  loading = temp$rotation[,1:npc]
  return(list(h_score, importance, loading))
  
}


#######################################################
# Internal functions called by getRobustPseudoWeights #
#######################################################




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


#' @title SmCCA Algorithm with Binary Phenotype
#' @description SmCCNet algorithm with multi-omics data and binary phenotype. This is a stepwise approach 
#' that (1) Use SmCCA to identify relationship between omics (exlude phenotype), (2) within highly connected omics features
#' selected in step 1, identify relationship between these selected omics features and phenotype of interest with 
#' sparse PLS. Sparse PLS algorithm for binary outcome first compute PLS by assuming outcome is continuous,
#' and extract multiple latent factors, then use latent factors to fit logistic regression, and weight latent factor by
#' regression parameters. 
#' 
#' @param X A list of training omics data matrix with same set and order of subjects.
#' @param Y Outcome binary variable, it is recommended that users set the level of this 
#' variable to 0 and 1.
#' @param Between_Discriminate_Ratio A vector with length 2 indicating the relative importance
#' of between-omics relationship and omics-phenotype relationship. For instance a ratio of 1:1
#' means between-omics relationship and omics-phenotype relationship contribute equally to the canonical weight
#' construction process.
#' @param SubsamplingPercent  A vector corresponds to the number of omics data, specifying for each omics data, what is the 
#' percentage of omics feature being subsampled at each iteration.
#' @param CCcoef A vector of scaling factors only indicates the relationship between-omics (exclude omics-phenotype).
#' @param LambdaBetween A vector of sparsity penalty value for each omics data for running the between-omics SmCCA, each 
#' penalty term should be within the range of 0 and 1.
#' @param LambdaPheno A real number between 0 and 1, a penalty term when running the sparse PLS wit phenotype of interest, recommend to set to 0 or lower 
#' value such as 0.1.
#' @param SubsamplingNum Total number of subsamples, the more the better in terms of accuracy, but at a cost of 
#' computation time.
#' @param ncomp_pls Number of latent components for PLS, default set to 3.
#' @param EvalClassifier Whether algorithm is at the phase of evaluating classification performance or construct network, if
#' TRUE, latent factors from SPLSDA will be returned, if FALSE, canonical weight will be returned. Default is FALSE. 
#' @param testData A list of testing omics data matrix with same set and order of subjects, only used when EvalClassifier is set to TRUE.
#' @return A canonical correlation weight matrix with combining all omics data. Each
#' column is the canonical correlation weights based on subsampled X features. The number of columns is \code{SubsamplingNum}.
#' @examples
#' 
#' 
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#' X1 <- matrix(rnorm(600,0,1), nrow = 60)
#' X2 <- matrix(rnorm(600,0,1), nrow = 60)
#' Y_binary <- rbinom(60,1,0.5)
#' 
#'  Ws <- getRobustWeightsMultiBinary(list(X1,X2), Y_binary, 
#'  SubsamplingPercent = c(0.8,0.8), CCcoef = NULL,
#'  LambdaBetween = c(0.5,0.5), LambdaPheno = 0.1, SubsamplingNum = 10)
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

