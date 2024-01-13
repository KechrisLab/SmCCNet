#' @title Run Sparse multiple Canonical Correlation Analysis and Obtain Canonical Weights (with Subsampling)
#' @description SmCCNet algorithm with multi-omics data and quantitative phenotype. Calculate the canonical weights for SmCCA.
#' 
#' @param X A list of omics data each with n subjects.
#' @param s A vector with length equals to the number of omics data (\eqn{X}), specifying the 
#' percentage of omics feature being subsampled at each subsampling iteration.
#' @param TraitWeight Whether to return canonical weight for trait (phenotype), default is set to \code{FALSE}.
#' @param Trait An \eqn{n\times 1} trait (phenotype) data matrix for the same n subjects.
#' @param Lambda Lasso penalty vector with length equals to the number of omics data (\eqn{X}). \code{Lambda} needs
#' to be between 0 and 1.
#' @param NoTrait Logical, default is \code{FALSE}. Whether trait information
#' is provided.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#' number leads to more accurate results, but at a higher computational cost.
#' @param CCcoef Optional scaling factors for the SmCCA pairwise canonical
#' correlations. If \code{CCcoef = NULL} (default), then the objective function
#' is the total sum of all pairwise canonical correlations. This 
#' coefficient vector follows the column order of \code{combn(T+1, 2)} assuming there are T omics data and a phenotype data.
#' @param trace Whether to display the CCA algorithm trace, default is set to \code{FALSE}.
#' @return A canonical correlation weight matrix with \eqn{p = \sum_{i} p_i} rows, where \eqn{p_i} is the number of features for the \eqn{i}th omics. Each
#' column is the canonical correlation weights based on subsampled features. The number of columns is \code{SubsamplingNum}.
#' @examples
#' 
#' 
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#' X1 <- matrix(rnorm(600,0,1), nrow = 60)
#' X2 <- matrix(rnorm(600,0,1), nrow = 60)
#' Y <- matrix(rnorm(60,0,1), nrow = 60)
#' # Unweighted SmCCA
#' result <- getRobustWeightsMulti(X = list(X1, X2), Trait = Y, NoTrait = FALSE,
#' Lambda = c(0.5, 0.5),s = c(0.7, 0.7), SubsamplingNum = 20)
#'   
#' @export
#' 

getRobustWeightsMulti <- function(X, Trait, Lambda,
                                   s = NULL, NoTrait = FALSE, 
                                   SubsamplingNum = 1000,
                                   CCcoef = NULL, trace = FALSE, TraitWeight = FALSE){
  
  if(sum(s == 0) > 1){
    stop("Subsampling proprotion needs to be greater than zero.")
  }else{
    if(sum(abs(s - 0.5) > 0.5) > 0){
      stop("Subsampling proportions can not exceed one.")}
  }
  
  if((sum(abs(Lambda - 0.5) > 0.5) > 0) | (sum(Lambda == 0) > 0)){
    stop("Invalid penalty parameter. Lambda1 needs to be between zero and one.")}
  
  p_data <- unlist(lapply(X, ncol))
  p <- sum(p_data)
  p_sub <- unlist(
    purrr::map(1:length(p_data), function(h){
      ceiling(p_data[h] * s[h])
    })
  )
  if (SubsamplingNum > 1)
  {  
  beta <- pbapply::pbsapply(seq_len(SubsamplingNum), function(x){
    
    # random subsampling
    samp <- purrr::map(1:length(p_data), function(h){
      sort(sample(seq_len(p_data[h]), p_sub[h], replace = FALSE))
    })
    # subset data
    x.par <- purrr::map(1:length(p_data), function(h){
      scale(X[[h]][ , samp[[h]]], center = TRUE, scale = TRUE)
    })
    # run SmCCA
    if (!is.null(Trait)){
      out <- getCanWeightsMulti(x.par, Trait, Lambda,
                             NoTrait = NoTrait,
                             trace = trace, CCcoef = CCcoef)
    }else{
      out <- getCanWeightsMulti(x.par,Trait = NULL, Lambda,
                             NoTrait = TRUE,
                             trace = trace, CCcoef = CCcoef)
    }
    
    
    w <- rep(0, p)
    # cumulative sum for features
    p_cum <- c(0, cumsum(p_data))
    for (h in 1:(length(p_cum) - 1))
    {
      w[samp[[h]] + p_cum[h]] <- out[[h]]
    }  

    coeff.avg <- w
    
    return(coeff.avg)
  })
  }else{  
    beta <- sapply(seq_len(SubsamplingNum), function(x){
      
      # random subsampling
      samp <- purrr::map(1:length(p_data), function(h){
        sort(sample(seq_len(p_data[h]), p_sub[h], replace = FALSE))
      })
      # subset data
      x.par <- purrr::map(1:length(p_data), function(h){
        scale(X[[h]][ , samp[[h]]], center = TRUE, scale = TRUE)
      })
      # run SmCCA
      if (!is.null(Trait)){
        out <- getCanWeightsMulti(x.par, Trait, Lambda,
                              NoTrait = NoTrait,
                              trace = trace, CCcoef = CCcoef)
      }else{
        out <- getCanWeightsMulti(x.par,Trait = NULL, Lambda,
                              NoTrait = TRUE,
                              trace = trace, CCcoef = CCcoef)
      }
      
      
      w <- rep(0, p)
      # cumulative sum for features
      p_cum <- c(0, cumsum(p_data))
      for (h in 1:(length(p_cum) - 1))
      {
        w[samp[[h]] + p_cum[h]] <- out[[h]]
      }  
      
      coeff.avg <- w
      
      return(coeff.avg)
    })
  }
  return(beta)
}



#' @title Get Canonical Weight SmCCA Algorithm (No Subsampling)
#' @description Run Sparse multiple Canonical Correlation Analysis (SmCCA) and return
#' canonical weight vectors.
#' 
#' @param X A list of omics data each with n subjects.
#' @param Trait An \eqn{n} by 1 trait (phenotype) data for the same samples.
#' @param Lambda Lasso penalty vector with length equals to the number of omics data (\eqn{X}). \code{Lambda} needs
#' to be between 0 and 1.
#' @param CCcoef Optional scaling factors for the SmCCA pairwise canonical
#' correlations. If \code{CCcoef = NULL} (default), then the objective function
#' is the total sum of all pairwise canonical correlations. It follows the column order of \code{combn(T+1, 2)}, where \code{T} is the total number of omics data.
#' @param NoTrait Whether or not trait (phenotype) information is provided, default is set to \code{TRUE}.
#' @param trace Whether to display CCA algorithm trace, default is set to \code{FALSE}.
#' @param TraitWeight Whether to return canonical weight for trait (phenotype), default is set to \code{FALSE}.
#' @return A canonical weight vector with size of \eqn{p} by 1.
#' @examples 
#' # This function is typically used as an internal function.
#' # It is also used when performing cross-validation,
#' # refer to multi-omics vignette for more detail.
#' # X <- list(X1,X2)
#' # result <- getCanWeightsMulti(X, Trait = as.matrix(Y), Lambda = c(0.5,0.5), NoTrait = FALSE)
#' # result <- getCanWeightsMulti(X, Trait = NULL, Lambda = c(0.5,0.5), NoTrait = TRUE)
#' # cccoef <- c(1,10,10)
#' # result <- getCanWeightsMulti(X, Trait = as.matrix(Y), CCcoef = cccoef, 
#' #                              Lambda = c(0.5,0.5), NoTrait = FALSE)
#' @export


# (INTERNAL)
getCanWeightsMulti <- function(X, Trait = NULL, Lambda, CCcoef = NULL,
                      NoTrait = TRUE, trace = FALSE, TraitWeight = FALSE){
  # Compute CCA weights.
  #
  # X: a list of omics data
  # Trait: An n by k trait data for the same samples (k >=1).
  # Lambda: LASSO pentalty parameters, need to be between 0 and 1.
  # CCcoef: A vector indicating weights for each pairwise canonical
  #   correlation.
  # NoTrait: Logical. Whether trait information is provided.
  #   to Trait will be assigned nonzero weights. The top 80% features are reserved.
  # trace: Logical. Whether to display CCA algorithm trace.
  
  for (j in 1:length(Lambda))
  {
    if(abs(Lambda[j] - 0.5) > 0.5){
      stop("Invalid penalty parameter. Lambda1 needs to be between zero and
           one.")}
  }  
  
  if(min(Lambda) == 0){
    stop("Invalid penalty parameter. Both Lambda1 and Lambda2 has to be
           greater than 0.")
  }
  if (NoTrait == TRUE)
  {
    # calculate penalty
    L <- unlist(purrr::map(1:length(X), function(i){
      max(1, sqrt(ncol(X[[i]])) * Lambda[[i]])
    }))
    # run SmCCA without phenotype
    out <- myMultiCCA(X, penalty = L,
                      CCcoef = CCcoef, trace = trace)
  }else{
    # calculate penalty
    L <- unlist(purrr::map(1:length(X), function(i){
      max(1, sqrt(ncol(X[[i]])) * Lambda[[i]])
    }))
    # add phenotype into data list
    X[[length(X)+1]] <- scale(Trait)
    # add penalty for phenotype
    L[length(X)] <- sqrt(ncol(Trait))
    out <- myMultiCCA(X, penalty = L,
                      CCcoef = CCcoef, trace = trace)
  }
  # store canonical weight
  if (TraitWeight == TRUE)
  {
    ws <- out$ws
  }else{
    if (NoTrait == TRUE)
    {
      ws <- out$ws
    }else{
      ws <- out$ws[-length(out$ws)]
    }
  } 
 
  
  return(ws)
}
