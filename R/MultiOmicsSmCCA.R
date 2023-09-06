#--------------------------------------------------------------------------
# Script Name: Multi-Omics SmCCA Source Code
# Version: 1.0
# Author: Weixuan Liu, Jenny Shi
# Date: 2023-09-02
# 
# Description:
# This R script is the source function of multi-omics Sparse multiple Canonical 
# Correlation Analysis (SmCCA) with/without subsampling.
#--------------------------------------------------------------------------


#' @title Calculate the canonical  weights for SmCCA.
#' @description Integrate multiple (can be more than 2) omics data type (and a quantitative phenotype), and calculate
#' the absolute canonical correlation weights for the omics features using
#' SmCCA. SmCCA takes into account a phenotype/trait.
#' SmCCA maximizes the total (weighted or unweighted) pairwise canonical
#' correlation weights between two omics data types and the trait. It requires
#' the trait to be quantitative. This function is the generalization of function
#' getRobustPseudoWeights().
#' 
#' @param X A list of omics data matrix with same set and order of subjects.
#' @param s A vector corresponds to the number of omics data, specifying for each omics data, what is the 
#' percentage of omics feature being subsampled at each iteration.
#' @param TraitWeight Whether to output the canonical weight for phenotype
#' @param Trait An \eqn{n\times 1} trait data matrix for the same n subjects.
#' @param Lambda A vector of LASSO penalty parameters for \code{X}. \code{Lambda} needs
#' to be between 0 and 1.
#' @param NoTrait Logical, default is \code{FALSE}. Whether trait information
#' is provided.
#' @param FilterByTrait Logical, default is \code{FALSE}. Whether only the top
#' (\eqn{80\%}) features with highest correlation to the trait will be assigned
#' nonzero weights. The choice of \eqn{80\%} is based on the PMA package.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#' number leads to more accurate results, but at a higher cost.
#' @param CCcoef Optional coefficients for the SmCCA pairwise canonical
#' correlations. If \code{CCcoef = NULL} (default), then the objective function
#' is the total sum of all pairwise canonical correlations. It can also be a
#' coefficient vector that follows the column order of \code{combn(K, 2)}.
#' @param trace Logical. Whether to display the CCA algorithm trace.
#' @return A canonical correlation weight matrix with \eqn{p_1+p_2} rows. Each
#' column is the canonical correlation weights based on subsampled \code{X1}
#' and \code{X2} features. The number of columns is \code{SubsamplingNum}.
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
                                   FilterByTrait = FALSE, SubsamplingNum = 1000,
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
#' @description Run Sparse Multiple Canonical Correlation Analysis (SmCCA) and return
#' canonical weight.
#' 
#' @param X A list of omics data each with n subjects.
#' @param Trait An n by k trait data for the same samples (k >=1).
#' @param Lambda A vector of LASSO pentalty parameters, need to be between 0 and 1,
#' should also have the same length as the length of X.
#' @param CCcoef A vector indicating weights for each pairwise canonical correlation.
#' @param NoTrait Logical. Whether trait information is provided to Trait will be assigned nonzero weights.
#' @param trace Logical. Whether to display CCA algorithm trace.
#' @param TraitWeight Logical. Whether to include canonical weight for trait, default is set to FALSE
#' 
#' @examples 
#' # not run
#' # X <- list(X1,X2)
#' # result <- getCanWeightsMulti(X, Trait = as.matrix(Y), Lambda = c(0.5,0.5), NoTrait = FALSE)
#' # result <- getCanWeightsMulti(X, Trait = NULL, Lambda = c(0.5,0.5), NoTrait = TRUE)
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
