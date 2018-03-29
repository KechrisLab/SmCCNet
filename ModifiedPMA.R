################################################################################
# Author: W. Jenny Shi 
# 
# About: This script builds upon the PMA package (by Daniela Witten et al) to 
#   allow coefficients for canonical correlation in the objective function and 
#   to include quantitative phenotype input in the sparse multiple canonical 
#   correlation analysis (SmCCA). 
#
# Note: The phenotype input is assumed to be continuous. 
#
################################################################################


myMultiCCA <- function(xlist, penalty=NULL, ws=NULL, niter=25, 
                       type="standard", ncomponents=1, trace=TRUE, 
                       standardize=TRUE, CCcoef = NULL){
  # Perform sparse multiple canonical correlation analysis (SmCCA) for a list 
  #   of data inputs.
  # 
  # xlist: A list of length K, where K is the number of data sets on which to 
  #   perform SmCCA. Dataset k should be a matrix of dimension $n \times p_k$,
  #   where $p_k$ is the number of features in data set k. If some quantitative 
  #   phenotype information is included in the list, it should be the last 
  #   element of the list and of dimension $n \times 1$.
  # penalty: The penalty terms to be used. Can be a single value (if the same 
  #   penalty term is to be applied to each data set) or a K-vector, indicating
  #   a different penalty term for each data set. There are 2 possible 
  #   interpretations for the penalty term: If type = "standard" then this is 
  #   an L1 bound on $w_k$, and it must be between 1 and $\sqrt(p_k)$ ($p_k$ is 
  #   the number of features in matrix $X_k$). If type = "ordered" then this is
  #   the parameter for the fussed lasso penalty on $w_k$.
  # type: Either "standard" or "ordered" to indicate if the columns are 
  #   unordered or ordered. If type = "standard", then a lasso penalty is aplied 
  #   to v, to enforce sparsity. If type = "ordered" (generally used for CGH
  #   data), then a fused lasso penalty is applied, to enforce both sparsity and
  #   soothness. This argument can be a vector of lenght K (if different data 
  #   sets are of different types) or it can be a single value "standard" or 
  #   "ordered" (if all data sets are of the same type).
  # ncomponents: Number of factors to output. Default is 1.
  # niter: Number of iterations to be perfromed. Default is 25.
  # ws: A list of length K. The kth element contains the first ncomponents 
  #   columns of the v matrix of the SVD of $X_k$. If NULL, then the SVD of 
  #   $X_1, ..., X_K$ will be computed inside the MultiCCA function. However,
  #   if you plan to run this function multiple times, then save a copy of the 
  #   argument so that it does not need to be re-computed.
  # trace: Logical. Whether to print out progress.
  # standardize: Whether to center and scale the columns of $X_1, ..., X_K$.
  #   Default is TRUE.
  # CCcoef: Optional coefficients for the pairwise canonical correlations (CC).  
  #   If CCcoef = NULL (default), then the objective function is the total sum  
  #   of all pairwise CC. It can also be a coefficient vector that follows the 
  #   column order of combn(K, 2).
  
  if(ncol(xlist[[length(xlist)]]) > 1){
    out <- MultiCCA(xlist, penalty=penalty, ws=ws, niter=niter, 
                    type=type, ncomponents=ncomponents, trace=ncomponents, 
                    standardize=standardize)
    
  }else{
    # The Kth data set in xlist is a one dimensional phenotype
    call <- match.call()
    K <- length(xlist)
    pair_CC <- combn(K, 2)
    num_CC <- ncol(combn(K, 2))
    
    # Check canonical correlation coefficients.
    if(is.null(CCcoef)){CCcoef <- rep(1, num_CC)}
    if(length(CCcoef) != num_CC){
      stop(paste0("Invalid coefficients for pairwise canonical correlations. 
                  Please provide ", num_CC, " numerical values for CCcoef. 
                  Be sure to match combn(K,2) column order."))
    }
    
    # Check data type.
    if(length(type)==1){# Expand a single type to a vector of length(xlist).
      if(type != "standard"){
        stop("The phenotype data must be continuous and follow the type 'standard'. ")
      }
      type <- rep(type, K)} 
    if(length(type)!=K){ 
      stop("Type must be a vector of length 1, or length(xlist)")}
    if(sum(type!="standard" & type!="ordered")>0){ 
      stop("Each element of type must be standard or ordered.")}
    
    # Standardize or not.
    if(standardize){xlist <- lapply(xlist, scale)}
    
    # Initialize weights.
    if(!is.null(ws)){
      makenull <- FALSE
      for(i in 1:(K-1)){
        if(ncol(ws[[i]])<ncomponents) makenull <- TRUE
      }
      if(makenull) ws <- NULL
    }
    if(is.null(ws)){
      ws <- list()
      for(i in 1:(K-1)){ 
        ws[[i]] <- matrix(svd(xlist[[i]])$v[,1:ncomponents], ncol=ncomponents)
      }
      ws[[K]] <- 1
    }
    ws.init <- ws
    
    # Check penalties.
    if(is.null(penalty)){
      penalty <- rep(NA, K)
      penalty[type=="standard"] <- 4 # this is the default value of sumabs
      for(k in 1:(K-1)){
        if(type[k]=="ordered"){
          v <- svd(xlist[[k]])$v[,1]
          penalty[k] <- ChooseLambda1Lambda2(v)
        }
      }
    }
    if(length(penalty)==1) penalty <- rep(penalty, K)
    if(sum(penalty<1 & type=="standard")){
      stop("Cannot constrain sum of absolute values of weights to be less than 
           1.")
    }
    for(i in 1:(K-1)){
      if(type[i]=="standard" && penalty[i]>sqrt(ncol(xlist[[i]]))){
        stop("L1 bound of weights should be no more than sqrt of the number of 
             columns of the corresponding data set.", fill=TRUE)
      }
    }
    
    ws.final <- ws.init
    for(i in 1:(K-1)){ 
      ws.final[[i]] <- matrix(0, nrow=ncol(xlist[[i]]), ncol=ncomponents)
    }
    cors <- NULL
    for(comp in 1:ncomponents){
      ws <- list()
      for(i in 1:(K-1)) ws[[i]] <- ws.init[[i]][,comp]
      ws[[K]] <- 1
      curiter <- 1
      crit.old <- -10
      crit <- -20
      storecrits <- NULL
      
      while(curiter<=niter && abs(crit.old-crit)/abs(crit.old)>.001 &&
            crit.old!=0){
        crit.old <- crit
        crit <- myGetCrit(xlist, ws, pair_CC, CCcoef)
        storecrits <- c(storecrits,crit)
        if(trace) cat(curiter, fill=FALSE)
        curiter <- curiter+1
        for(i in 1:(K-1)){
          ws[[i]] <- myUpdateW(xlist, i, K, penalty[i], ws, type[i], ws.final,
                             pair_CC, CCcoef)
        }
      }
      
      for(i in 1:(K-1)) ws.final[[i]][,comp] <- ws[[i]]
      cors <- c(cors, myGetCors(xlist, ws, pair_CC, CCcoef))
    }
    
    out <- list(ws=ws.final, ws.init=ws.init, K=K, call=call, type=type, 
                penalty=penalty, cors=cors)
    class(out) <- "MultiCCA"
  }
  return(out)
}





############################################
# Internal functions called by myMultiCCA. #
############################################

myUpdateW <- function(xlist, i, K, sumabsthis, ws, type="standard", ws.final,
                    pair_CC, CCcoef){
  # Update canonical weights for the ith data set.
  # 
  # xlist: Data list.
  # i: Data set index.
  # K: Total number of data sets.
  # sumabsthis: Penalty for the ith data set.
  # ws: First ncomponents columns of the v matrix of the SVD of $X_k$'s.
  # type: Type of data sets.
  # ws.final: Current weight matrix. 
  # pair_CC: The output of combn(K, 2). A two-row table that include indices for
  #   pairwise canonical correlaltions between members of xlist.
  # CCcoef: Optional coefficients for the pairwise canonical correlations (CC).  

  
  tots0 <- sapply(1:length(CCcoef), function(x){
    pairx <- pair_CC[ , x]
    
    if(pairx[1] != i & pairx[[2]] != i){
      y <- 0
    }else{
      if(pairx[1] == i){j <- pairx[2]}
      if(pairx[2] == i){j <- pairx[1]}
      Xi <- xlist[[i]]
      Xj <- xlist[[j]]
      
      # diagmat is the diagonal correlation matrix calculated using previous 
      #   canonical directions.
      # If phenotype is included, only the first canonical direction is used. 
      #   diagmat is therefore a zero matrix.
      diagmat <- (t(ws.final[[i]])%*%t(Xi))%*%(Xj%*%ws.final[[j]])
      diagmat[row(diagmat)!=col(diagmat)] <- 0
      y <- t(Xi)%*%(Xj%*%ws[[j]]) - 
        ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
      y <- y * CCcoef[x]
    }
    
    return(y)
  })
  tots <- rowSums(tots0)
  
  if(type=="standard"){
    sumabsthis <- BinarySearch(tots, sumabsthis)
    w <- soft(tots, sumabsthis)/l2n(soft(tots, sumabsthis))
  } else {
    tots <- as.numeric(tots)
    tots <- tots/mean(abs(tots)) 
    w <- FLSA(tots,lambda1=sumabsthis,lambda2=sumabsthis)[1,1,]
    #    flsa.out <- diag.fused.lasso.new(tots,lam1=sumabsthis)
    #    lam2ind <- which.min(abs(flsa.out$lam2-sumabsthis))
    #    w <- flsa.out$coef[,lam2ind]
    w <- w/l2n(w)
    w[is.na(w)] <- 0
  }
  return(w)
}


myGetCrit <- function(xlist, ws, pair_CC, CCcoef){
  # Compute the matrix form SmCCA objective function value for given weights.
  
  crits <- apply(pair_CC, 2, function(x){
    i <- x[1]
    j <- x[2]
    y <- t(ws[[i]])%*%t(xlist[[i]])%*%xlist[[j]]%*%ws[[j]]
    return(y)
  })
  crit <- sum(crits * CCcoef)
  return(crit)
}


myGetCors <- function(xlist, ws, pair_CC, CCcoef){
  # Compute total weighted canonical correlations for given weights.
  
  CCs <- apply(pair_CC, 2, function(x){
    i <- x[1]
    j <- x[2]
    y <- cor(xlist[[i]]%*%ws[[i]], xlist[[j]]%*%ws[[j]]) 
    if(is.na(y)){y <- 0}
    return(y)
  })
  
  Cors <- sum(CCs * CCcoef)
  return(Cors)
}


BinarySearch <- function(argu,sumabs){
  #  Update sumabs so that the L1 norm of argu equals given penalty.
  
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 150){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}


l2n <- function(vec){
  # Computes the L2 norm. If the norm is zero, set it to 0.05.
  
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}


soft <- function(x,d){
  # Soft thresholding. 
  
  return(sign(x)*pmax(0, abs(x)-d))
}


##################################################
# Additional internal functions for ordered data #
##################################################
# ChooseLambda1Lambda2( )
# FLSA()




