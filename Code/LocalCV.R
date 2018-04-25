##############################################################################################################################
# Author: W. Jenny Shi
#
# About: This script allows to run K-fold CV in parallel locally
#
# Variables needed:
#   K: number of folds in CV.
#   resultDir: result repository.
#   dataF: data file. Must contain feature datasets X1 (n by p1) and X2 (n by p2), and a quantitative phenotype Y (n by 1).
#   s1, s2: subsampling proprotions for X1 and X2, respectively.
#   subSample: number of subsamples.
#   P1P2: optional penalty pairs.
#   CCcoef: weights for the cannonical correlations CCor(X1, X2), CCor(X1, Y), and CCor(X2, Y) in SmCCA.
##############################################################################################################################

set.seed(12345)
library(parallel)


# Create a CV folder.
CVDir <- paste0(resultDir, K, "foldCV/" )
dir.create(CVDir)


# Split data into training and test sets.
load(dataF)
N <- nrow(X1); p1 <- ncol(X1); p2 <- ncol(X2)
foldIdx <- split(1:N, sample(1:N, K))
for(i in 1:K){
  iIdx <- foldIdx[[i]]
  x1.train <- scale(X1[-iIdx, ])
  x2.train <- scale(X2[-iIdx, ])
  yy.train <- scale(Y[-iIdx, ])
  x1.test <- scale(X1[iIdx, ])
  x2.test <- scale(X2[iIdx, ])
  yy.test <- scale(Y[iIdx, ])

  if(is.na(min(min(x1.train), min(x2.train), min(yy.train), min(x1.test), min(x2.test), min(yy.test)))){
    stop("Invalid scaled data.")
  }

  subD <- paste0(CVDir, "CV_", i, "/")
  dir.create(subD)
  save(x1.train, x2.train, yy.train, x1.test, x2.test, yy.test,
       s1, s2, P1P2, p1, p2, subSamp, CCcoef,
       file = paste0(subD, "Data.Rdata"))
}


# Run each fold in parallel. Need to ensure that at least K threads are available on the machine.
cl <- makeCluster(K, type = "FORK")
clusterExport(cl = cl, "CVDir")
parSapply(cl, 1:K, function(CVidx){
  source("../Code/ModifiedPMA.R")
  source("../Code/SmCCNetSource.R")
  
  subD <- paste0(CVDir, "CV_", CVidx, "/")
  load(paste0(subD, "Data.Rdata"))
  dir.create(paste0(subD, "SmCCA/"))
  
  RhoTrain <- RhoTest <- DeltaCor <- rep(0, nrow(P1P2))
  for(idx in 1:nrow(P1P2)){
    print(paste0("Running SmCCA on CV_", CVidx, " idx=", idx))
    
    l1 <- P1P2[idx, 1]; l2 <- P1P2[idx, 2]
    Ws <- getRobustPseudoWeights(x1.train, x2.train, yy.train, l1, l2, s1, s2, NoTrait = FALSE,
                  FilterByTrait = FALSE, Bipartite = TRUE, 
                  SubsamplingNum = subSamp, CCcoef = CCcoef)
    meanW <- rowMeans(Ws)
    v <- meanW[1:p1]; u <- meanW[p1 + 1:p2]

    rho.train <- cor(x1.train %*% v, x2.train %*% u) * CCcoef[1] + 
      cor(x1.train %*% v, yy.train) * CCcoef[2] + 
      cor(x2.train %*% u, yy.train) * CCcoef[3]
    rho.test <- cor(x1.test %*% v, x2.test %*% u) * CCcoef[1] + 
      cor(x1.test %*% v, yy.test) * CCcoef[2] + 
      cor(x2.test %*% u, yy.test) * CCcoef[3]
    
    RhoTrain[idx] <- round(rho.train, digits = 5)
    RhoTest[idx] <- round(rho.test, digits = 5)
    DeltaCor[idx] <- abs(rho.train - rho.test)
    
    if(idx %% 10 == 0){
      save(P1P2, RhoTrain, RhoTest, DeltaCor, idx, 
           file = paste0(subD, "temp.Rdata"))
    }
    
  }
  
  DeltaCor.all <- cbind(P1P2, RhoTrain, RhoTest, DeltaCor)
  colnames(DeltaCor.all) <- c("l1", "l2", "Training CC", "Test CC", "CC Pred. Error")
  write.csv(DeltaCor.all, 
            file = paste0(subD, "SmCCA/SCCA_", subSamp,"_allDeltaCor.csv"))
  
  system(paste0("rm ", subD, "temp.Rdata"))
  
  return(CVidx)

  }
)

# Close cluster
stopCluster(cl)

