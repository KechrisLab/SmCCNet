#' Canonical Correlation Value for SmCCA
#' 
#' Get canonical correlation value for SmCCA given canonical weight vectors and scaling factor
#'
#'@param X list. A list of data each with n subjects.
#'@param CCcoef vector. A vector indicating weights for each pairwise canonical correlation.
#'@param CCWeight list. A list consists of canonical weight vectors corresponds to 
#'each element of X.
#'@param Y matrix. A matrix consists of phenotype.
#'
#'@examples 
#'# getCanCorMulti()
#'@export

getCanCorMulti <- function(X, CCcoef, CCWeight, Y)
{
  # sort out all combinations
  paircomb <- utils::combn(length(X), 2)
  # projection
  omics_projection <- purrr::map(1:length(X), function(xx){
    X[[xx]] %*% CCWeight[[xx]]
  })
  omics_projection <- do.call(cbind, omics_projection)
  # calculate correlation between omics projection
  omics_cor <- stats::cor(omics_projection)
  # only extract the upper triangle of correlation matrix
  omics_cor <- omics_cor[upper.tri(omics_cor)]
  # calculate canonical correlation (pairwise between-omics)
  CCBetween <- sum(CCcoef[1:length(X)] * omics_cor)
  # calculate canonical correlation (pairwise omics-phenotype)
  pheno_cor <- sum(stats::cor(omics_projection, Y) * CCcoef[(length(X) + 1):length(CCcoef)])
  # adding two components together
  OverallCC <- CCBetween + pheno_cor
  
  return(OverallCC)
  
}


#' @title Scaling Factor Input
#' 
#' @description
#' Input the vector of the type of dataset, and return prompts that ask the user to supply the scaling 
#' factor intended for SmCCNet algorithm to prioritize the correlation structure of 
#' interest. All scaling factor values supplied should be numeric and nonnegative.
#' @param DataType A vector that contains the the data type of each omics data.
#' 
#' @examples
#' # not run
#' # scalingFactorInput(c('gene','mirna', 'phenotype'))
#' 
#' @export
#' 

scalingFactorInput <- function(DataType = NULL)
{
  # find all possible combinations
  combs <- utils::combn(length(DataType), 2)
  scalingfactors <- rep(0, ncol(combs))
  # ask user to supply the scaling factor for each combination:
  for (i in 1:ncol(combs))
  {
    # Initialize condition as not satisfied
    satisfied <- FALSE 
    
    suppressWarnings(while (!satisfied) {
      # Ask for user input
      input <- readline(prompt = paste0("Please enter scaling factor value for correlation between ",
                                        DataType[combs[1,i]], ' and ',DataType[combs[2,i]], ":"))
      
      # Check if the input is numeric
      if (is.numeric(as.numeric(input)) && !is.na(as.numeric(input))) {
        # Check if the input is nonnegative
        if (as.numeric(input) >= 0) {
          satisfied <- TRUE  # Set condition to satisfied
        } else {
          cat("Invalid input. Please enter a nonnegative number.\n")
        }
      } else {
        cat("Invalid input. Please enter a number.\n")
      }
    })
    scalingfactors[i] <- as.numeric(input)
    
  }  
  
  for (j in 1:ncol(combs))
  {
    cat( paste0("The scaling factor value for correlation between ",
                DataType[combs[1,j]], ' and ',DataType[combs[2,j]], " is: ", 
                scalingfactors[j], '\n'))
  }
  return(scalingfactors)
  
}

#'@title Evaluation of Classifier with Different Evaluation Metrics
#'@description  evaluate binary classifier's performance with respect to user-selected
#'metric (accuracy, auc score, precision, recall, f1).
#'@param obs Observed phenotype, vector consists of 0, 1.
#'@param pred Predicted probability of the phenotype, vector consists of any value between 0 and 1
#'@param EvalMethod Binary classifier evaluation method, should be one of the following:
#''accuracy', 'auc', 'precision', 'recall', and 'f1', default is set to accuracy.
#'@param BinarizeThreshold Cutoff thresold used to binarize the predicted probability, default is set 
#'to 0.5.
#'@param print_score Whether to print the evaluation score or not, default is set to 0.5.
#'
#'@examples 
#'# simulate observed binary phenotype
#'obs <- rbinom(100,1,0.5)
#'# simulate predicted probability
#'pred <- runif(100, 0,1)
#'# calculate the score
#'pred_score <- classifierEval(obs, pred, EvalMethod = 'f1', print_score = FALSE)
#'
#'@export
#'

classifierEval <- function(obs, pred, EvalMethod = 'accuracy', 
                           BinarizeThreshold = 0.5, print_score = TRUE)
{
  
  # specifying which method to use as evaluation metric
  if (EvalMethod == 'accuracy')
  {
    # binarize prediction
    pred_binary <- ifelse(pred > BinarizeThreshold, 1, 0)
    # calculate prediction accuracy:
    eval_score <- mean(as.factor(pred_binary) == as.factor(obs)) 
    
  }else if (EvalMethod == 'auc')
  {
    # calculate auc score
    suppressMessages(eval_score <- pROC::auc(as.factor(obs), pred))
    
  } else if (EvalMethod == 'precision')
  {
    # binarize prediction
    pred_binary <- ifelse(pred > BinarizeThreshold, 1, 0)
    # calculate precision
    TP <- sum(obs == 1 & pred_binary == 1)
    FP <- sum(obs == 0 & pred_binary == 1)
    eval_score <- TP / (TP + FP)
    
  } else if (EvalMethod == 'recall')
  {
    # binarize prediction
    pred_binary <- ifelse(pred > BinarizeThreshold, 1, 0)
    # calculate precision
    TP <- sum(obs == 1 & pred_binary == 1)
    FN <- sum(obs == 1 & pred_binary == 0)
    eval_score <- TP / (TP + FN)
    
  } else if (EvalMethod == 'f1')
  {
    # binarize prediction
    pred_binary <- ifelse(pred > BinarizeThreshold, 1, 0)
    # calculate f1 score
    TP <- sum(obs == 1 & pred_binary == 1)
    TN <- sum(obs == 0 & pred_binary == 0)
    FP <- sum(obs == 0 & pred_binary == 1)
    FN <- sum(obs == 1 & pred_binary == 0)
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    eval_score <- 2 * (precision * recall) / (precision + recall)
  } else{
    stop('correct evaluation metric should be suuplied, selections are among 
                 accuracy, auc, precision, recall, f1 score.')
  }
  # print out the evaluation score
  if (print_score == TRUE)
    cat(paste0('The ', EvalMethod, ' score is: ', eval_score, '\n'))
  # return the score
  return(eval_score)
}


#' @title Automated SmCCNet
#' 
#' @description Automated SmCCNet that automatically identify project problem (single-omics vs multi-omics),
#' and analysis method (CCA vs. PLS) based on the input data that is provided.
#' This method automatically preprocess data, choose scaling factors, subsampling percentage, and optimal penalty 
#' terms, then run through the complete SmCCNet pipeline without the requirement for
#' users to provide any information. This function will store all the subnetwork information to local directory 
#' user is providing, as well as return all the global network and evaluation information.  
#' 
#' 
#' @param X A list of matrices with same set and order of subjects 
#' @param Y Phenotype variable of either numeric or binary, for binary variable, for binary Y, it should be binarized to 0,1 before running this function.
#' @param AdjustedCovar A data frame of covariates of interest to be adjusted for through regressing-out approach
#' @param Kfold Number of folds for cross-validation
#' @param EvalMethod Selections among 'accuracy', 'auc', 'precision', 'recall', and 'f1', indicating for evaluating binary outcome, what's the metric to use
#' @param subSampNum Number of subsampling to run, the higher the better in terms of accuracy, but at a cost of computational time
#' @param BetweenShrinkage A real number > 0 that helps shrink the importance of omics-omics correlation component, the larger this number
#' is, the greater the shrinkage it is.
#' @param ScalingPen A numeric vector of length 2 used as the penalty terms for scaling factor determination method: default set to 0.1, and 
#' should be between 0 and 1.
#' @param DataType A vector indicating what type of data is each element of X, example would be c('gene', 'miRNA')
#' @param CutHeight A numeric value specifying the cut height for hierarchical clustering, should be between 0 and 1
#' @param cluster Determine if clustering algorithm should be applied, default is TRUE.
#' @param min_size Minimally possible subnetwork size after network pruning, default set to 10.
#' @param max_size Maximally possible subnetwork size after network pruning, default set to 100.
#' @param summarization Summarization method used for network pruning and summarization, should be either 'NetSHy' or 'PCA'.
#' @param saving_dir Directory where user would like to store the subnetwork results.
#' @param preprocess Whether the data preprocessing step should be conducted
#' @param ncomp_pls Number of components for PLS algorithm, only used when binary phenotype is given, default is set to 3.
#' @param tuneLength The total number of candidate penalty term values for each omics data, default is set to 5.
#' @param tuneRangeCCA A vector of length 2 that represents the range of candidate penalty term values for each omics data based on canonical correlation analysis, 
#' default is set to c(0.1,0.5).
#' @param tuneRangePLS A vector of length 2 that represents the range of candidate penalty term values for each omics data based on partial least squared discriminant analysis, 
#' default is set to c(0.5,0.9).
#' @param seed Random seed for result reproducibility, default is set to 123.
#' @return subnetwork modules are stored in local directory specified by user, this function will return the global network information, which 
#'  include global adjacency matrix, data correlation matrix, hierarchical clustering result, Omics Abundance data, and Cross-Validation Result.  
#' @examples
#' 
#' 
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#' X1 <- matrix(rnorm(60000,0,1), nrow = 200)
#' colnames(X1) <- paste('gene',1:300)
#' X2 <- matrix(rnorm(60000,0,1), nrow = 200)
#' colnames(X2) <- paste('protein',1:300)
#' Y <- matrix(rnorm(200,0,1), nrow = 200)
#' Y_binary <- rbinom(200,1,0.5)
#' ### single-omics PLS
#' # result <- fastAutoSmCCNet(X = list(X1), Y = as.factor(Y_binary), 
#' # Kfold = 3, subSampNum = 20, DataType = c('Gene'))
#' ### single-omics CCA
#' # result <- fastAutoSmCCNet(X = list(X1), Y = Y, Kfold = 3, subSampNum = 20, 
#' # DataType = c('Gene'))
#' ### multi-omics PLS
#' # result <- fastAutoSmCCNet(X = list(X1,X2), Y = as.factor(Y_binary), Kfold = 3, 
#' # subSampNum = 20, DataType = c('Gene', 'miRNA'), 
#' # saving_dir = getwd())
#' ### multi-omics CCA
#' # result <- fastAutoSmCCNet(X = list(X1,X2), Y = Y, Kfold = 3, subSampNum = 20, 
#' # DataType = c('Gene', 'miRNA'), saving_dir = getwd())
#' 
#' @export

fastAutoSmCCNet <- function(X, Y, AdjustedCovar = NULL, preprocess = FALSE, Kfold = 5, 
                            EvalMethod = 'accuracy', subSampNum, DataType, 
                            BetweenShrinkage = 2, ScalingPen = c(0.1,0.1),
                            CutHeight = 1 - 0.1^10,
                            cluster = TRUE, min_size = 10,
                            max_size = 100, summarization = 'NetSHy', saving_dir = getwd(), 
                            ncomp_pls = 3,
                            tuneLength = 5,
                            tuneRangeCCA = c(0.1,0.5),
                            tuneRangePLS = c(0.5,0.9),
                            seed = 123)
{
  
  # set random seed
  seed <- set.seed(seed)
  cat("\n")
  cat("**********************************\n")
  cat("* Welcome to Automated SmCCNet! *\n")
  cat("**********************************\n")
  cat("\n")
  
  # Preprocess?
  if (preprocess == TRUE)
  {
    cat("\n")
    cat("--------------------------------------------------\n")
    cat(">> Starting data preprocessing...\n")
    cat("--------------------------------------------------\n")
    cat("\n")
    if (!is.null(AdjustedCovar))
    {
      cat('Covariate(s) are provided, now regressing out covariate effects from omics data.', '\n')
      AdjustedCovar <- as.data.frame(AdjustedCovar)
    }
    # preprocess data
    X <- purrr::map(1:length(X), function(xx){
      dataPreprocess(X = as.data.frame(X[[xx]]), covariates = AdjustedCovar, is_cv = FALSE, 
      center = TRUE, scale = TRUE)
    })
    X <- lapply(X, as.matrix)
  }
  # Define problem: Single-Omics? Multi-Omics? Categorical (Binary? MultiLevel?)? Continuous?
  if (length(X) > 1)
  {
    AnalysisType <- 'multiomics'
  }else{
    if (length(X) == 1)
    {
      AnalysisType <- 'singleomics'
    }else{
      stop('User must provide a list of dataframe.')
    }
  }
  
  # check the outcome
  if (ncol(as.matrix(Y)) > 1)
  {
    method <- 'CCA'
  }else{
    if (is.factor(Y))
    {
      method <- 'PLS'
    }
    if (is.numeric(Y))
    {
      method <- 'CCA'
    }
  }
  cat(paste0('This project uses ', AnalysisType, ' ', method), '\n')
 
  
  # Automatically set CC coefficients: check pairwise canonical correlation
  if (AnalysisType == 'multiomics')
  {
    cat("\n")
    cat("--------------------------------------------------\n")
    cat(">> Now determining the scaling factor for multi-omics analysis...\n")
    cat("--------------------------------------------------\n")
    cat("\n")
    AllComb <- utils::combn(length(X), 2)
    ScalingFactor <- rep(0, ncol(AllComb))
    for (i in 1:ncol(AllComb))
    {
      # define the pair of matrices
      X_pair <- X[AllComb[,i]]
      # extract canonical weight
      CC_weight <- getCanWeightsMulti(X_pair, Trait = NULL, Lambda = ScalingPen, NoTrait = TRUE)
      # define scaling factor
      ScalingFactor[i] <- abs(stats::cor(X_pair[[1]]%*% CC_weight[[1]], X_pair[[2]]%*% CC_weight[[2]]))
      
    }
    # shrink between-omics scaling factor based on shrinkage parameter
    ScalingFactor <- ScalingFactor/BetweenShrinkage
    if (method == 'PLS')
    {
      # set between-discriminated ratio
      BDRatio <- c(mean(ScalingFactor), 1)
      cat('The between-omics and omics-phenotype importance ratio is: ', BDRatio[1], ':', BDRatio[2], '\n')
    }
    # include scaling factors for multi-omics-phenotype relationship
    ScalingFactorCC <- c(ScalingFactor, rep(1, length(X)))
    CombPheno <- utils::combn(length(X) + 1,2)
    NonPhenoIndex <- apply(CombPheno, 2, function(x){x[2] != (length(X) + 1)})
    ScalingFactorTemp <- rep(1, ncol(CombPheno))
    ScalingFactorTemp[NonPhenoIndex] <- ScalingFactor
    ScalingFactor <- ScalingFactorTemp
    # if multi-omics PLS is used
    if (method == 'PLS')
    {
      ScalingFactor <- ScalingFactor[NonPhenoIndex]
    }
    
    cat("\n")
    cat('The scaling factor selection is: ', ScalingFactor, '\n')
 
  } else{
    cat('single omics analysis, skip scaling factor.', '\n')
  }
  
  
  
  
  cat("\n")
  cat("--------------------------------------------------\n")
  cat(">> Determining the best penalty selection through cross-validation...\n")
  cat("--------------------------------------------------\n")
  cat("\n")
  # Automatically select a penalty candidates 
  if (method == 'CCA')
  {
    # define penalty columns
    pen <- matrix(0, nrow = tuneLength, ncol = length(X))
    # fill in penalty candidate
    for (i in 1:ncol(pen))
    {
      pen[,i] <- seq(from = tuneRangeCCA[1], 
                             to = tuneRangeCCA[2], 
                             length.out = tuneLength)
    }
    # convert matrix to a list of columns
    list_cols <- as.list(as.data.frame(pen))
    # generate all possible combinations
    PenComb <- do.call(expand.grid, list_cols)
  }
  if (method == 'PLS')
  {
    # check if single-omics or multiomics
    if (AnalysisType == 'singleomics')
    {
      PenComb <- seq(from = tuneRangePLS[1], 
                     to = tuneRangePLS[2], 
                     length.out = tuneLength)
    }
    if (AnalysisType == 'multiomics')
    {
      # define penalty columns
      pen <- matrix(0, nrow = tuneLength, ncol = length(X) + 1)
      # fill in penalty candidate
      for (i in 1:ncol(pen))
      {
        pen[,i] <- seq(from = tuneRangePLS[1], 
                       to = tuneRangePLS[2], 
                       length.out = tuneLength)
      }
      # convert matrix to a list of columns
      list_cols <- as.list(as.data.frame(pen))
      # generate all possible combinations
      PenComb <- do.call(expand.grid, list_cols)
    }
  }
  
  # Automatically select percent of subsampling
  SubsamplingPercent <- rep(0, length(X))
  for (i in 1:length(X))
  {
    if(ncol(X[[i]]) < 300)
    {
      SubsamplingPercent[i] <- 0.9
    }else{
      SubsamplingPercent[i] <- 0.7
    }
  }
  
  
  # split data into folds
  X <- lapply(X, scale)
  # check to see whether Y need to be scaled or not
  if (method == 'CCA')
  {
    Y <- scale(Y)
  }
  
  foldIdx <- split(1:nrow(X[[1]]), sample(1:nrow(X[[1]]), Kfold))
  folddata <- purrr::map(1:length(foldIdx), function(x){
    
    Y <- as.matrix(Y)
    X_train <- list()
    X_test <- list()
    Y_train <- list()
    Y_test <- list()
    for (i in 1:length(X))
    {
      X_train[[i]] <- X[[i]][-foldIdx[[x]],]
      X_test[[i]] <- X[[i]][foldIdx[[x]],]
    }
    Y_train <- Y[-foldIdx[[x]],]
    Y_test <- Y[foldIdx[[x]],]
    return(list(X_train = X_train, X_test = X_test,Y_train = Y_train, Y_test = Y_test))
  })
  names(folddata) <- paste0('fold_', 1:Kfold)
  
  
  
  # Run SmCCNet (for single-omics)
  if (AnalysisType == 'singleomics')
  {
    
    # case 1: continuous outcome
    if (method == 'CCA')
    {
      # define parallel computing multisession
      future::plan(future::multisession, workers = Kfold)
      CVResult <- furrr::future_map(1:Kfold, function(xx) {
        # select the kth fold
        omicsdata <- folddata[[xx]]
        PenComb <- PenComb[,1]
        RhoTrain <- rep(0, length(PenComb))
        RhoTest <- rep(0, length(PenComb))
        DeltaCor <- rep(0, length(PenComb))
        for(idx in 1:length(PenComb)){
          # choose the penalty term
          l1 <- PenComb[idx]
          # run single omics SmCCA with continuous outcome
          Ws <- getRobustWeightsSingle(omicsdata[[1]][[1]], as.matrix(omicsdata[[3]]), l1, 1,
                                             SubsamplingNum = 1)
          rho.train <-  stats::cor(omicsdata[[1]][[1]] %*% Ws, as.matrix(omicsdata[[3]]))
          
          rho.test <-  stats::cor(omicsdata[[2]][[1]] %*% Ws, as.matrix(omicsdata[[4]]))
          
          # store cross-validation result
          RhoTrain[idx] <- round(rho.train, digits = 5)
          RhoTest[idx] <- round(rho.test, digits = 5)
          DeltaCor[idx] <- abs(rho.train - rho.test)
          
          
          
        }
        # combine cross-validation result
        CVResult <- as.data.frame(cbind(RhoTrain, RhoTest, DeltaCor))
        return(CVResult)
      },.progress = TRUE,.options = furrr::furrr_options(seed = TRUE))
      
      # aggregate CV result and select the best penalty term
      AggregatedCVResult <- Reduce("+", CVResult) / length(CVResult)
      EvalMetric <- apply(AggregatedCVResult, 1, function(x) {x[3]/abs(x[2])})
      BestPen <- PenComb[which.min(EvalMetric),1]
      cat(paste0('\n','The best penalty term after ', Kfold,'-fold cross-validation is: ', BestPen), '\n')
      cat(paste0('with testing canonical correlation = ', round(AggregatedCVResult$RhoTest[which.min(EvalMetric)],3), 
                 ', and prediction error = ', round(AggregatedCVResult$DeltaCor[which.min(EvalMetric)],3)), '\n')
      
      cat('Now running single-omics CCA with best penalty term on the complete dataset.', '\n')
      # run single-omics CCA
      Ws <- getRobustWeightsSingle(X1 = X[[1]], Trait = as.matrix(as.numeric(Y)), Lambda1 = BestPen, 
                                         s1 = SubsamplingPercent, SubsamplingNum = subSampNum)
      # construct global network module
      Abar <- getAbar(Ws, P1 = ncol(X[[1]]), FeatureLabel = colnames(X[[1]]))
    }
    # case when there is categorical outcome
    if (method == 'PLS')
    {
      # check if it is non-binary
      if (length(summary(as.factor(Y))) > 2)
        stop('Currently not support non-binary phenotype.')
      
      
      # define parallel computing multisession
      future::plan(future::multisession, workers = Kfold)
      CVResult <- furrr::future_map(1:Kfold, function(xx) {
        # select the kth fold
        omicsdata <- folddata[[xx]]
        TrainMetric <- rep(0, length(PenComb))
        TestMetric <- rep(0, length(PenComb))
        for(idx in 1:length(PenComb)){
          # choose the penalty term
          l1 <- PenComb[idx]
          # run single omics SmCCA with continuous outcome
          Ws <- spls::splsda(omicsdata[[1]][[1]], omicsdata[[3]], K = ncomp_pls, eta = l1, kappa=0.5,
                             classifier= 'logistic', scale.x=FALSE)
          
          weight <- matrix(0,nrow = ncol(omicsdata[[1]][[1]]), ncol = ncomp_pls)
          weight[Ws[["A"]],] <- Ws[["W"]]
          
          # create training and testing data, and fit logistic regression model
          train_data <- data.frame(x = (omicsdata[[1]][[1]] %*% weight)[,1:ncomp_pls], y = as.factor(omicsdata[[3]]))
          test_data <- data.frame(x = (omicsdata[[2]][[1]] %*% weight)[,1:ncomp_pls])
          logisticFit <- stats::glm(y~., family = 'binomial',data = train_data)
          # make prediction for train/test set
          train_pred <- stats::predict(logisticFit, train_data, type = 'response')
          test_pred <- stats::predict(logisticFit, test_data, type = 'response')
          # specifying which method to use as evaluation metric
          if (EvalMethod == 'accuracy')
          {
            # binarize prediction
            train_pred <- ifelse(train_pred > 0.5, 1, 0)
            test_pred <- ifelse(test_pred > 0.5, 1, 0)
            # calculate prediction accuracy:
            acc.train <- mean(as.factor(train_pred) == as.factor(omicsdata[[3]])) 
            acc.test <- mean(as.factor(test_pred) == as.factor(omicsdata[[4]]))
            # store prediction metric result
            TrainMetric[idx] <- acc.train
            TestMetric[idx] <- acc.test
          }
          if (EvalMethod == 'auc')
          {
            # calculate auc score
            auc.train <- pROC::auc(as.factor(omicsdata[[3]]), train_pred)
            auc.test <- pROC::auc(as.factor(omicsdata[[4]]), test_pred)
            # store prediction metric result
            TrainMetric[idx] <- auc.train
            TestMetric[idx] <- auc.test
          }
          if (EvalMethod == 'precision')
          {
            # binarize prediction
            train_pred <- ifelse(train_pred > 0.5, 1, 0)
            test_pred <- ifelse(test_pred > 0.5, 1, 0)
            # calculate precision
            TP_train <- sum(omicsdata[[3]] == 1 & train_pred == 1)
            FP_train <- sum(omicsdata[[3]] == 0 & train_pred == 1)
            precision.train <- TP_train / (TP_train + FP_train)
            TP_test <- sum(omicsdata[[4]] == 1 & test_pred == 1)
            FP_test <- sum(omicsdata[[4]] == 0 & test_pred == 1)
            precision.test <- TP_test / (TP_test + FP_test)
            # store prediction metric result
            TrainMetric[idx] <- precision.train
            TestMetric[idx] <- precision.test
          }
          if (EvalMethod == 'recall')
          {
            # binarize prediction
            train_pred <- ifelse(train_pred > 0.5, 1, 0)
            test_pred <- ifelse(test_pred > 0.5, 1, 0)
            # calculate precision
            TP_train <- sum(omicsdata[[3]] == 1 & train_pred == 1)
            FN_train <- sum(omicsdata[[3]] == 1 & train_pred == 0)
            recall.train <- TP_train / (TP_train + FN_train)
            TP_test <- sum(omicsdata[[4]] == 1 & test_pred == 1)
            FN_test <- sum(omicsdata[[4]] == 1 & test_pred == 0)
            recall.test <- TP_test / (TP_test + FN_test)
            # store prediction metric result
            TrainMetric[idx] <- recall.train
            TestMetric[idx] <- recall.test
          }
          if (EvalMethod == 'f1')
          {
            # binarize prediction
            train_pred <- ifelse(train_pred > 0.5, 1, 0)
            test_pred <- ifelse(test_pred > 0.5, 1, 0)
            # calculate precision
            TP_train <- sum(omicsdata[[3]] == 1 & train_pred == 1)
            TN_train <- sum(omicsdata[[3]] == 0 & train_pred == 0)
            FP_train <- sum(omicsdata[[3]] == 0 & train_pred == 1)
            FN_train <- sum(omicsdata[[3]] == 1 & train_pred == 0)
            precision.train <- TP_train / (TP_train + FP_train)
            recall.train <- TP_train / (TP_train + FN_train)
            f1.train <- 2 * (precision.train * recall.train) / (precision.train + recall.train)
            TP_test <- sum(omicsdata[[4]] == 1 & test_pred == 1)
            TN_test <- sum(omicsdata[[4]] == 0 & test_pred == 0)
            FP_test <- sum(omicsdata[[4]] == 0 & test_pred == 1)
            FN_test <- sum(omicsdata[[4]] == 1 & test_pred == 0)
            precision.test <- TP_test / (TP_test + FP_test)
            recall.test <- TP_test / (TP_test + FN_test)
            f1.test <- 2 * (precision.test * recall.test) / (precision.test + recall.test)
            # store prediction metric result
            TrainMetric[idx] <- f1.train
            TestMetric[idx] <- f1.test
          }
         
          
          
          
        }
        # combine cross-validation result
        CVResult <- as.data.frame(cbind(TrainMetric, TestMetric))
        return(CVResult)
      },.progress = TRUE,.options = furrr::furrr_options(seed = TRUE))
      
      # aggregate CV result and select the best penalty term
      AggregatedCVResult <- Reduce("+", CVResult) / length(CVResult)
      BestPen <- PenComb[which.max(AggregatedCVResult$TestMetric)]
      cat(paste0('\n'))
      for (xx in 1:length(BestPen))
      {  
        cat(paste0('The best penalty term for omics ', xx, ' after ', Kfold,'-fold cross-validation is: ', BestPen[xx]), '\n')
      }
      cat(paste0('with testing ', EvalMethod,' score = ', round(AggregatedCVResult$TestMetric[which.max(AggregatedCVResult$TestMetric)],3)), '\n')
      cat('Now running single-omics PLS with best penalty term on the complete dataset.', '\n')
      
      # run single-omics PLS
      Ws <- getRobustWeightsSingleBinary(X1 = X[[1]], Trait = as.numeric(Y) - 1, Lambda1 = BestPen, K = ncomp_pls,
                                         s1 = SubsamplingPercent, SubsamplingNum = subSampNum)
      # construct global network module
      Abar <- getAbar(Ws, P1 = ncol(X[[1]]), FeatureLabel = colnames(X[[1]]))
      
      
    }
    
    
  }
  
  # Run SmCCNet (for multi-omics)
  if (AnalysisType == 'multiomics')
  {
    
    # case 1: continuous outcome
    if (method == 'CCA')
    {
      # define parallel computing multisession
      future::plan(future::multisession, workers = Kfold)
      CVResult <- furrr::future_map(1:Kfold, function(xx) {
        # select the kth fold
        omicsdata <- folddata[[xx]]
        RhoTrain <- rep(0, length(PenComb))
        RhoTest <- rep(0, length(PenComb))
        DeltaCor <- rep(0, length(PenComb))
        for(idx in 1:nrow(PenComb)){
          # choose the penalty term
          l1 <- PenComb[idx, ]
          # run single omics SmCCA with continuous outcome
          Ws <- getCanWeightsMulti(omicsdata[[1]], Trait = as.matrix(omicsdata[[3]]), Lambda = l1, NoTrait = FALSE, CCcoef = ScalingFactor)
          rho.train <-  getCanCorMulti(X = omicsdata[[1]], Y = omicsdata[[3]],CCWeight =  Ws, CCcoef = ScalingFactorCC)
          
          rho.test <-  getCanCorMulti(X = omicsdata[[2]], Y = omicsdata[[4]],CCWeight =  Ws, CCcoef = ScalingFactorCC)
          
          # store cross-validation result
          RhoTrain[idx] <- round(rho.train, digits = 5)
          RhoTest[idx] <- round(rho.test, digits = 5)
          DeltaCor[idx] <- abs(rho.train - rho.test)
          
          
          
        }
        # combine cross-validation result
        CVResult <- as.data.frame(cbind(RhoTrain, RhoTest, DeltaCor))
        return(CVResult)
      },.progress = TRUE,.options = furrr::furrr_options(seed = TRUE))
      # aggregate CV result and select the best penalty term
      AggregatedCVResult <- Reduce("+", CVResult) / length(CVResult)
      EvalMetric <- apply(AggregatedCVResult, 1, function(x) {x[3]/abs(x[2])})
      BestPen <- PenComb[which.min(EvalMetric),]
      cat(paste0('\n'))
      for (xx in 1:length(BestPen))
      {  
        cat(paste0('The best penalty term for omics ', xx, ' after ', Kfold,'-fold cross-validation is: ', BestPen[xx]), '\n')
      }
      cat(paste0('with testing canonical correlation = ', round(AggregatedCVResult$RhoTest[which.min(EvalMetric)],3), 
                 ', and prediction error = ', round(AggregatedCVResult$DeltaCor[which.min(EvalMetric)],3)), '\n')
      
      cat('Now running multi-omics CCA with best penalty term on the complete dataset.', '\n')
      # run multi-omics CCA
      Ws <- getRobustWeightsMulti(X, Trait = as.matrix(Y), NoTrait = FALSE,CCcoef = ScalingFactor,
                                   Lambda = BestPen,s = SubsamplingPercent, SubsamplingNum = subSampNum)
      # construct global network
      featurelabel <- unlist(lapply(X, colnames))
      Abar <- getAbar(Ws, P1 = ncol(X[[1]]), FeatureLabel = featurelabel)
    }
    # case when there is categorical outcome
    if (method == 'PLS')
    {
      # check if it is non-binary
      if (length(summary(as.factor(Y))) > 2)
        stop('Currently not support non-binary outcome.')
      
      
      # define parallel computing multisession
      future::plan(future::multisession, workers = Kfold)
      # running multisession parallel computing with multi-block PLS
      CVResult <- furrr::future_map(1:Kfold, function(xx) {
        # select the kth fold
        omicsdata <- folddata[[xx]]
        TrainMetric <- rep(0, nrow(PenComb))
        TestMetric <- rep(0, nrow(PenComb))
        for(idx in 1:nrow(PenComb)){
          # choose the penalty term
          l1 <- PenComb[idx,]
          # run multi-omics SmCCA with continuous outcome
          # Ws <- getRobustWeightsMultiBinary(omicsdata[[1]], as.numeric(omicsdata[[3]]), SubsamplingPercent = rep(1, length(omicsdata[[1]])), Between_Discriminate_Ratio = BDRatio,
          #                                        LambdaBetween = l1[1,1:(length(l1) - 1)], LambdaPheno = l1[1,length(l1)], SubsamplingNum = 1, CCcoef = ScalingFactor, ncomp_pls = ncomp_pls)
          
          suppressMessages(projection <- getRobustWeightsMultiBinary(omicsdata[[1]], 
                                                                   as.numeric(omicsdata[[3]]), 
                                                                   SubsamplingPercent=rep(1, length(omicsdata[[1]])),
                                                                   Between_Discriminate_Ratio = BDRatio,
                                                                   LambdaBetween = l1[1,1:(length(l1) - 1)], 
                                                                   LambdaPheno = l1[1,length(l1)], 
                                                                   SubsamplingNum = 1, 
                                                                   CCcoef = ScalingFactor,
                                                                   ncomp_pls = ncomp_pls, EvalClassifier = TRUE,
                                                                   testData = omicsdata[[2]]))
          # determine the type of dataset each feature belongs to
          # feature_size <- unlist(lapply(X, ncol))
          # types <- unlist(purrr::map((1:length(feature_size)), function(x){
          #  rep(DataType[x], feature_size[x])
          #}))
          # spliting the canonical weight vector by type
          #Ws_split <- split(Ws[,1], types)
          # project each omics dataset into the lower dimensional embedding (for both training and testing data)
          #omics_projection_train <- purrr::map(1:length(X), function(xx){
          #  omicsdata[[1]][[xx]] %*% Ws_split[[xx]]
          #})
          #omics_projection_test <- purrr::map(1:length(X), function(xx){
          #  omicsdata[[2]][[xx]] %*% Ws_split[[xx]]
          #})

          # create training and testing data, and fit logistic regression model
          #train_data <- data.frame(x = do.call(cbind, omics_projection_train), y = as.factor(omicsdata[[3]]))
          #test_data <- data.frame(x =  do.call(cbind, omics_projection_test))
          train_data <- data.frame(x = projection[[1]], y = as.factor(omicsdata[[3]]))
          test_data <- data.frame(x =  projection[[2]])
          # catching error when performing the logistic regression
          has_error <- FALSE
          tryCatch({
            # fit logistic regression model
            logisticFit <- stats::glm(y ~ ., family = 'binomial', data = train_data)
            # make prediction for train/test set
            train_pred <- stats::predict(logisticFit, train_data, type = 'response')
            test_pred <- stats::predict(logisticFit, test_data, type = 'response')
            # specifying which method to use as evaluation metric
            if (EvalMethod == 'accuracy')
            {
              # binarize prediction
              train_pred <- ifelse(train_pred > 0.5, 1, 0)
              test_pred <- ifelse(test_pred > 0.5, 1, 0)
              # calculate prediction accuracy:
              acc.train <- sum(train_pred == omicsdata[[3]])/length(train_pred)
              acc.test <- sum(test_pred == omicsdata[[4]])/length(test_pred)
              # store prediction metric result
              TrainMetric[idx] <- acc.train
              TestMetric[idx] <- acc.test
            }else if (EvalMethod == 'auc')
            {
              # calculate auc score
              auc.train <- pROC::auc(as.factor(omicsdata[[3]]), train_pred)
              auc.test <- pROC::auc(as.factor(omicsdata[[4]]), test_pred)
              # store prediction metric result
              TrainMetric[idx] <- auc.train
              TestMetric[idx] <- auc.test
            }else if (EvalMethod == 'precision')
            {
              # binarize prediction
              train_pred <- ifelse(train_pred > 0.5, 1, 0)
              test_pred <- ifelse(test_pred > 0.5, 1, 0)
              # calculate precision
              TP_train <- sum(omicsdata[[3]] == 1 & train_pred == 1)
              FP_train <- sum(omicsdata[[3]] == 0 & train_pred == 1)
              precision.train <- TP_train / (TP_train + FP_train)
              TP_test <- sum(omicsdata[[4]] == 1 & test_pred == 1)
              FP_test <- sum(omicsdata[[4]] == 0 & test_pred == 1)
              precision.test <- TP_test / (TP_test + FP_test)
              # store prediction metric result
              TrainMetric[idx] <- precision.train
              TestMetric[idx] <- precision.test
            }else if (EvalMethod == 'recall')
            {
              # binarize prediction
              train_pred <- ifelse(train_pred > 0.5, 1, 0)
              test_pred <- ifelse(test_pred > 0.5, 1, 0)
              # calculate precision
              TP_train <- sum(omicsdata[[3]] == 1 & train_pred == 1)
              FN_train <- sum(omicsdata[[3]] == 1 & train_pred == 0)
              recall.train <- TP_train / (TP_train + FN_train)
              TP_test <- sum(omicsdata[[4]] == 1 & test_pred == 1)
              FN_test <- sum(omicsdata[[4]] == 1 & test_pred == 0)
              recall.test <- TP_test / (TP_test + FN_test)
              # store prediction metric result
              TrainMetric[idx] <- recall.train
              TestMetric[idx] <- recall.test
            }else if (EvalMethod == 'f1')
            {
              # binarize prediction
              train_pred <- ifelse(train_pred > 0.5, 1, 0)
              test_pred <- ifelse(test_pred > 0.5, 1, 0)
              # calculate precision
              TP_train <- sum(omicsdata[[3]] == 1 & train_pred == 1)
              TN_train <- sum(omicsdata[[3]] == 0 & train_pred == 0)
              FP_train <- sum(omicsdata[[3]] == 0 & train_pred == 1)
              FN_train <- sum(omicsdata[[3]] == 1 & train_pred == 0)
              precision.train <- TP_train / (TP_train + FP_train)
              recall.train <- TP_train / (TP_train + FN_train)
              f1.train <- 2 * (precision.train * recall.train) / (precision.train + recall.train)
              TP_test <- sum(omicsdata[[4]] == 1 & test_pred == 1)
              TN_test <- sum(omicsdata[[4]] == 0 & test_pred == 0)
              FP_test <- sum(omicsdata[[4]] == 0 & test_pred == 1)
              FN_test <- sum(omicsdata[[4]] == 1 & test_pred == 0)
              precision.test <- TP_test / (TP_test + FP_test)
              recall.test <- TP_test / (TP_test + FN_test)
              f1.test <- 2 * (precision.test * recall.test) / (precision.test + recall.test)
              # store prediction metric result
              TrainMetric[idx] <- f1.train
              TestMetric[idx] <- f1.test
            }else{
              stop('correct evaluation metric should be suuplied, selections are among 
                 accuracy, auc, precision, recall, f1 score.')
            }
            
          },
          error = function(e) {
            cat("Caught an error:", e$message, "on iteration", i, "\n")
            has_error <- TRUE
          })
          
          if (has_error) {
            next  # Skip to the next iteration
          }
          
          

          
          
        }
        # combine cross-validation result
        CVResult <- as.data.frame(cbind(TrainMetric, TestMetric))
        return(CVResult)
      },.progress = TRUE,.options = furrr::furrr_options(seed = TRUE))
      
      # aggregate CV result and select the best penalty term
      AggregatedCVResult <- Reduce("+", CVResult) / length(CVResult)
      BestPen <- PenComb[which.max(AggregatedCVResult$TestMetric),]
      cat('\n')
      for (xx in 1:length(X))
      {  
        cat(paste0('The best penalty term for omics ', xx, ' after ', Kfold,'-fold cross-validation is: ', BestPen[xx]), '\n')
      }
      cat(paste0('and the best penalty term on classifier is: ', BestPen[length(X) + 1], '\n'))
      cat(paste0('with testing ', EvalMethod,' score = ', round(AggregatedCVResult$TestMetric[which.max(AggregatedCVResult$TestMetric)],3)), '\n')
      cat('Now running multi-omics PLS with best penalty term on the complete dataset.', '\n')
      
      # run single-omics PLS
      outcome <- as.matrix(as.numeric(Y) - 1)
      Ws <- getRobustWeightsMultiBinary(X, 
                                        outcome, 
                                        SubsamplingPercent= SubsamplingPercent,
                                        Between_Discriminate_Ratio = BDRatio,
                                        LambdaBetween = BestPen[1:length(X)], 
                                        LambdaPheno = BestPen[length(X)+1], 
                                        SubsamplingNum = subSampNum, 
                                        CCcoef = ScalingFactor,
                                        ncomp_pls = ncomp_pls, EvalClassifier = FALSE)
      
      # construct global network module
      featurelabel <- unlist(lapply(X, colnames))
      Abar <- getAbar(Ws, P1 = ncol(X[[1]]), FeatureLabel = featurelabel)
      
      
    }
    
    
  }

 
  # save the cross-validation result 
  PenComb <- as.data.frame(PenComb)
  colnames(PenComb) <- paste0('l', 1:ncol(PenComb))
  AggregatedCVResultPen <- as.data.frame(cbind(PenComb, AggregatedCVResult))
  utils::write.csv(AggregatedCVResultPen, paste0(saving_dir, '/', 'CVResult.csv'))
  # save cross-validation data
  save(folddata, file = paste0(saving_dir, '/', 'CVFold.Rdata'))
  # save global network result
  globalNetwork <- list(AdjacencyMatrix = Abar, Data = X, Phenotype = Y)
  save(globalNetwork,
       file = paste0(saving_dir, '/', 'globalNetwork.Rdata'))
  
  

  
  cat("\n")
  cat("--------------------------------------------------\n")
  cat(">> Now starting network clustering...\n")
  cat("--------------------------------------------------\n")
  cat("\n")
  if (cluster == TRUE)
  {
    # run hierarchical clustering algorithm
    OmicsModule <- getOmicsModules(Abar, CutHeight = CutHeight, PlotTree = FALSE)
  }else{
    cat(">> Cluster = FALSE, skip clustering...\n")
    OmicsModule <- list(1:nrow(Abar))
  }  
  # store feature label
  AbarLabel <- unlist(lapply(X, colnames))
  # calculate correlation matrix
  bigCor2 <- stats::cor(do.call(cbind, X))
  # Combined data
  combined_data <- do.call(cbind, X)
  
  cat('Clustering completed...')
  cat("\n")
  cat("\n")
  cat("--------------------------------------------------\n")
  cat(">> Now starting network pruning and summarization score extraction...\n")
  cat("--------------------------------------------------\n")
  cat("\n")
  # filter out network modules with insufficient number of nodes
  module_length <- unlist(lapply(OmicsModule, length))
  network_modules <- OmicsModule[module_length > min_size]
  cat(paste0('There are ', length(network_modules), ' network modules before pruning'))
  # define data type
  feature_size <- unlist(lapply(X, ncol))
  types <- unlist(purrr::map((1:length(feature_size)), function(x){
    rep(DataType[x], feature_size[x])
  }))
  
  if (length(network_modules) == 0)
  {
    stop('No network is identify after clustering, consider changing the tree-cut parameter or changing the clustering algorithm.')
  }
  
  # extract pruned network modules
  for(i in 1:length(network_modules))
  {
    cat('\n','Now extracting subnetwork for network module ', i,'\n')
    # define subnetwork
    abar_sub <- Abar[network_modules[[i]],network_modules[[i]]]
    cor_sub <- bigCor2[network_modules[[i]],network_modules[[i]]]
    # prune network module
    cat(paste0('Network module ', i, ' Result:'))
    networkPruning(Abar = abar_sub,CorrMatrix = cor_sub, type = types[network_modules[[i]]], data = combined_data[,network_modules[[i]]], 
                                     Pheno = as.numeric(Y), ModuleIdx = i, min_mod_size = min_size, 
                                     max_mod_size = min(nrow(abar_sub), max_size), method = summarization, 
                                     saving_dir = saving_dir)
  
  }
 
  cat("\n")
  cat("\n")
  cat("************************************\n")
  cat("*   Execution Finished!            *\n")
  cat("*   Automated SmCCNet Completed!   *\n")
  cat("************************************\n")
  cat("\n")
  
  cat("There are 4 files stored in the local directory: \n")
  cat("1. CVResult.csv: cross-validation result. \n")
  cat("2. CVFold.Rdata: omics and phenotype data splited into different folds. \n")
  cat("3. globalNetwork.Rdata: global network adjacency matrix, omics data, phenotype. \n")
  cat("4. size_a_net_b.Rdata: subnetworks result file after clustering and network pruning. \n")
  
  
  
  cat("SUBNETWORK RESULT FILE CONTAINS:\n")
  cat("1. correlation_sub: correlation matrix for the subnetwork.\n")
  cat("2. M: adjacency matrix for the subnetwork.\n")
  cat("3. omics_corelation_data: individual molecular feature correlation with phenotype.\n")
  cat("4. pc_correlation: first 3 PCs correlation with phenotype.\n")
  cat("5. pc_loading: principal component loadings.\n")
  cat("6. pca_x1_score: principal component score and phenotype data.\n")
  cat("7. mod_size: number of molecular features in the subnetwork.\n")
  cat("8. sub_type: type of feature for each molecular features.\n")
  cat("\n")
  
  cat("NOTE:\n")
  cat("1. The results are saved in the user-specified directory.\n")
  cat("   If not specified, use 'getwd()' to find the current working directory.\n")
  cat("\n")
  cat("2. The results are in 'size_a_net_b.Rdata', where:\n")
  cat("   'a' = network size, 'b' = module number.\n")
  cat("\n")
  cat("3. To prevent overwriting result files across projects, rename the result file after each run.\n")
  cat("\n")
  cat("************************************\n")
  return(list(AdjacencyMatrix = Abar, Data = X, ClusteringModules = OmicsModule, CVResult = AggregatedCVResult))
  
  
  
}

