# Run Sparse multiple Canonical Correlation Analysis and Obtain Canonical Weights (with Subsampling)

SmCCNet algorithm with multi-omics data and binary phenotype. This is a
stepwise approach (1) use SmCCA to identify relationship between omics
(exlude phenotype), (2) within highly connected omics features selected
in step 1, identify relationship between these selected omics features
and phenotype of interest with sparse PLS. First, it computes PLSDA by
assuming outcome is continuous to extract multiple latent factors, then
uses latent factors to fit logistic regression, and weight latent factor
by regression parameters. Refer to multi-omics vignette for more detail.

## Usage

``` r
getRobustWeightsMultiBinary(
  X,
  Y,
  Between_Discriminate_Ratio = c(1, 1),
  SubsamplingPercent = NULL,
  CCcoef = NULL,
  LambdaBetween,
  LambdaPheno = NULL,
  SubsamplingNum = 1000,
  ncomp_pls = 3,
  EvalClassifier = FALSE,
  testData = NULL
)
```

## Arguments

- X:

  A list of omics data each with n subjects.

- Y:

  A vector of binary variable, user needs to set the level of this
  variable to 0 and 1.

- Between_Discriminate_Ratio:

  A vector with length 2 specifying the relative importance of
  between-omics relationship and omics-phenotype relationship. For
  instance a ratio of 1:1 (c(1,1) in the argument) means between-omics
  relationship and omics-phenotype relationship contribute equally to
  the canonical weights extraction.

- SubsamplingPercent:

  A vector with length equal to the number of omics data (`X`),
  specifying the percentage of omics feature being subsampled at each
  subsampling iteration.

- CCcoef:

  A vector of scaling factors only for between-omics relationship
  (exclude omics-phenotype). This coefficient vector follows the column
  order of `combn(T, 2)` when there are `T` omics data.

- LambdaBetween:

  A vector of sparsity penalty value for each omics data to run the
  between-omics SmCCA, each penalty term should be within the range of 0
  and 1.

- LambdaPheno:

  A penalty term when running the sparse PLS with phenotype, penalty
  term should be within the range of 0 and 1.

- SubsamplingNum:

  Number of feature subsamples. Default is 1000. Larger number leads to
  more accurate results, but at a higher computational cost, default is
  set to 1000.

- ncomp_pls:

  Number of latent components for PLS, default set to 3.

- EvalClassifier:

  If `TRUE`, the algorithm is at the phase of evaluating classification
  performance, and the latent factors from SPLSDA will be returned; if
  FALSE, the algorithm is at the phase of constructing multi-omics
  network, canonical weight will be returned. Default is set to `FALSE`.

- testData:

  A list of testing omics data matrix, should have the exact same order
  as data list X, only used when EvalClassifier is set to `TRUE` for
  performing cross-validation, refer to multi-omics vignette for detail.

## Value

If `EvalClassifier` is set to `FALSE`, a canonical correlation weight
matrix is returned with combined omics data. Each column is the
canonical correlation weights based on subsampled X features. The number
of columns is `SubsamplingNum`. If `EvalClassifier` is set to `TRUE`,
then latent factors from training and testing data will be returned for
classifier evaluation.

## Examples

``` r


## For illustration, we only subsample 5 times.
set.seed(123)
X1 <- matrix(rnorm(600,0,1), nrow = 60)
X2 <- matrix(rnorm(600,0,1), nrow = 60)
Y_binary <- rbinom(60,1,0.5)

Ws <- getRobustWeightsMultiBinary(list(X1,X2), Y_binary, 
      SubsamplingPercent = c(0.8,0.8), CCcoef = NULL,
      LambdaBetween = c(0.5,0.5), LambdaPheno = 0.1, SubsamplingNum = 10)
  
```
