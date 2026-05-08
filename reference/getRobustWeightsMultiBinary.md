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
  testData = NULL,
  verbose = FALSE
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
  between-omics relationship and omics-phenotype relationship.

- SubsamplingPercent:

  A vector with length equal to the number of omics data (`X`).

- CCcoef:

  A vector of scaling factors only for between-omics relationship.

- LambdaBetween:

  A vector of sparsity penalty value for each omics data.

- LambdaPheno:

  A penalty term when running the sparse PLS with phenotype.

- SubsamplingNum:

  Number of feature subsamples.

- ncomp_pls:

  Number of latent components for PLS.

- EvalClassifier:

  If TRUE, return latent factors for classification.

- testData:

  A list of testing omics data matrix.

- verbose:

  Logical; if TRUE, print progress/error messages, otherwise run
  silently.

## Value

Canonical weight matrix or latent projections depending on
EvalClassifier.
