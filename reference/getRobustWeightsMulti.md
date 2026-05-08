# Run Sparse multiple Canonical Correlation Analysis and Obtain Canonical Weights (with Subsampling)

SmCCNet algorithm with multi-omics data and quantitative phenotype.
Calculate the canonical weights for SmCCA.

## Usage

``` r
getRobustWeightsMulti(
  X,
  Trait,
  Lambda,
  s = NULL,
  NoTrait = FALSE,
  SubsamplingNum = 1000,
  CCcoef = NULL,
  trace = FALSE,
  TraitWeight = FALSE
)
```

## Arguments

- X:

  A list of omics data each with n subjects.

- Trait:

  An \\n\times 1\\ trait (phenotype) data matrix for the same n
  subjects.

- Lambda:

  Lasso penalty vector with length equals to the number of omics data
  (\\X\\). `Lambda` needs to be between 0 and 1.

- s:

  A vector with length equals to the number of omics data (\\X\\),
  specifying the percentage of omics feature being subsampled at each
  subsampling iteration.

- NoTrait:

  Logical, default is `FALSE`. Whether trait information is provided.

- SubsamplingNum:

  Number of feature subsamples. Default is 1000. Larger number leads to
  more accurate results, but at a higher computational cost.

- CCcoef:

  Optional scaling factors for the SmCCA pairwise canonical
  correlations. If `CCcoef = NULL` (default), then the objective
  function is the total sum of all pairwise canonical correlations. This
  coefficient vector follows the column order of `combn(T+1, 2)`
  assuming there are T omics data and a phenotype data.

- trace:

  Whether to display the CCA algorithm trace, default is set to `FALSE`.

- TraitWeight:

  Whether to return canonical weight for trait (phenotype), default is
  set to `FALSE`.

## Value

A canonical correlation weight matrix with \\p = \sum\_{i} p_i\\ rows,
where \\p_i\\ is the number of features for the \\i\\th omics. Each
column is the canonical correlation weights based on subsampled
features. The number of columns is `SubsamplingNum`.

## Examples

``` r


## For illustration, we only subsample 5 times.
set.seed(123)
X1 <- matrix(rnorm(600,0,1), nrow = 60)
X2 <- matrix(rnorm(600,0,1), nrow = 60)
Y <- matrix(rnorm(60,0,1), nrow = 60)
# Unweighted SmCCA
result <- getRobustWeightsMulti(X = list(X1, X2), Trait = Y, NoTrait = FALSE,
Lambda = c(0.5, 0.5),s = c(0.7, 0.7), SubsamplingNum = 20)
  
```
