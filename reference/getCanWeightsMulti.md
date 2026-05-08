# Get Canonical Weight SmCCA Algorithm (No Subsampling)

Run Sparse multiple Canonical Correlation Analysis (SmCCA) and return
canonical weight vectors.

## Usage

``` r
getCanWeightsMulti(
  X,
  Trait = NULL,
  Lambda,
  CCcoef = NULL,
  NoTrait = TRUE,
  trace = FALSE,
  TraitWeight = FALSE
)
```

## Arguments

- X:

  A list of omics data each with n subjects.

- Trait:

  An \\n\\ by 1 trait (phenotype) data for the same samples.

- Lambda:

  Lasso penalty vector with length equals to the number of omics data
  (\\X\\). `Lambda` needs to be between 0 and 1.

- CCcoef:

  Optional scaling factors for the SmCCA pairwise canonical
  correlations. If `CCcoef = NULL` (default), then the objective
  function is the total sum of all pairwise canonical correlations. It
  follows the column order of `combn(T+1, 2)`, where `T` is the total
  number of omics data.

- NoTrait:

  Whether or not trait (phenotype) information is provided, default is
  set to `TRUE`.

- trace:

  Whether to display CCA algorithm trace, default is set to `FALSE`.

- TraitWeight:

  Whether to return canonical weight for trait (phenotype), default is
  set to `FALSE`.

## Value

A canonical weight vector with size of \\p\\ by 1.

## Examples

``` r
 # \donttest{
 x1 <- matrix(rnorm(1000), nrow = 50)
 x2 <- matrix(rnorm(1000), nrow = 50)
 y <- matrix(rnorm(50), nrow = 50)
 X <- list(x1,x2)
 result <- getCanWeightsMulti(X, Trait = y, Lambda = c(0.5,0.5), NoTrait = FALSE)
 result <- getCanWeightsMulti(X, Trait = NULL, Lambda = c(0.5,0.5), NoTrait = TRUE)
 cccoef <- c(1,10,10)
 result <- getCanWeightsMulti(X, Trait = y, CCcoef = cccoef, 
                              Lambda = c(0.5,0.5), NoTrait = FALSE)
                              # }
```
