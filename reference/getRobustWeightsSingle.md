# Single-omics SmCCA with Quantitative Phenotype

Compute aggregated (SmCCA) canonical weights for single omics data with
quantitative phenotype (subampling enabled).

## Usage

``` r
getRobustWeightsSingle(
  X1,
  Trait,
  Lambda1,
  s1 = 0.7,
  SubsamplingNum = 1000,
  trace = FALSE
)
```

## Arguments

- X1:

  An \\n\times p_1\\ data matrix (e.g. mRNA) with \\p_1\\ features and
  \\n\\ subjects.

- Trait:

  An \\n\times 1\\ trait (phenotype) data matrix for the same \\n\\
  subjects.

- Lambda1:

  LASSO penalty parameter for `X1`. `Lambda1` needs to be between 0 and
  1.

- s1:

  Proportion of features in `X1` to be included, default at `s1 = 0.7`.
  `s1` needs to be between 0 and 1, default is set to 0.7.

- SubsamplingNum:

  Number of feature subsamples. Default is 1000. Larger number leads to
  more accurate results, but at a higher computational cost.

- trace:

  Whether to display the CCA algorithm trace, default is set to FALSE.

## Value

A canonical correlation weight matrix with \\p_1\\ rows. Each column is
the canonical correlation weights based on subsampled `X1` features. The
number of columns is `SubsamplingNum`.

## Examples

``` r


## For illustration, we only subsample 5 times.
set.seed(123)

# Single Omics SmCCA
W1 <- getRobustWeightsSingle(X1, Trait = Y, Lambda1 = 0.05,
  s1 = 0.7, 
  SubsamplingNum = 5, trace = FALSE)
  
  
```
