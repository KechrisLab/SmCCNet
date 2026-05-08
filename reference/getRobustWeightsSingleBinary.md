# Single-omics SmCCA with Binary Phenotype

Compute aggregated (SmCCA) canonical weights for single omics data with
quantitative phenotype (subampling enabled).

## Usage

``` r
getRobustWeightsSingleBinary(
  X1,
  Trait,
  Lambda1,
  s1 = 0.7,
  SubsamplingNum = 1000,
  K = 3
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

  Proportion of mRNA features to be included, default at `s1 = 0.7`.
  `s1` needs to be between 0 and 1, default is set to 0.7.

- SubsamplingNum:

  Number of feature subsamples. Default is 1000. Larger number leads to
  more accurate results, but at a higher computational cost.

- K:

  Number of hidden components for PLSDA, default is set to 3.

## Value

A partial least squared weight matrix with \\p_1\\ rows. Each column is
the canonical correlation weights based on subsampled `X1` features. The
number of columns is `SubsamplingNum`.

## Examples

``` r


X <- matrix(rnorm(600,0,1), nrow = 60)
Y <- rbinom(60,1,0.5)
Ws <- getRobustWeightsSingleBinary(X1 = X, Trait = as.matrix(Y), Lambda1 = 0.8, 
0.7, SubsamplingNum = 10)
```
