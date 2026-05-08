# Canonical Correlation Value for SmCCA

Calculate canonical correlation value for SmCCA given canonical weight
vectors and scaling factor

## Usage

``` r
getCanCorMulti(X, CCcoef, CCWeight, Y)
```

## Arguments

- X:

  A list of data each with same number of subjects.

- CCcoef:

  A vector of scaling factors indicating weights for each pairwise
  canonical correlation.

- CCWeight:

  A list of canonical weight vectors corresponds to each data in \\X\\.

- Y:

  A phenotype matrix, should have only one column.

## Value

A numeric value of the total canonical correlation

## Examples

``` r
library(SmCCNet)
data("ExampleData")
getCanCorMulti(list(X1,X2), CCcoef = c(1,1,1), 
CCWeight = list(rnorm(500,0,1), rnorm(100,0,1)), Y = Y)
#> [1] 0.4186895
```
