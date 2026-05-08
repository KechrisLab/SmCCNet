# Calculate similarity matrix based on canonical weights.

Compute the similarity matrix based on the outer products of absolute
canonical correlation weights, can be used for both single and
multi-omics setting.

## Usage

``` r
getAbar(Ws, FeatureLabel = NULL)
```

## Arguments

- Ws:

  A canonical correlation weight vector or matrix. If `Ws` is a matrix,
  then each column corresponds to one weight vector.

- FeatureLabel:

  A vector of feature labels for each feature in the adjacency matrix

## Value

A \\p\times p\\ symmetric non-negative matrix.

## Examples

``` r

w <- matrix(rnorm(6), nrow = 3)
Ws <- apply(w, 2, function(x)return(x/sqrt(sum(x^2))))
abar <- getAbar(Ws,  FeatureLabel = c('omics1', 'omics2', 'omics3'))
```
