# NetSHy Summarization Score

Implement NetSHy network summarization via a hybrid approach (Vu et
al.,) to summarize network by considering the network topology with
Laplacian matrix.

## Usage

``` r
summarizeNetSHy(X, A, npc = 1)
```

## Arguments

- X:

  An \\n\times m\\ data matrix with \\m\\ features and \\n\\ subjects.

- A:

  Corresponding adjacency matrix of size \\p\\ by \\p\\.

- npc:

  Number of principal components used to summarize the network, default
  is set to 1.

## Value

A list consists of (1) subject-level network summarization score, (2)
principal component importance information: standard deviation, percent
of variance explained, and cumulative proportion of variance explained,
and (3) principal component feature-level loadings.

## References

Vu, Thao, Elizabeth M. Litkowski, Weixuan Liu, Katherine A. Pratte,
Leslie Lange, Russell P. Bowler, Farnoush Banaei-Kashani, and Katerina
J. Kechris. "NetSHy: network summarization via a hybrid approach
leveraging topological properties." Bioinformatics 39, no. 1 (2023):
btac818.

## Examples

``` r
# simulate omics data
OmicsData <- matrix(rnorm(200,0,1), nrow = 10, ncol = 20)
# simulate omics adjacency matrix
set.seed(123)
w <- rnorm(20)
w <- w/sqrt(sum(w^2))
featurelabel <- paste0('omics',1:20)
abar <- getAbar(w, FeatureLabel = featurelabel)
# extract NetSHy summarization score
netshy_score <- summarizeNetSHy(OmicsData, abar)
```
