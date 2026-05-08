# preprocess a omics dataset before running omics SmCCNet

Data preprocess pipeline to: (1) filter by coefficient of variation
(cv), (2) center or scale data and (3) adjust for clinical covariates.

## Usage

``` r
dataPreprocess(
  X,
  covariates = NULL,
  is_cv = FALSE,
  cv_quantile = 0,
  center = TRUE,
  scale = TRUE
)
```

## Arguments

- X:

  dataframe with the size of \\n\\ by \\p\\, where \\n\\ is the sample
  size and \\p\\ is the feature size.

- covariates:

  dataframe with covariates to be adjusted for.

- is_cv:

  Whether to use coefficient of variation filter (small cv filter out).

- cv_quantile:

  CV filtering quantile.

- center:

  Whether to center the dataset X.

- scale:

  Whether to scale the dataset X.

## Value

Processed omics data with the size of nxp.

## Examples

``` r

X1 <- as.data.frame(matrix(rnorm(600, 0, 1), nrow = 60))
covar <- as.data.frame(matrix(rnorm(120, 0, 1), nrow = 60))
processed_data <- dataPreprocess(X = X1, covariates = covar, is_cv = TRUE, 
cv_quantile = 0.5, center = TRUE, scale = TRUE)
```
