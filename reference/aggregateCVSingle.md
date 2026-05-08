# Aggregate and Save Cross-validation Result for Single-omics Analysis

Saves cross-validation results in a table with the user-defined
directory and outputs penalty term with the highest testing canonical
correlation, lowest prediction error, and lowest scaled prediction
error.

## Usage

``` r
aggregateCVSingle(CVDir, SCCAmethod = "SmCCA", K = 5, NumSubsamp = 500)
```

## Arguments

- CVDir:

  A directory where the result is stored.

- SCCAmethod:

  The canonical correlation analysis method that is used in the model,
  used to name cross-validation table file, default is set to 'SmCCA'.

- K:

  number of folds for cross-validation.

- NumSubsamp:

  Number of subsampling used.

## Value

A vector of length 3 with indices of the penalty term that (1) maximize
the testing canonical correlation, (2) minimize the prediction error and
(3) minimize the scaled prediction error.
