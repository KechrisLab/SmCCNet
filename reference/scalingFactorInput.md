# Scaling Factor Input Prompt

Input the vector of the annotation of each type of dataset in the data
list X (e.g., `c('gene', 'protein')`), and return prompt ask the user to
supply the scaling factor for SmCCNet algorithm to prioritize the
correlation structures of interest. All scaling factor values supplied
should be numeric and nonnegative.

## Usage

``` r
scalingFactorInput(DataType = NULL)
```

## Arguments

- DataType:

  A character vector that contains the annotation of each type of omics
  dataset in X.

## Value

A numeric vector of scaling factors.

## Examples

``` r
if(interactive()){scalingFactorInput(c('gene','mirna', 'phenotype'))}
```
