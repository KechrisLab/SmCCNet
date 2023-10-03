library(testthat)
library(SmCCNet)

test_that("data preprocessing is done correctly", {
  set.seed(123)
  X1 <- matrix(rnorm(60000,0,1), nrow = 200)
  processed_data <- dataPreprocess(X = as.data.frame(X1), covariates = NULL, is_cv = TRUE, 
                                    cv_quantile = 0.2, center = FALSE, scale = FALSE)
  processed_data <- as.matrix(processed_data)
  cv_value = apply(X1, 2, EnvStats::cv)
  cv_value = abs(cv_value)
  data <- X1[, which(cv_value > stats::quantile(cv_value, 0.2))]
  expect_equal(as.vector(processed_data), as.vector(data))
})
