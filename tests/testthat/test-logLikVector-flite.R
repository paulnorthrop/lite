#context("logLikVector.lite")

# Check that logLik(logLikVector) and logLik() agree

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Bernoulli
  # Set up data
  cdata <- c(exdex::cheeseboro)
  u <- 45
  exc <- cdata > u

  # Fit a Bernoulli distribution
  fit <- fitBernoulli(exc)

  # The logLik method sums the individual log-likelihood contributions.
  res1 <- logLik(logLikVector(fit))
  res2 <- logLik(fit)

  test_that("logLik vs logLik(logLikVector)", {
    testthat::expect_equal(res1, res2, ignore_attr = TRUE)
  })
}
