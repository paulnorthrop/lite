#context("estfun")

# Check that estfun produces a score vector that sums to (very close to) 0

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Fit model
  cdata <- exdex::cheeseboro
  cfit <- flite(cdata, u = 45, k = 3)

  # GP inferences
  gp_object <- attr(cfit, "gp")
  gp_fit <- attr(gp_object, "original_fit")
  close_to_zero <- colSums(estfun(gp_fit))

  test_that("estfun.GP", {
    testthat::expect_equal(close_to_zero, c(0, 0), tolerance = 1e-5,
                           ignore_attr = TRUE)
  })

  # Bernoulli inferences
  bernoulli_object <- attr(cfit, "Bernoulli")
  bernoulli_fit <- attr(bernoulli_object, "original_fit")
  close_to_zero <- colSums(estfun(bernoulli_fit))

  test_that("estfun.Bernoulli", {
    testthat::expect_equal(close_to_zero, 0, tolerance = 1e-5,
                           ignore_attr = TRUE)
  })
}
