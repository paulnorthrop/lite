#context("estfun")

# Check that


# Check that calculating the value of the maximised log-likelihood by summing
# the values of max_loglik for the 3 constituent models gives the same value
# as the value of the (adjusted) log-likelihood evaluated at the MLE

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Make inferences
  cdata <- exdex::cheeseboro
  cfit <- flite(cdata, u = 45, k = 3)

  gp_object <- attr(cfit, "gp")
  gp_fit <- attr(gp_object, "original_fit")
  close_to_zero <- colSums(estfun(gp_fit))

  test_that("Fitted object and summary() agree", {
    testthat::expect_equal(close_to_zero, c(0, 0), tolerance = 1e-5,
                           ignore_attr = TRUE)
  })
}
