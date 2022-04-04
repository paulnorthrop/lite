#context("logLik.lite")

# Check that calculating the value of the maximised log-likelihood by summing
# the values of max_loglik for the 3 constituent models gives the same value
# as the value of the (adjusted) log-likelihood evaluated at the MLE

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Make inferences
  cdata <- exdex::cheeseboro
  cfit <- flite(cdata, u = 45, k = 3)

  # 2 ways to find the maximised log-likelihood value
  res1 <- logLik(cfit)
  res2 <- check_logLik_flite(cfit)

  test_that("Fitted object and summary() agree", {
    testthat::expect_equal(res1, res2, ignore_attr = TRUE)
  })
}
