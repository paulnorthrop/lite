#context("thetahat = 1 when k = 0")

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Make inferences
  cdata <- exdex::cheeseboro
  cfit <- flite(cdata, u = 45, k = 0)

  test_that("thetahat = 1 when k = 0", {
    testthat::expect_equal(coef(cfit)[4], 1, ignore_attr = TRUE)
  })
}
