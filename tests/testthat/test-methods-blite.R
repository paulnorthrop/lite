#context("blite methods")

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Fit model
  cdata <- exdex::cheeseboro
  u <- 45
  cpost <- blite(cdata, u = u, k = 3)

  # nobs.blite
  temp <- nobs(cpost)
  # Number of non-missing observations
  pnobs <- sum(!is.na(cdata), na.rm = TRUE)
  # Number of threshold exceedances
  gpnobs <- sum(cdata > u, na.rm = TRUE)
  test_that("nobs.blite", {
    testthat::expect_equal(nobs(cpost)[1:2], c(pnobs, gpnobs), ignore_attr = TRUE)
  })

  # coef.blite
  test_that("coef.blite", {
    testthat::expect_equal(dim(coef(cpost)), c(1, 4), ignore_attr = TRUE)
  })

  # vcov.blite
  test_that("vcov.blite", {
    testthat::expect_equal(dim(vcov(cpost)), c(4, 4), ignore_attr = TRUE)
  })

}
