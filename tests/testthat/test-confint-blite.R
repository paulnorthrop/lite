#context("confint.blite")

# Check that confint.blite() behaves as expected

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Fit model
  cdata <- exdex::cheeseboro
  cpost <- blite(cdata, u = 45, k = 3)

  # Perform the tests
  # Confidence intervals for all parameters
  res <- confint(cpost)
  test_that("confint: correct dim, all parameters", {
    testthat::expect_equal(dim(res), c(4, 2), ignore_attr = TRUE)
  })
  # Confidence intervals for GP parameters only
  res <- confint(cpost, parm = "gp")
  test_that("confint: correct dim, GP parameters", {
    testthat::expect_equal(dim(res), c(2, 2), ignore_attr = TRUE)
  })
  # Confidence interval for pu only
  res <- confint(cpost, parm = "pu")
  test_that("confint: correct dim, pu parameter", {
    testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
  })
  # Confidence interval for sigmau only
  res <- confint(cpost, parm = "sigmau")
  test_that("confint: correct dim, sigmau parameter", {
    testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
  })
  # Confidence interval for xi only
  res <- confint(cpost, parm = "xi")
  test_that("confint: correct dim, xi parameter", {
    testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
  })
  # Confidence interval for theta only
  res <- confint(cpost, parm = "theta")
  test_that("confint: correct dim, theta parameter", {
    testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
  })
  # Confidence interval for pu and sigmau
  res <- confint(cpost, parm = c("pu", "sigmau"))
  test_that("confint: correct dim, pu and sigma parameters", {
    testthat::expect_equal(dim(res), c(2, 2), ignore_attr = TRUE)
  })
  # Confidence interval for pu and theta
  res <- confint(cpost, parm = c("pu", "theta"))
  test_that("confint: correct dim, pu and theta parameters", {
    testthat::expect_equal(dim(res), c(2, 2), ignore_attr = TRUE)
  })
  # Confidence interval for pu, theta and sigmau
  res <- confint(cpost, parm = c("pu", "theta", "sigmau"))
  test_that("confint: correct dim, pu, theta and sigmau parameters", {
    testthat::expect_equal(dim(res), c(3, 2), ignore_attr = TRUE)
  })

  # Supplying parm incorrectly
  test_that("confint: invalid parameter name", {
    testthat::expect_error(confint(cpost, parm = "invalid"))
  })
}
