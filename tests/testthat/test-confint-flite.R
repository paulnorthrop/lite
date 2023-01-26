#context("confint.flite")

# Check that confint.flite() behaves as expected

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Make inferences
  cdata <- exdex::cheeseboro
  cfit <- flite(cdata, u = 45, k = 3)

  # Perform the tests for profile = TRUE and profile = FALSE
  for (profile in c(TRUE, FALSE)) {
    # Confidence intervals for all parameters
    res <- confint(cfit, profile = profile)
    test_that("confint: correct dim, all parameters", {
      testthat::expect_equal(dim(res), c(4, 2), ignore_attr = TRUE)
    })
    # Confidence intervals for GP parameters only
    res <- confint(cfit, parm = "gp", profile = profile)
    test_that("confint: correct dim, GP parameters", {
      testthat::expect_equal(dim(res), c(2, 2), ignore_attr = TRUE)
    })
    # Confidence interval for pu only
    res <- confint(cfit, parm = "pu", profile = profile)
    test_that("confint: correct dim, pu parameter", {
      testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
    })
    # Confidence interval for sigmau only
    res <- confint(cfit, parm = "sigmau", profile = profile)
    test_that("confint: correct dim, sigmau parameter", {
      testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
    })
    # Confidence interval for xi only
    res <- confint(cfit, parm = "xi", profile = profile)
    test_that("confint: correct dim, xi parameter", {
      testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
    })
    # Confidence interval for theta only
    res <- confint(cfit, parm = "theta", profile = profile)
    test_that("confint: correct dim, theta parameter", {
      testthat::expect_equal(dim(res), c(1, 2), ignore_attr = TRUE)
    })
    # Confidence interval for pu and sigmau
    res <- confint(cfit, parm = c("pu", "sigmau"), profile = profile)
    test_that("confint: correct dim, pu and sigma parameters", {
      testthat::expect_equal(dim(res), c(2, 2), ignore_attr = TRUE)
    })
    # Confidence interval for pu and theta
    res <- confint(cfit, parm = c("pu", "theta"), profile = profile)
    test_that("confint: correct dim, pu and theta parameters", {
      testthat::expect_equal(dim(res), c(2, 2), ignore_attr = TRUE)
    })
    # Confidence interval for pu, theta and sigmau
    res <- confint(cfit, parm = c("pu", "theta", "sigmau"), profile = profile)
    test_that("confint: correct dim, pu, theta and sigmau parameters", {
      testthat::expect_equal(dim(res), c(3, 2), ignore_attr = TRUE)
    })
  }

  # Supplying parm incorrectly
  test_that("confint: invalid parameter name", {
    testthat::expect_error(confint(cfit, parm = "invalid"))
  })
}
