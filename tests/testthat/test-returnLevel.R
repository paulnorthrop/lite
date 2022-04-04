#context("returnLevel")

# Check that an error is thrown if we ask plot.returnLevel to find a confidence
# interval that has a higher level of confidence than one already calculated.
# Also check that if we ask for a lower level then the confidence interval
# is shorter.

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Fit model
  cdata <- exdex::cheeseboro
  cfit <- flite(cdata, u = 45, k = 3)
  rl <- returnLevel(cfit, inc = 5, ny = 31 * 24)

  # Throw error
  test_that("plot.returnLevel error case", {
    testthat::expect_error(plot(rl, level = 0.99))
  })

  # Shorter interval
  oldrl <- rl$rl_prof
  newrl <- plot(rl, level = 0.9)
  test_that("old vs new RL: MLEs agree", {
    testthat::expect_equal(oldrl["mle"], newrl["mle"])
  })
  test_that("old vs new RL: new lower is higher", {
    testthat::expect_gt(newrl["lower"], oldrl["lower"])
  })
  test_that("old vs new RL: new upper is lower", {
    testthat::expect_lt(newrl["upper"], oldrl["upper"])
  })

  # Not supplying ny
  test_that("returnLevel no ny error case", {
    testthat::expect_error(returnLevel(cfit))
  })

  # Suppling ny twice
  cfit2 <- flite(cdata, u = 45, k = 3, ny = 1)
  rl2 <- returnLevel(cfit2, prof = FALSE, ny = 31 * 24)
  test_that("returnLevel ny twice", {
    testthat::expect_equal(rl2$ny, 31 * 24)
  })
}
