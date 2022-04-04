#context("returnLevel")

# Check that an error is thrown if we ask plot.returnLevel to find a confidence
# interval that has a higher level of confidence than one already calculated

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Fit model
  cdata <- exdex::cheeseboro
  cfit <- flite(cdata, u = 45, k = 3)
  rl <- returnLevel(cfit, inc = 5, ny = 31 * 24)

  test_that("plot.returnLevel error case", {
    testthat::expect_error(plot(rl, level = 0.99))
  })
}
