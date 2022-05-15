#context("predict.blite")

# Check that
# (a) if we do not supply ny to blite nor predict.blite then we get an error
# (b) if we supply ny to both then the value supplied to predict.blite is used
# (c) hpd intervals are shorter than the corresponding equi-tailed intervals
# (d) if we replace all the simulated values of theta with 1 then the
#     estimated predictive distribution shifts towards larger values

got_exdex <- requireNamespace("exdex", quietly = TRUE)

if (got_exdex) {
  ### Cheeseboro wind gusts

  # Fit model
  cdata <- exdex::cheeseboro
  cpost <- blite(cdata, u = 45, k = 3)

  # (a) Not supplying ny
  test_that("predict.blite() no ny error case", {
    testthat::expect_error(predict(cpost))
  })

  # (b) Supplying ny twice
  cpost <- blite(cdata, u = 45, k = 3, ny = 31 * 24)
  # Note: 800 is incorrect for these data!
  pred2 <- predict(cpost, ny = 800)
  test_that("predict.blite() ny twice", {
    testthat::expect_equal(pred2$ny, 800)
  })

  # (c) HPD intervals are shorter than equi-tailed intervals
  equi_tailed <- predict(cpost)$long
  equi_tailed_length <- equi_tailed[1, 2] - equi_tailed[1, 1]
  hpd <- predict(cpost, hpd = TRUE)$short
  hpd_length <- hpd[1, 2] - hpd[1, 1]
  test_that("HPD shorter than equi-tailed", {
    testthat::expect_lt(hpd_length, equi_tailed_length)
  })

  # Reveal these tests once revdbayes 1.5.9 is on CRAN

#  # (d) Replace simulated values of theta with 1s
#  # (i) Predictive intervals
#  theta_eq_1 <- cpost
#  theta_eq_1[, "theta"] <- 1
#  hpd2 <- predict(theta_eq_1, hpd = TRUE)$short
#  test_that("If theta = 1 then lower predictive interval limit increases", {
#    testthat::expect_gt(hpd2[1, 1], hpd[1, 1])
#  })
#  test_that("If theta = 1 then upper predictive interval limit increases", {
#    testthat::expect_gt(hpd2[1, 2], hpd[1, 2])
#  })
#  # (ii) Predictive quantiles
#  q1 <- predict(cpost, type = "q", n_years = c(100, 1000))$y
#  q2 <- predict(temp, type = "q", n_years = c(100, 1000))$y
#  test_that("If theta = 1 then predictive quantiles increase", {
#    testthat::expect_true(all(q2 > q1))
#  })

}
