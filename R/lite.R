#' lite: Likelihood-based Inference for Time Series Extremes
#'
#' Explain what the package does
#' @details
#' Add details
#'
#' The main function is \code{\link{flite}}, which performs frequentist
#' inference for time series extremes.
#'
#' See \code{vignette("lite-vignette", package = "lite")} for an overview of the
#' package.
#' @references Northrop, P. J. and Chandler, R. E. (2021).
#'   chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment. R package
#'   version 1.1.5. \url{https://CRAN.R-project.org/package=chandwich}.
#' @references Northrop, P. J. and Christodoulides, C. (2022). exdex:
#' Estimation of the Extremal Index. R package version 1.1.1.
#' \url{https://CRAN.R-project.org/package=exdex/}.
#' @docType package
#' @name lite
#' @import sandwich
#' @importFrom stats nobs vcov coef logLik
#' @importFrom graphics plot
NULL
