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
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered.
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \doi{10.1093/biomet/asm015}
#' @references Fawcett, L. and Walshaw, D. (2012), Estimating return levels
#'   from serially dependent extremes. \emph{Environmetrics}, \strong{23},
#'   272-283. \doi{10.1002/env.2133}
#' @references Northrop, P. J. and Chandler, R. E. (2021).
#'   chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment. R package
#'   version 1.1.5. \url{https://CRAN.R-project.org/package=chandwich}.
#' @references Northrop, P. J. and Christodoulides, C. (2022). exdex:
#' Estimation of the Extremal Index. R package version 1.1.1.
#' \url{https://CRAN.R-project.org/package=exdex/}.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @seealso \code{\link{flite}} for frequentist threshold-based inference for
#'   time series extremes.
#' @seealso \code{\link{returnLevel}} for frequentist threshold-based inference
#'   for return levels.
#' @docType package
#' @name lite
#' @import sandwich
#' @importFrom stats nobs vcov coef logLik confint
#' @importFrom graphics plot
NULL
