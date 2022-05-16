#' lite: Likelihood-Based Inference for Time Series Extremes
#'
#' Performs likelihood-Based inference for stationary time series extremes.
#' The general approach follows Fawcett and Walshaw (2012). Marginal extreme
#' value inferences are adjusted for cluster dependence in the data using the
#' methodology in Chandler and Bate (2007), producing an adjusted
#' log-likelihood for the model parameters.  A log-likelihood for the extremal
#' index is produced using the K-gaps model of Suveges and Davison (2010).
#' These log-likelihoods are combined to make inferences about return levels.
#'
#' @details The main functions are
#' \itemize{
#'   \item{\code{\link{flite}}: makes frequentist threshold-based inference for
#'   time series extremes to produce an adjusted log-likelihood for the model
#'   parameters.}
#'   \item{\code{\link{returnLevel}}: performs inference for return levels
#'   using the adjusted log-likelihood.}
#' }
#'
#' The main functions are \code{\link{flite}} and \code{\link{blite}}, which
#' perform frequentist and Bayesian inference for time series extremes,
#' respectively.
#'
#' See the vignettes
#' \code{vignette("lite-1-frequentist", package = "lite")} and
#' \code{vignette("lite-1-bayesian", package = "lite")}
#' for an overview of the package.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered.
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \doi{10.1093/biomet/asm015}
#' @references Fawcett, L. and Walshaw, D. (2012), Estimating return levels
#'   from serially dependent extremes. \emph{Environmetrics}, \strong{23},
#'   272-283. \doi{10.1002/env.2133}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @seealso \code{\link{flite}} for frequentist threshold-based inference for
#'   time series extremes.
#' @seealso \code{\link{returnLevel}} for frequentist threshold-based inference
#'   for return levels.
#' @seealso \code{\link{blite}} for Bayesian threshold-based inference for
#'   time series extremes.
#' @seealso \code{\link{predict.blite}} for predictive inference for the
#'   largest value observed in \eqn{N} years.
#' @docType package
#' @name lite
#' @import sandwich
#' @import revdbayes
#' @importFrom stats nobs vcov coef logLik confint predict
#' @importFrom graphics plot
NULL
