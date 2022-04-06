#' Functions for the \code{estfun} method
#'
#' Functions to calculate contributions to the score vector from individual
#' observations for a fitted model object.
#' @param x A fitted model object.
#' @param ... Further arguments.  None are currently used for
#'   \code{estfun.Bernoulli} or \code{estfun.GP}.
#' @details A \code{\link[sandwich]{estfun}} method is required by
#'   \code{\link[sandwich:vcovCL]{meatCL}} to calculate the
#'   \code{\link[sandwich]{meat}} in the sandwich covariance estimator on which
#'   the log-likelihood adjustments in \code{\link{flite}} are based.
#'   Specifically, \code{\link[sandwich:vcovCL]{meatCL}} is used to calculate
#'   the argument \code{V} passed to \code{\link[chandwich]{adjust_loglik}}.
#' @return An \eqn{n \times k}{n x k} numeric matrix containing contributions
#'   to the score function from \eqn{n} observations for each of the \eqn{k}
#'   parameters.
#'
#'   \strong{estfun.Bernoulli}: an \eqn{n \times 2}{n x 1} matrix, where
#'   \eqn{n} is the sample size, the length of the input \code{data} to
#'   \code{\link{fitBernoulli}}.  The column is named \code{prob}.
#'
#'   \strong{estfun.generalisedPareto}: an \eqn{n \times 2}{n x 2} matrix, where
#'   \eqn{n} is the sample size, the length of the input \code{data} to
#'   \code{\link{fitGP}}.  The columns are named \code{sigma[u]} and \code{xi}
#' @seealso \code{\link{Bernoulli}} for maximum likelihood inference for the
#'   Bernoulli distribution.
#' @seealso \code{\link{generalisedPareto}} for maximum likelihood inference
#'   for the generalised Pareto distribution.
#' @examples
#' library(lite)
#' got_exdex <- requireNamespace("exdex", quietly = TRUE)
#' if (got_exdex) {
#'
#'  # estfun.Bernoulli
#'  bfit <- fitBernoulli(exdex::cheeseboro > 45)
#' head(estfun(bfit))
#'
#'  # estfun.generalisedPareto
#'  gpfit <- fitGP(c(exdex::cheeseboro), u = 45)
#' # head(estfun(gpfit))
#'
#' }
#' @name estfun
NULL
## NULL

# Bernoulli

#' @rdname estfun
#' @export
estfun.Bernoulli <- function(x, ...) {
  U <- x$data / x$mle - (1 - x$data) / (1 - x$mle)
  dim(U)<- c(length(U), 1)
  colnames(U) <- names(coef(x))
  return(U)
}

# GP

#' @rdname estfun
#' @param eps,m These arguments control the estimation of the observed
#'   information in \code{gpObsInfo} when the GP shape parameter \eqn{\xi} is
#'   very close to zero.  In these cases, direct calculation is unreliable.
#'   \code{eps} is a (small, positive) numeric scalar.  If the absolute value
#'   of the input value of \eqn{\xi}, that is, \code{pars[2]}, is smaller than
#'   \code{eps} then we approximate the \code{[2, 2]} element using a Taylor
#'   series expansion in \eqn{\xi}, evaluated up to and including the
#'   \code{m}th term.
#' @export
estfun.GP <- function(x, eps = 1e-5, m = 3, ...) {
  if (eps <= 0) {
    stop("'eps' must be positive")
  }
  if (m < 0) {
    stop("'m' must be non-negative")
  }
  pars <- coef(x)
  sigma <- pars[1]
  xi <- pars[2]
  z <- xi / sigma
  # Threshold excesses
  y <- x$exceedances - x$threshold
  zy <- z * y
  t0 <- 1 + zy
  U <- matrix(NA, nrow = length(y), ncol = 2)
  U[, 1] <- -1 / sigma + (xi + 1) * y / (t0 * sigma ^ 2)
  if (abs(xi) < eps) {
    i <- 0:m
    sum_fn <- function(zy) {
      return(sum((-1) ^ i * (i + 1) * zy ^ i / (i + 2)))
    }
    tsum <- vapply(zy, sum_fn, 0.0)
    U[, 2] <- y ^ 2 * tsum / sigma ^ 2 - y / (sigma * t0)
  } else {
    t1 <- log1p(zy) / z ^ 2
    t2 <- y / (z * t0)
    U[, 2] <- (t1 - t2) / sigma ^ 2 - y / (sigma * t0)
  }
  colnames(U) <- names(coef(x))
  return(U)
}
