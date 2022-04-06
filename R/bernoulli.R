#' Frequentist inference for the Bernoulli distribution
#'
#' Functions involved in making inferences about the probability of success
#' in a Bernoulli distribution using maximum likelihood estimation.
#'
#' @param data A numeric vector of outcomes from Bernoulli trials: 0 for a
#'    failure, 1 for a success.  Alternatively, a logical vector with FALSE
#'    for a failure and TRUE for a success. Missing values are removed using
#'   \code{\link[stats:na.fail]{na.omit}}.
#' @param object A fitted model object returned from \code{fitBernoulli()}.
#' @param ... Further arguments. None are used currently.
#' @details
#' \code{fitBernoulli}: fit a Bernoulli distribution using maximum likelihood
#'   estimation.
#'
#' \code{nobs, coef, vcov} and \code{logLik} methods are provided.
#' @return
#' \code{fitBernoulli} returns an object of class \code{"Bernoulli"}, a list
#' with components: \code{maxLogLik}, \code{mle}, \code{nobs}, \code{vcov},
#' \code{n0}, \code{n1}, \code{data}, \code{obs_data}, where \code{data} are
#' the input data and, \code{obs_data} are the input data after any missing
#' values have been removed, using \code{\link[stats:na.fail]{na.omit}} and
#' \code{n0} and \code{n1} are, respectively, the number of failures and the
#' number of successes.
#' @examples
#' got_exdex <- requireNamespace("exdex", quietly = TRUE)
#' if (got_exdex) {
#'  # Set up data
#'  cdata <- c(exdex::cheeseboro)
#'  u <- 45
#'  exc <- cdata > u
#'
#'  # Fit a Bernoulli distribution
#'  fit <- fitBernoulli(exc)
#'
#'  # Calculate the log-likelihood at the MLE
#'  res <- logLikVector(fit)
#'
#'  # The logLik method sums the individual log-likelihood contributions.
#'  logLik(res)
#'
#'  # nobs, coef, vcov, logLik methods for objects returned from fitBernoulli()
#'  nobs(fit)
#'  coef(fit)
#'  vcov(fit)
#'  logLik(fit)
#' }
#' @name Bernoulli
NULL
## NULL

#' @rdname Bernoulli
#' @export
fitBernoulli <- function(data) {
  if (!is.vector(data)) {
    stop("'data' must be a vector")
  }
  res <- list()
  res$data <- data
  obs_data <- stats::na.omit(data)
  res$obs_data <- obs_data
  res$mle <- mean(obs_data)
  res$nobs <- length(obs_data)
  res$vcov <- as.matrix(res$mle * (1 - res$mle) / res$nobs)
  n1 <- sum(obs_data)
  n0 <- res$nobs - n1
  res$maxLogLik <- n1 * log(res$mle) + n0 * log(1 - res$mle)
  res$n0 <- n0
  res$n1 <- n1
  class(res) <- "Bernoulli"
  return(res)
}

# Methods for class "Bernoulli"

#' @rdname Bernoulli
#' @export
nobs.Bernoulli <- function(object, ...) {
  return(object$nobs)
}

#' @rdname Bernoulli
#' @export
coef.Bernoulli <- function(object, ...) {
  val <- object$mle
  names(val) <- "prob"
  return(val)
}

#' @rdname Bernoulli
#' @export
vcov.Bernoulli <- function(object, ...) {
  vc <- object$vcov
  dimnames(vc) <- list("prob", "prob")
  return(vc)
}

#' @rdname Bernoulli
#' @export
logLik.Bernoulli <- function(object, ...) {
  val <- object$maxLogLik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
