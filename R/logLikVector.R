#' Functions for log-likelihood contributions
#'
#' Generic function for calculating log-likelihood contributions from
#' individual observations for a fitted model object.
#'
#' @param object A fitted model object.
#' @param pars A numeric parameter vector.
#'
#'   For \code{logLikVector.Bernoulli} this is a vector of length 1 containing
#'   a value of the Bernoulli success probability.
#'
#'   For \code{logLikVector.GP} this is a numeric vector of length 2 containing
#'   the values of the generalised Pareto scale (\eqn{\sigma_u}) and shape
#'   (\eqn{\xi}) parameters.
#' @param ... Further arguments. None are currently used for
#'   \code{logLikVector.Bernoulli} or \code{logLikVector.GP}.
#' @details A \code{logLikVector} method is used to construct a log-likelihood
#'   function to supply as the argument \code{loglik} to the function
#'   \code{\link[chandwich]{adjust_loglik}}, which performs log-likelihood
#'   adjustment for parts 1 and 2 of the inferences performed by
#'   \code{\link{flite}}.
#'
#'   The \code{logLik} method \code{logLik.LogLikVector} sums the
#'   log-likelihood contributions from individual observations.
#' @return For \code{logLikVector}: an object of class \code{logLikVec}.
#'   This is a numeric vector of length \eqn{n} containing contributions to the
#'   the independence log-likelihood from \eqn{n} observations, with attributes
#'   \code{"df"} (degrees of freedom), giving the number of estimated
#'   parameters in the model, and \code{"nobs"}, giving the number observations
#'   used to perform the estimation.
#'
#'   For \code{logLik.logLikVector}: an object of class \code{logLik}.  This is
#'   a number with the attributes \code{"df"} and \code{"nobs"} as described
#'   above.
#' @seealso \code{\link{Bernoulli}} for maximum likelihood inference for the
#'   Bernoulli distribution.
#' @seealso \code{\link{generalisedPareto}} for maximum likelihood inference
#'   for the generalised Pareto distribution.
#' @examples
#' # logLikVector.Bernoulli
#' bfit <- fitBernoulli(c(exdex::cheeseboro) > 45)
#' bvec <- logLikVector(bfit)
#' head(bvec)
#' logLik(bvec)
#' logLik(bfit)
#'
#' # estfun.generalisedPareto
#' gpfit <- fitGP(c(exdex::cheeseboro), u = 45)
#' gpvec <- logLikVector(gpfit)
#' head(gpvec)
#' logLik(gpvec)
#' logLik(gpfit)
#' @name logLikVector
NULL
## NULL

#' @rdname logLikVector
#' @export
logLikVector <- function(object, ...) {
  UseMethod("logLikVector")
}

#' @rdname logLikVector
#' @export
logLikVector.Bernoulli <- function(object, pars = NULL, ...) {
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  prob <- pars[1]
  if (prob < 0 || prob > 1) {
    val <- -Inf
  } else {
    val <- stats::dbinom(object$obs_data, 1, prob, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVector"
  return(val)
}

#' @rdname logLikVector
#' @export
logLikVector.GP <- function(object, pars = NULL, ...) {
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  sigma <- pars[1]
  xi <- pars[2]
  # Calculate the log-likelihood contributions
  if (any(sigma <= 0)) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(x = object$exceedances, loc = object$threshold,
                          scale = sigma, shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVector"
  return(val)
}

#' @rdname logLikVector
#' @export
logLik.logLikVector <- function(object, ...) {
  save_attributes <- attributes(object)
  object <- sum(object)
  attributes(object) <- save_attributes
  class(object) <- "logLik"
  return(object)
}
