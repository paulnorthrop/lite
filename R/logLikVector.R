#' Functions for log-likelihoood contributions
#'
#' @name logLikVector
NULL
## NULL

#' Evaluate log-likelihood contributions from specific observations
#'
#' Generic function to calculate log-likelihood contributions from
#' individual observations for a fitted model object.
#'
#' @param x A fitted model object.
#' @param ... Further arguments.
#' @rdname logLikVector
#' @export
logLikVector <- function(x, ...) {
  UseMethod("logLikVector")
}

#' Sum log-likelihood contributions from individual observations
#'
#' @param object An object of class \code{"logLikVector"} returned from a
#'   \code{logLikVector} method.
#' @param ... Further arguments.
#' @rdname logLikVector
#' @export
logLik.logLikVector <- function(object, ...) {
  save_attributes <- attributes(object)
  object <- sum(object)
  attributes(object) <- save_attributes
  class(object) <- "logLik"
  return(object)
}
