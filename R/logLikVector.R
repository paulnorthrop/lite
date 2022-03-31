#' Evaluate log-likelihood contributions from specific observations
#'
#' Generic function for calculating log-likelihood contributions from
#' individual observations for a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Further arguments.
#' @export
logLikVector <- function(object, ...) {
  UseMethod("logLikVector")
}

#' Sum log-likelihood contributions from individual observations
#'
#' S3 logLik method for logLikVector objects
#'
#' @param object An object of class \code{"logLikVector"} return from a
#'   \code{logLikVector} method.
#' @param ... Further arguments.
#' @export
logLik.logLikVector <- function(object, ...) {
  save_attributes <- attributes(object)
  object <- sum(object)
  attributes(object) <- save_attributes
  class(object) <- "logLik"
  return(object)
}
