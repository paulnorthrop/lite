# Methods for class "GP"

#' @export
logLikVector.GP <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  #
  # Threshold exceedances (values that lie above the threshold)
  response_data <- object$data
  sigma <- pars[1]
  xi <- pars[2]
  # Calculate the log-likelihood contributions
  if (any(sigma <= 0)) {
    val <- -Inf
  } else {
    val <- revdbayes::dgp(response_data, loc = object$threshold, scale = sigma,
                          shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVector"
  return(val)
}

#' @export
nobs.GP <- function(object, ...) {
  return(object$nexc)
}

#' @export
coef.GP <- function(object, ...) {
  val <- object$mle
  return(val)
}

#' @export
vcov.GP <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(c("sigma_u", "xi"), c("sigma_u", "xi"))
  return(vc)
}

#' @export
logLik.GP <- function(object, ...) {
  val <- -object$nllh
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
