#' Inference for the generalised Pareto distribution
#'
#' Functions involved in making inferences about the scale and shape
#' parameters of a generalised Pareto distribution using maximum likelihood
#' estimation.
#'
#' @param data A numeric vector of raw data.  Missing values are removed using
#'   \code{\link[stats:na.fail]{na.omit}}.
#' @param u A numeric scalar.  The extremal value threshold.
#' @param pars A numeric parameter vector of length 2 containing the values of
#'   the generalised Pareto scale and shape parameters.
#' @param excesses A numeric vector of threshold excesses, that is, amounts
#'   by which exceedances of \code{u} exceed \code{u}.
#' @param object A fitted model object returned from \code{fitGP()}.
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}} (if \code{cluster = NULL}),
#'   or \code{\link[sandwich:vcovCL]{meatCL}} (if \code{cluster} is not
#'   \code{NULL}).
#' @details
#' \code{fitGP}: fit a generalised Pareto distribution using maximum likelihood
#'   estimation.  This function calls \code{\link[revdbayes]{grimshaw_gp_mle}}.
#'
#' \code{gpObsInfo}: calculates the observed information matrix for a random
#' sample \code{excesses} from the generalized Pareto distribution, that is,
#' the negated Hessian matrix of the generalized Pareto independence
#' log-likelihood, evaluated at \code{pars}.
#'
#' \code{logLikVector.GP}: calculates contributions to an independence
#' log-likelihood based on the generalised Pareto distribution.
#'
#' \code{gpObsInfo} returns a 2 by 2 matrix with row and columns names
#' \code{c("sigma[u]", "xi")}.
#' \code{nobs, coef, vcov} and \code{logLik} methods are provided.
#' @return
#' \code{fitGP} returns an object of class \code{"GP"}, a list
#' with components: \code{maxLogLik, threshold, mle, vcov, exceedances, nexc},
#' where \code{exceedances} is a vector containing the values that exceed the
#' threshold \code{threshold} and \code{nexc} is the length of this vector.
#'
#' \code{logLikVector.GP} returns an object of class \code{"logLikVector"}, a
#' vector length \code{length(data)} containing the likelihood contributions
#' from the individual observations in \code{data}.
#' @examples
#' # Set up data and set a threshold
#' cdata <- c(exdex::cheeseboro)
#' u <- quantile(cdata, probs = 0.9, na.rm = TRUE)
#'
#' # Fit a generalised Pareto distribution
#' fit <- fitGP(cdata, u)
#'
#' # Calculate the log-likelihood at the MLE
#' res <- logLikVector(fit)
#'
#' # The logLik method sums the individual log-likelihood contributions.
#' logLik(res)
#'
#' # nobs, coef, vcov, logLik methods for objects returned from fitGP()
#' nobs(fit)
#' coef(fit)
#' vcov(fit)
#' logLik(fit)
#' @name generalisedPareto
NULL
## NULL


# ================================= fitGP =================================== #

#' @rdname generalisedPareto
#' @export
fitGP <- function(data, u) {
  if (!is.vector(data)) {
    stop("'data' must be a vector")
  }
  if (length(u) != 1) {
    stop("'u' must have length 1")
  }
  # Remove missing values
  data <- na.omit(data)
  # Create excesses of threshold u
  excesses <- data[data > u] - u
  # Call Grimshaw (1993) function, note: k is -xi, a is sigma
  grimshaw_fit <- revdbayes::grimshaw_gp_mle(excesses)
  res <- list()
  # mle for (sigma, xi)
  res$mle <- c("sigma[u]" = grimshaw_fit$a, xi = -grimshaw_fit$k)
  # number of threshold excesses
  res$nexc <- length(excesses)
  sc <- rep_len(res$mle[1], res$nexc)
  xi <- res$mle[2]
  res$vcov <- solve(gpObsInfo(res$mle, excesses))
  res$maxLogLik <- -sum(log(sc)) -
    sum(log(1 + xi * excesses / sc) * (1 / xi + 1))
  # Return the exceedances (values that lie above the threshold)
  res$exceedances <- data[data > u]
  res$threshold <- u
  class(res) <- "GP"
  return(res)
}

# ============================== gpObsInfo ================================== #

#' @rdname generalisedPareto
#' @export
gpObsInfo <- function(pars, excesses) {
  y <- excesses
  s <- pars[1]
  x <- pars[2]
  i <- matrix(NA, 2, 2)
  i[1, 1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1, 2] <- i[2, 1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  i[2, 2] <- sum(2 * log(1 + x * y / s) / x ^ 3 - 2 * y / (s + x * y) / x ^ 2 -
                  (1 + 1 / x) * y ^ 2 / (s + x * y) ^ 2)
  dimnames(i) <- list(c("sigma[u]", "xi"), c("sigma[u]", "xi"))
  return(i)
}

# Methods for class "GP"

#' @rdname generalisedPareto
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

#' @rdname generalisedPareto
#' @export
nobs.GP <- function(object, ...) {
  return(object$nexc)
}

#' @rdname generalisedPareto
#' @export
coef.GP <- function(object, ...) {
  val <- object$mle
  names(val) <- c("sigma[u]", "xi")
  return(val)
}

#' @rdname generalisedPareto
#' @export
vcov.GP <- function(object, ...) {
  vc <- object$vcov
  dimnames(vc) <- list(c("sigma[u]", "xi"), c("sigma[u]", "xi"))
  return(vc)
}

#' @rdname generalisedPareto
#' @export
logLik.GP <- function(object, ...) {
  val <- object$maxLogLik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- length(coef(object))
  class(val) <- "logLik"
  return(val)
}
