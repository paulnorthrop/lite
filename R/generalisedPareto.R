#' Frequentist inference for the generalised Pareto distribution
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
#' @param eps,m These arguments control the estimation of the observed
#'   information in \code{gpObsInfo} when the GP shape parameter \eqn{\xi} is
#'   very close to zero.  In these cases, direct calculation is unreliable.
#'   \code{eps} is a (small, positive) numeric scalar.  If the absolute value
#'   of the input value of \eqn{\xi}, that is, \code{pars[2]}, is smaller than
#'   \code{eps} then we approximate the \code{[2, 2]} element using a Taylor
#'   series expansion in \eqn{\xi}, evaluated up to and including the
#'   \code{m}th term.
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}} (if \code{cluster = NULL}),
#'   or \code{\link[sandwich:vcovCL]{meatCL}} (if \code{cluster} is not
#'   \code{NULL}).
#' @details
#' \code{fitGP}: fit a generalised Pareto distribution using maximum likelihood
#'   estimation, using an \strong{independence} log-likelihood formed by
#'   summing contributions from individual observations. No adjustment for
#'   cluster dependence has been made in estimating the variance-covariance
#'   matrix stored as component in \code{vcov} in the returned object. This
#'   function calls \code{\link[revdbayes]{grimshaw_gp_mle}}.
#'
#' \code{coef}, \code{vcov}, \code{nobs} and \code{logLik} methods are
#' provided for objects of class \code{"GP"} returned from \code{fitGP}.
#'
#' \code{gpObsInfo}: calculates the observed information matrix for a random
#' sample \code{excesses} from the generalized Pareto distribution, that is,
#' the negated Hessian matrix of the generalized Pareto independence
#' log-likelihood, evaluated at \code{pars}.
#' @return
#' \code{fitGP} returns an object of class \code{"GP"}, a list
#' with components: \code{maxLogLik}, \code{threshold}, \code{mle},
#' \code{vcov}, \code{exceedances}, \code{nexc},
#' where \code{exceedances} is a vector containing the values that exceed the
#' threshold \code{threshold} and \code{nexc} is the length of this vector.
#'
#'   \code{coef.GP}: a numeric vector of length 2 with names
#'     \code{c("sigma[u]", "xi")}.  The MLEs of the GP parameters
#'     \eqn{\sigma_u} and \eqn{\xi}.
#'
#'   \code{vcov.GP}: a 2 by 2 numeric matrix with row and column names
#'     \code{c("sigma[u]", "xi")}.  The estimated variance-covariance matrix
#'     for the model parameters. No adjustment for cluster dependence has been
#'     made.
#'
#'   \code{nobs.GP}: a numeric vector of length 1.  The number of
#'     observations used to estimate (\eqn{\sigma_u}, \eqn{\xi}).
#'
#'   \code{logLik.GP}: an object of class \code{"logLik"}: a numeric scalar
#'     with value equal to the maximised log-likelihood. The returned object
#'     also has attributes \code{nobs}, the numbers of observations used in
#'     each of these model fits, and \code{"df"} (degrees of freedom), which is
#'     equal to the number of total number of parameters estimated (2).
#'
#' \code{gpObsInfo} returns a 2 by 2 matrix with row and columns names
#' \code{c("sigma[u]", "xi")}.
#' @examples
#' # Set up data and set a threshold
#' cdata <- c(exdex::cheeseboro)
#'
#' # Fit a generalised Pareto distribution
#' fit <- fitGP(cdata, 45)
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
  data <- stats::na.omit(data)
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
    sum(log1p(xi * excesses / sc) * (1 / xi + 1))
  # Return the exceedances (values that lie above the threshold)
  res$exceedances <- data[data > u]
  res$threshold <- u
  class(res) <- "GP"
  return(res)
}

# Methods for class "GP"

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
nobs.GP <- function(object, ...) {
  return(object$nexc)
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

# ============================== gpObsInfo ================================== #

#' @rdname generalisedPareto
#' @export
gpObsInfo <- function(pars, excesses, eps = 1e-5, m = 3) {
  if (eps <= 0) {
    stop("'eps' must be positive")
  }
  if (m < 0) {
    stop("'m' must be non-negative")
  }
  y <- excesses
  # sigma
  s <- pars[1]
  # xi
  x <- pars[2]
  i <- matrix(NA, 2, 2)
  i[1, 1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1, 2] <- i[2, 1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  # Direct calculation of i22 is unreliable for x close to zero.
  # If abs(x) < eps then we expand the problematic terms (all but t4 below)
  # in powers of z up to z ^ 2. The terms in 1/z and 1/z^2 cancel leaving
  # only a quadratic in z.
  z <- x / s
  zy <- z * y
  t0 <- 1 + zy
  t4 <- y ^ 2 / t0 ^ 2
  if (any(t0 <= 0)) {
    stop("The log-likelihood is 0 for this combination of data and parameters")
  }
  if (abs(x) < eps) {
    j <- 0:m
    sum_fn <- function(zy) {
      return(sum((-1) ^ j * (j ^ 2 + 3 * j + 2) * zy ^ j / (j + 3)))
    }
    tsum <- vapply(zy, sum_fn, 0.0)
    i[2, 2] <- sum(y ^ 3 * tsum / s ^ 3 - t4 / s ^ 2)
  } else {
    t1 <- 2 * log1p(zy) / z ^ 3
    t2 <- 2 * y / (z ^ 2 * t0)
    t3 <- y ^ 2 / (z * t0 ^ 2)
    i[2, 2] <- sum((t1 - t2 - t3) / s ^ 3 - t4 / s ^ 2)
  }
  dimnames(i) <- list(c("sigma[u]", "xi"), c("sigma[u]", "xi"))
  return(i)
}
