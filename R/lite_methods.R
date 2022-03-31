# ================================ coef.lite =============================== #

#' Extract model coefficients method for objects of class \code{"lite"}
#'
#' \code{coef} method for class \code{"lite"}.
#'
#' @param object an object inheriting from class \code{"lite"}, a result of a
#'   call to \code{\link{flite}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return A numeric vector of length 4 with names
#'   \code{c("p_u", "sigma_u", "xi", "theta")}.
#' @export
coef.lite <- function(object, ...) {
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  bfit <- attr(object, "bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "kgaps")
  cf <- c(attr(bfit, "MLE"), attr(gfit, "MLE"), kfit$theta)
  names(cf) <- c("p_u", "sigma_u", "xi", "theta")
  return(cf)
}

# ================================ vcov.lite =============================== #

#' Calculate the variance-covariance matrix for an object of class
#' \code{"lite"}
#'
#' \code{vcov} method for class \code{"lite"}.
#'
#' @param object an object inheriting from class \code{"lite"}, a result of a
#'   call to \code{\link{flite}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return A 4 by 4 numeric matrix with row and column names
#'   \code{c("p_u", "sigma_u", "xi", "theta")}.
#' @export
vcov.lite <- function(object, ...) {
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  bfit <- attr(object, "bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "kgaps")
  vc <- matrix(0, 4, 4)
  vc[1, 1] <- attr(bfit, "adjVC")
  vc[2:3, 2:3] <- attr(gfit, "adjVC")
  vc[4, 4] <- vcov(kfit)
  vc_names <- c("p_u", "sigma_u", "xi", "theta")
  dim(vc) <- c(length(vc_names), length(vc_names))
  dimnames(vc) <- list(vc_names, vc_names)
  return(vc)
}

# ================================ nobs.lite =============================== #

#' Extract the number of observations from a fit for class \code{"lite"}
#'
#' \code{nobs} method for class \code{"lite"}.
#'
#' @param object an object inheriting from class \code{"lite"}, a result of a
#'   call to \code{\link{flite}}.
#' @param model One of \code{c("gp", "kgaps", "bernoulli")}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return A numeric vector of length 3 with names
#'   \code{c("p_u", "gp", "theta")}.
#'
#' @export
nobs.lite <- function(object,  ...) {
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  bnobs <- attr(attr(object, "bernoulli"), "nobs")
  gnobs <- attr(attr(object, "gp"), "nobs")
  knobs <-  attr(object, "kgaps")$ss$n_kgaps
  n <- c(bnobs, gnobs, knobs)
  names(n) <- c("p_u", "gp", "theta")
  return(n)
}

# ================================ logLik.lite ============================== #

#' Extract log-likelihood for objects of class \code{"lite"}
#'
#' \code{logLik} method for class \code{"lite"}.
#'
#' @param object an object of class \code{"lite"}, a result of a call to
#'   \code{\link{flite}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return An object of class \code{"logLik"}: a numeric scalar with
#' value equal to the maximised log-likelihood.  This is the sum of
#' contributions from three fitted models, from a Bernoulli model for
#' occurrences of threshold exceedances, a generalised Pareto model for
#' threshold excesses and a \eqn{K}-gaps model for the extremal index.
#' The returned object also has attributes \code{nobs}, the numbers of
#' observations used in each of these model fits, and \code{"df"}, which is
#' equal to the number of total number of parameters estimated (4).
#' @export
logLik.lite <- function(object, ...) {
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  bfit <- attr(object, "bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "kgaps")
  bloglik <- attr(bfit, "max_loglik")
  gloglik <- attr(gfit, "max_loglik")
  kloglik <- kfit$max_loglik
  val <- bloglik + gloglik + kloglik
  attr(val, "nobs") <- nobs(y)
  attr(val, "df") <- 4
  class(val) <- "logLik"
  return(val)
}

# =============================== summary.lite ============================== #

#' Summarising times series extreme fits
#'
#' \code{summary} method for class \code{"lite"}
#'
#' @param object an object inheriting from class \code{"lite"}, a result of a
#'   call to \code{\link{flite}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return An object containing the original function call and a matrix of
#'   estimates and estimated standard errors with row names
#'   \code{c("p_u", "sigma_u", "xi", "theta")}.  The object is printed by
#'   \code{\link{print.summary.lite}}.
#' @section Examples:
#' See the examples in \code{\link{flite}}.
#' @export
summary.lite <- function(object, digits = max(3, getOption("digits") - 3L),
                          ...) {
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  res <- attributes(object)["call"]
  mles <- signif(coef(object), digits = digits)
  ses <- signif(sqrt(diag(vcov(object))), digits = digits)
  res$matrix <- cbind(`Estimate` = mles, `Std. Error` = ses)
  rownames(res$matrix) <- names(mles)
  class(res) <- "summary.lite"
  return(res)
}

# ============================ print.summary.lite =========================== #

#' Print method for objects of class \code{"summary.kgaps"}
#'
#' \code{print} method for an object \code{x} of class \code{"summary.kgaps"}.
#'
#' @param x An object of class "summary.lite", a result of a call to
#'   \code{\link{summary.lite}}.
#' @param ... Additional arguments passed on to \code{\link{print.default}}.
#' @return Prints the numeric matrix returned from \code{\link{summary.lite}}.
#' @section Examples:
#' See the examples in \code{\link{flite}}.
#' @export
print.summary.lite <- function(x, ...) {
  if (!inherits(x, "summary.lite")) {
    stop("use only with \"summary.lite\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")
  print(x$matrix, ...)
  invisible(x)
}

