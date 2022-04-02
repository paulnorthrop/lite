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
#'   \code{c("p[u]", "sigma[u]", "xi", "theta")}.
#' @export
coef.lite <- function(object, ...) {
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  bfit <- attr(object, "bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "kgaps")
  cf <- c(attr(bfit, "MLE"), attr(gfit, "MLE"), kfit$theta)
  names(cf) <- c("p[u]", "sigma[u]", "xi", "theta")
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
#'   \code{c("p[u]", "sigma[u]", "xi", "theta")}.
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
  vc_names <- c("p[u]", "sigma[u]", "xi", "theta")
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
#'   \code{c("p[u]", "gp", "theta")}.
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
  names(n) <- c("p[u]", "gp", "theta")
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
  # Note: the value of type does not affect the value of logLik
  val <- object(coef(object))
  attr(val, "nobs") <- nobs(object)
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
#'   \code{c("p[u]", "sigma[u]", "xi", "theta")}.  The object is printed by
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

#' Print method for objects of class \code{"summary.lite"}
#'
#' \code{print} method for an object \code{x} of class \code{"summary.lite"}.
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


# ================================= plot.lite =============================== #

#' Plot method for objects of class \code{"lite"}
#'
#' \code{print} method for an object \code{x} of class \code{"lite"}.
#'
#' @param x An object of class "lite", a result of a call to
#'   \code{\link{flite}}.
#' @param which A character scalar indicating which plot(s) to produce.
#'   If \code{which = "all"} then all 4 plots described in \strong{Details}
#'   are produced.  Otherwise, only one of these plots is produced, with the
#'   possible names of the arguments being in the order that the plots are
#'   described in \strong{Details}.
#' @param adj_type A character scalar passed to
#'   \code{\link[chandwich]{conf_intervals}} and
#'   \code{\link[chandwich]{conf_region}} as the argument \code{type} to select
#'   the type of adjustment applied to the independence log-likelihood.  Of the
#'   3 adjustments, \code{"vertical"} is preferred because it preserves
#'   constraints on the parameters, whereas the \code{"cholesky"} and
#'   \code{"spectral"} adjustment do not.  In the generalised Pareto case the
#'   constraint that \eqn{\xi > - \sigma_u / x_{(n)}}{\xi > \sigma_u / x_(n)},
#'   where \eqn{x_{(n)}}{x_(n)} is the largest excesses of the threshold \eqn{u},
#'   is preserved.
#' @param ... Arguments to be passed to \code{\link[stats]{plot}}, such as
#'   graphical parameters.
#' @details If \code{which = "all"} then 4 plots are produced.
#'     \itemize{
#'       \item{Top left: (adjusted) log-likelihood for the threshold exceedence
#'         probability \eqn{p_u}, with a horizontal line indicating a
#'         95\% confidence interval for \eqn{p_u}.}
#'       \item{Top right: contour plot of the (adjusted) log-likelihood for the
#'         GP parameters \eqn{(\sigma_u, \xi)}, showing
#'         (25, 50, 75, 90, 95)\% confidence regions. The linear constraint
#'         \eqn{\xi > - \sigma_u / x_{(n)}}{\xi > \sigma_u / x_(n)} is drawn
#'         on the plot.}
#'       \item{Bottom left: (adjusted) log-likelihood for \eqn{\xi}, with a
#'         horizontal line indicating a 95\% confidence interval for \eqn{\xi}.}
#'       \item{Bottom right: log-likelihood for the extremal index \eqn{\theta},
#'         with a horizontal line indicating a 95\% confidence interval for
#'         \eqn{\theta}.}
#'     }
#' @section Examples:
#' See the examples in \code{\link{flite}}.
#' @export
plot.lite <- function(x, which = c("all", "pu", "gp", "xi", "theta"),
                      adj_type = c("vertical", "none", "cholesky", "spectral"),
                      ...) {
  if (!inherits(x, "lite")) {
    stop("use only with \"lite\" objects")
  }
  adj_type <- match.arg(adj_type)
  which <- match.arg(which)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(4, 4, 1, 1))
  if (which == "all") {
    which <- c("pu", "gp", "xi", "theta")
    graphics::layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  }
  # Bernoulli (p[u])
  if ("pu" %in% which) {
    ci <- chandwich::conf_intervals(attr(x, "bernoulli"), type = adj_type)
    bplot <- function(obj, ..., xlab = expression(p[u]),
                      ylab = "log-likelihood") {
      plot(obj, ..., xlab = xlab, ylab = ylab)
    }
    bplot(ci, ...)
  }
  # GP - to do, perhaps plot.confreg
  if ("gp" %in% which) {
    gp <- attr(x, "gp")
    cr <- chandwich::conf_region(gp, type = adj_type)
    gpplot <- function(obj, ..., conf = c(25, 50, 75, 90, 95)) {
      plot(obj, ..., conf = conf)
    }
    gpplot(cr, ...)
    ofit <- attr(gp, "original_fit")
    m <- max(ofit$exceedances) - ofit$threshold
    abline(a = 0, b = -1 / m)
  }
  # GP shape parameter xi
  if ("xi" %in% which) {
    ci_xi <- chandwich::conf_intervals(gp, which_pars = "xi", type = adj_type)
    xiplot <- function(obj, ..., ylab = "profile log-likelihood") {
      plot(obj, ..., ylab = ylab)
    }
    xiplot(ci_xi, ...)
  }
  # K-gaps for theta
  if ("theta" %in% which) {
    tplot <- function(obj, ..., main = "") {
      plot(obj, ..., main = main)
    }
    tplot(confint(attr(x, "kgaps")), ...)
  }
  return(invisible())
}

