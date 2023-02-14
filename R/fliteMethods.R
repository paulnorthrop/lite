#' Methods for objects of class \code{"flite"}
#'
#' Methods for objects of class \code{"flite"} returned from
#' \code{\link{flite}}.
#' @param x An object inheriting from class \code{"flite"}, a result of a
#'   call to \code{\link{flite}}.
#' @param object An object inheriting from class \code{"flite"}, a result of a
#'   call to \code{\link{flite}}.
#' @param ... For \code{plot.flite}: arguments passed to
#'   \code{\link[graphics:plot.default]{plot}}, such as graphical parameters.
#'
#'   For \code{print.summary.flite}: additional arguments passed to
#'   \code{\link{print.default}}.
#'
#'   For \code{confint.flite}: additional arguments passed to
#'   \code{\link[chandwich]{conf_intervals}}.
#'
#'   Otherwise \code{...} is unused.
#' @return
#'   \code{plot.flite}: No return value, only the plot is produced.
#'
#'   \code{coef.flite}: a numeric vector of length 4 with names
#'     \code{c("p[u]", "sigma[u]", "xi", "theta")}.  The MLEs of the parameters
#'     \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'     \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'     \eqn{\xi} and \eqn{\theta}.
#'
#'   \code{vcov.flite}: a \eqn{4 \times 4}{4 x 4} matrix with row and
#'     column names \code{c("p[u]", "sigma[u]", "xi", "theta")}.  The estimated
#'     variance-covariance matrix for the model parameters.  If
#'     \code{adjust = TRUE} then the elements corresponding to
#'     \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'     \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'     and \eqn{\xi} are adjusted for cluster dependence using
#'     a sandwich estimator; otherwise they are not adjusted.
#'
#'   \code{nobs.flite}: a numeric vector of length 3 with names
#'     \code{c("p[u]", "gp", "theta")}.  The respective number of observations
#'     used to estimate \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'     (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'     \eqn{\xi}) and \eqn{\theta}.
#'
#'   \code{logLik.flite}: an object of class \code{"logLik"}: a numeric scalar
#'     with value equal to the maximised log-likelihood.  This is the sum of
#'     contributions from three fitted models, from a Bernoulli model for
#'     occurrences of threshold exceedances, a generalised Pareto model for
#'     threshold excesses and a \eqn{K}-gaps model for the extremal index.
#'     The returned object also has attributes \code{nobs}, the numbers of
#'     observations used in each of these model fits, and \code{"df"} (degrees
#'     of freedom), which is equal to the number of total number of parameters
#'     estimated (4).
#'
#'   \code{summary.flite}: an object containing the original function call and
#'     a matrix of estimates and estimated standard errors with row names
#'     \code{c("p[u]", "sigma[u]", "xi", "theta")}.  The object is printed by
#'     \code{\link{print.summary.flite}}.
#'
#'   \code{print.summary.flite}: the argument \code{x} is returned, invisibly.
#'
#'   \code{confint.flite}: a numeric matrix with 2 columns giving the lower and
#'     upper confidence limits for each parameter. These columns are labelled
#'     as \code{(1-level)/2} and \code{1-(1-level)/2}, expressed as a
#'     percentage, by default \code{2.5\%} and \code{97.5\%}.  The row names
#'     are the names of the parameters supplied in \code{parm}.
#' @seealso \code{\link{flite}} to perform frequentist threshold-based
#'   inference for time series extremes.
#' @name fliteMethods
NULL
## NULL

# ================================= plot.flite ============================== #

#' Plot method for objects of class \code{"flite"}
#'
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
#'   constraint that
#'   \ifelse{html}{\eqn{\xi} > -\eqn{\sigma}\out{<sub>u</sub>} /
#'   \eqn{x}\out{<sub>(n)</sub>}}{\eqn{\xi > -\sigma_u / x_{(n)}}}
#'   where \ifelse{html}{\eqn{x}\out{<sub>(n)</sub>}}{\eqn{x_{(n)}}}
#'   is the largest excesses of the threshold \eqn{u}, is preserved.
#' @details For \code{plot.flite}, if \code{which = "all"} then 4 plots are
#'   produced.
#'     \itemize{
#'       \item{Top left: (adjusted) log-likelihood for the threshold exceedence
#'         probability \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'         with a horizontal line indicating a 95\% confidence interval for
#'         \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}}.}
#'       \item{Top right: contour plot of the (adjusted) log-likelihood for the
#'         GP parameters
#'         (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'         \eqn{\xi}),
#'         showing (25, 50, 75, 90, 95)\% confidence regions. The linear
#'         constraint
#'         \ifelse{html}{\eqn{\xi} > -\eqn{\sigma}\out{<sub>u</sub>} / \eqn{x}
#'         \out{<sub>(n)</sub>}}{\eqn{\xi > -\sigma_u / x_{(n)}}} is drawn on
#'         the plot.}
#'       \item{Bottom left: (adjusted) log-likelihood for \eqn{\xi}, with a
#'         horizontal line indicating a 95\% confidence interval for
#'         \eqn{\xi}.}
#'       \item{Bottom right: log-likelihood for the extremal index
#'         \eqn{\theta}, with a horizontal line indicating a 95\% confidence
#'         interval for \eqn{\theta}.}
#'     }
#' @rdname fliteMethods
#' @export
plot.flite <- function(x, which = c("all", "pu", "gp", "xi", "theta"),
                       adj_type = c("vertical", "none", "cholesky",
                                    "spectral"),
                       ...) {
  if (!inherits(x, "flite")) {
    stop("use only with \"flite\" objects")
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
    ci <- chandwich::conf_intervals(attr(x, "Bernoulli"), type = adj_type)
    bplot <- function(obj, ..., xlab = expression(p[u]),
                      ylab = "log-likelihood") {
      plot(obj, ..., xlab = xlab, ylab = ylab)
    }
    bplot(ci, ...)
  }
  # GP (sigma[u], xi)
  if ("gp" %in% which) {
    gp <- attr(x, "gp")
    cr <- chandwich::conf_region(gp, type = adj_type)
    gpplot <- function(obj, ..., conf = c(25, 50, 75, 90, 95)) {
      plot(obj, ..., conf = conf)
    }
    gpplot(cr, ...)
    ofit <- attr(gp, "original_fit")
    m <- max(ofit$exceedances) - ofit$threshold
    graphics::abline(a = 0, b = -1 / m)
  }
  # GP shape parameter xi
  if ("xi" %in% which) {
    gp <- attr(x, "gp")
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
    tplot(confint(attr(x, "theta")), ...)
  }
  return(invisible())
}

# ================================ coef.flite =============================== #

#' Extract model coefficients method for objects of class \code{"flite"}
#'
#' @rdname fliteMethods
#' @export
coef.flite <- function(object, ...) {
  if (!inherits(object, "flite")) {
    stop("use only with \"flite\" objects")
  }
  bfit <- attr(object, "Bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "theta")
  cf <- c(attr(bfit, "MLE"), attr(gfit, "MLE"), kfit$theta)
  names(cf) <- c("p[u]", "sigma[u]", "xi", "theta")
  return(cf)
}

# ================================ vcov.flite =============================== #

#' Calculate the variance-covariance matrix for an object of class
#' \code{"flite"}
#'
#' @param adjust A logical scalar.  If \code{adjust = TRUE} then the elements
#'   of the variance-covariance matrix corresponding to
#'   (\ifelse{html}{\eqn{p}\out{<sub>u</sub>},
#'   \eqn{\sigma}\out{<sub>u</sub>}}{\eqn{p_u}, \eqn{\sigma_u}}, \eqn{\xi}),
#'   are estimated using a sandwich estimator. See \code{\link{flite}}.
#'   Otherwise, this matrix is the inverse of the observed information matrix.
#' @rdname fliteMethods
#' @export
vcov.flite <- function(object, adjust = TRUE, ...) {
  if (!inherits(object, "flite")) {
    stop("use only with \"flite\" objects")
  }
  bfit <- attr(object, "Bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "theta")
  vc <- matrix(0, 4, 4)
  if (adjust) {
    vc[1, 1] <- attr(bfit, "adjVC")
    vc[2:3, 2:3] <- attr(gfit, "adjVC")
  } else {
    vc[1, 1] <- attr(bfit, "VC")
    vc[2:3, 2:3] <- attr(gfit, "VC")
  }
  vc[4, 4] <- vcov(kfit)
  vc_names <- c("p[u]", "sigma[u]", "xi", "theta")
  dim(vc) <- c(length(vc_names), length(vc_names))
  dimnames(vc) <- list(vc_names, vc_names)
  return(vc)
}

# ================================ nobs.flite ============================== #

#' Extract the number of observations from a fit for class \code{"flite"}
#'
#' @rdname fliteMethods
#' @export
nobs.flite <- function(object,  ...) {
  if (!inherits(object, "flite")) {
    stop("use only with \"flite\" objects")
  }
  bnobs <- attr(attr(object, "Bernoulli"), "nobs")
  gnobs <- attr(attr(object, "gp"), "nobs")
  knobs <-  attr(object, "theta")$ss$n_kgaps
  n <- c(bnobs, gnobs, knobs)
  names(n) <- c("p[u]", "gp", "theta")
  return(n)
}

# ================================ logLik.flite ============================= #

#' Extract log-likelihood for objects of class \code{"flite"}
#'
#' @rdname fliteMethods
#' @export
logLik.flite <- function(object, ...) {
  if (!inherits(object, "flite")) {
    stop("use only with \"flite\" objects")
  }
  # Note: the value of type does not affect the value of logLik
  val <- object(coef(object))
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- 4
  class(val) <- "logLik"
  return(val)
}

# =============================== summary.flite ============================= #

#' Summarising times series extreme fits
#'
#' @param digits An integer. Passed to \code{\link[base:Round]{signif}} to
#'   round the values in the summary.
#' @rdname fliteMethods
#' @export
summary.flite <- function(object, adjust = TRUE,
                          digits = max(3, getOption("digits") - 3L), ...) {
  if (!inherits(object, "flite")) {
    stop("use only with \"flite\" objects")
  }
  res <- attributes(object)["call"]
  mles <- signif(coef(object), digits = digits)
  ses <- signif(sqrt(diag(vcov(object, adjust = adjust))), digits = digits)
  res$matrix <- cbind(`Estimate` = mles, `Std. Error` = ses)
  rownames(res$matrix) <- names(mles)
  class(res) <- "summary.flite"
  return(res)
}

# ============================ print.summary.flite ========================== #

#' Print method for objects of class \code{"summary.flite"}
#'
#' @rdname fliteMethods
#' @export
print.summary.flite <- function(x, ...) {
  if (!inherits(x, "summary.flite")) {
    stop("use only with \"summary.flite\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")
  print(x$matrix, ...)
  invisible(x)
}

# ================================ confint.flite ============================ #

#' Confidence intervals for \code{"flite"} objects
#'
#' @param object An object of class \code{"flite"}, returned by
#'   \code{\link{flite}}.
#' @param parm A character vector specifying the parameters for which
#'   confidence intervals are required. The default, \code{which = "all"},
#'   produces confidence intervals for all the parameters, that is,
#'   \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'   \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'   \eqn{\xi} and \eqn{\theta}. If \code{which = "gp"} then intervals are
#'   produced only for
#'   \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}} and
#'   \eqn{\xi}. Otherwise, \code{parm} must be a subset of
#'   \code{c("pu", "sigmau", "xi", "theta")}.
#' @param level The confidence level required.  A numeric scalar in (0, 1).
#' @param profile A logical scalar. If \code{TRUE} then confidence intervals
#'   based on an (adjusted) profile loglikelihood are returned.  If
#'   \code{FALSE} then intervals based on approximate large sample normal
#'   theory, which are symmetric about the MLE, are returned.
#' @rdname fliteMethods
#' @export
confint.flite <- function(object, parm = "all", level = 0.95,
                          adj_type = c("vertical", "none", "cholesky",
                                       "spectral"), profile = TRUE, ...) {
  if (!inherits(object, "flite")) {
    stop("use only with \"flite\" objects")
  }
  parm_values <- c("pu", "sigmau", "xi", "theta")
  check_values <- c("pu", "sigmau", "xi", "theta", "all", "gp")
  p_message <- "c(''pu'', ''sigmau'', ''xi'', ''theta'')"
  if (!all(is.element(parm, check_values))) {
    stop(paste("''parm'' must be ''all'', ''gp'' or a subset of", p_message))
  }
  if (length(parm) == 1) {
    if (parm == "all") {
      parm <- c("pu", "sigmau", "xi", "theta")
    } else if (parm == "gp") {
      parm <- c("sigmau", "xi")
    }
  }
  if (level <= 0 | level >= 1) {
    stop("''level'' must be in (0, 1)")
  }
  adj_type <- match.arg(adj_type)
  # Set up a matrix to store the results
  ci_mat <- matrix(NA, nrow = length(parm), ncol = 2)
  # Set a counter to keep track of the row in which to store results
  therow <- 1L
  # Bernoulli (p[u])
  if ("pu" %in% parm) {
    ci_pu <- chandwich::conf_intervals(attr(object, "Bernoulli"),
                                       conf = 100 * level, type = adj_type,
                                       profile = profile, ...)
    if (profile) {
      ci_mat[therow, ] <- ci_pu$prof_CI
    } else {
      ci_mat[therow, ] <- ci_pu$sym_CI
    }
    therow <- therow + 1L
  }
  # GP scale parameter sigma[u]
  if ("sigmau" %in% parm) {
    gp <- attr(object, "gp")
    ci_sigmau <- chandwich::conf_intervals(gp, which_pars = "sigma[u]",
                                           conf = 100 * level, type = adj_type,
                                           profile = profile, ...)
    if (profile) {
      ci_mat[therow, ] <- ci_sigmau$prof_CI
    } else {
      ci_mat[therow, ] <- ci_sigmau$sym_CI
    }
    therow <- therow + 1L
  }
  # GP shape parameter xi
  if ("xi" %in% parm) {
    gp <- attr(object, "gp")
    ci_xi <- chandwich::conf_intervals(gp, which_pars = "xi",
                                       conf = 100 * level, type = adj_type,
                                       profile = profile)
    if (profile) {
      ci_mat[therow, ] <- ci_xi$prof_CI
    } else {
      ci_mat[therow, ] <- ci_xi$sym_CI
    }
    therow <- therow + 1L
  }
  # K-gaps for theta
  if ("theta" %in% parm) {
    if (profile) {
      interval_type <- "lik"
    } else {
      interval_type <- "norm"
    }
    ci_theta <- confint(attr(object, "theta"), level = level,
                        interval_type = interval_type)
    ci_mat[therow, ] <- ci_theta$cis
  }
  low <- paste0(100 * (1 - level)/ 2, "%")
  up <- paste0(100 - 100 * (1 - level)/ 2, "%")
  colnames(ci_mat) <- c(low, up)
  rownames(ci_mat) <- parm
  return(ci_mat)
}
