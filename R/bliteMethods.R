#' Methods for objects of class \code{"blite"}
#'
#' Methods for objects of class \code{"blite"} returned from
#' \code{\link{blite}}. \code{confint.blite} is a misnomer: it returns
#' Bayesian credible intervals.
#' @param x An object inheriting from class \code{"blite"}, a result of a
#'   call to \code{\link{blite}}.
#' @param object An object inheriting from class \code{"blite"}, a result of a
#'   call to \code{\link{blite}}.
#' @param ... For \code{plot.blite}: arguments passed to
#'   \code{\link[graphics:plot.default]{plot}}, such as graphical parameters.
#'
#'   For \code{coef.blite}: additional arguments passed to \code{fun}.
#'
#'   For \code{print.summary.blite}: additional arguments passed to
#'   \code{\link{print.default}}.
#'
#'   Otherwise \code{...} is unused.
#' @return
#'   \code{plot.blite}: No return value, only the plot is produced.
#'
#'   \code{coef.blite}: a numeric vector of length 4 with names
#'     \code{c("p[u]", "sigma[u]", "xi", "theta")}.  The values of summary
#'     statistics calculated using the function \code{fun}.
#'
#'   \code{vcov.blite}: a \eqn{4 \times 4}{4 x 4} matrix with row and
#'     column names \code{c("p[u]", "sigma[u]", "xi", "theta")}.  An estimate
#'     of the posterior covariance matrix, calculated using
#'     \code{\link[stats]{cov}}.
#'
#'   \code{nobs.blite}: a numeric vector of length 3 with names
#'     \code{c("p[u]", "gp", "theta")}.  The respective number of observations
#'     used to infer \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'     (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'     \eqn{\xi}) and \eqn{\theta}.
#'
#'   \code{summary.blite}: an object containing the original function call and
#'     a matrix of summaries of the posterior samples for each of the
#'     parameters.  If \code{short = TRUE} then there are 2 columns, containing
#'     either the sample posterior means and standard deviations
#'     (\code{mean = TRUE}) or the sample posterior medians and inter-quartile
#'     ranges (\code{mean = FALSE}).  If \code{short = FALSE} then there are 4
#'     columns, with each column containing the usual 6-number summary produced
#'     by \code{\link[base]{summary}}. The object is printed by
#'     \code{\link{print.summary.blite}}.
#'
#'   \code{print.summary.blite}: the argument \code{x} is returned, invisibly.
#'
#'   \code{confint.blite}: a numeric matrix with 2 columns giving the lower and
#'     upper credible limits for each parameter. These columns are labelled
#'     as \code{(1-level)/2} and \code{1-(1-level)/2}, expressed as a
#'     percentage, by default \code{2.5\%} and \code{97.5\%}.  The row names
#'     are the names of the parameters supplied in \code{parm}.
#' @seealso \code{\link{blite}} to perform frequentist threshold-based
#'   inference for time series extremes.
#' @seealso \code{\link{predict.blite}}: for predictive inference for the
#'   largest value observed in \eqn{N} years.
#' @name bliteMethods
NULL
## NULL

# ================================= plot.blite =============================== #

#' Plot method for objects of class \code{"blite"}
#'
#' @param which A character scalar indicating which plot(s) to produce.
#'   If \code{which = "all"} then all 4 plots described in \strong{Details}
#'   are produced.  Otherwise, only one of these plots is produced, with the
#'   possible names of the arguments being in the order that the plots are
#'   described in \strong{Details}.
#' @details For \code{plot.blite}, if \code{which = "all"} then 4 plots are produced.
#'     \itemize{
#'       \item{Top left: histogram of the posterior sample for the threshold
#'         exceedance probability
#'         \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}}.}
#'       \item{Top right: scatter plot of posterior sample for the GP
#'         parameters
#'         (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'         \eqn{\xi}).
#'         The linear constraint
#'         \ifelse{html}{\eqn{\xi} > -\eqn{\sigma}\out{<sub>u</sub>} / \eqn{x}
#'         \out{<sub>(n)</sub>}}{\eqn{\xi > -\sigma_u / x_{(n)}}}
#'         is drawn on the plot.}
#'       \item{Bottom left: histogram of the posterior sample for the GP shape
#'         parameter \eqn{\xi}.}
#'       \item{Bottom right: histogram of the posterior sample for the extremal
#'         index \eqn{\theta}.}
#'     }
#' @rdname bliteMethods
#' @export
plot.blite <- function(x, which = c("all", "pu", "gp", "xi", "theta"), ...) {
  if (!inherits(x, "blite")) {
    stop("use only with \"blite\" objects")
  }
  which <- match.arg(which)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(4, 4, 1, 1))
  if (which == "all") {
    which <- c("pu", "gp", "xi", "theta")
    graphics::layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  }
  my_hist <- function(x, ..., xlab = parse(text = colnames(x)[1]),
                      main = "") {
    graphics::hist(x, probability = TRUE, ..., xlab = xlab, main = main)
  }
  my_scatterplot <- function(x, ..., xlab = parse(text = colnames(x)[1]),
                             ylab = parse(text = colnames(x)[2])) {
    graphics::plot(x, ..., xlab = xlab, ylab = ylab)
  }
  # Bernoulli (p[u])
  if ("pu" %in% which) {
    temp <- attr(x, "Bernoulli")
    my_hist(temp$bin_sim_vals, ...)
  }
  # GP (sigma[u], xi)
  if ("gp" %in% which) {
    temp <- attr(x, "gp")
    my_scatterplot(temp$sim_vals, ...)
    ofit <- attr(attr(attr(x, "flite_object"), "gp"), "original_fit")
    m <- max(ofit$exceedances) - ofit$threshold
    graphics::abline(a = 0, b = -1 / m)
  }
  # GP shape parameter xi
  if ("xi" %in% which) {
    temp <- attr(x, "gp")
    my_hist(temp$sim_vals[, "xi", drop = FALSE], ...)
  }
  # K-gaps for theta
  if ("theta" %in% which) {
    temp <- attr(x, "theta")
    my_hist(temp$sim_vals, ...)
  }
  return(invisible())
}

# ================================ coef.blite =============================== #

#' Posterior point estimate summary for objects of class \code{"blite"}
#'
#' @rdname bliteMethods
#' @param fun A summary function to be applied to each column of the simulated
#'   values in \code{object}.  If \code{fun} is missing then
#'   \code{\link[base]{mean}} is used.
#' @export
coef.blite <- function(object, fun, ...) {
  if (!inherits(object, "blite")) {
    stop("use only with \"blite\" objects")
  }
  if (missing(fun)) {
    fun <- mean
  }
  cf <- apply(object, 2, fun, ...)
  cf <- matrix(cf, ncol = 4)
  colnames(cf) <- c("p[u]", "sigma[u]", "xi", "theta")
  return(cf)
}

# ================================ vcov.blite =============================== #

#' Calculate the variance-covariance matrix for an object of class
#' \code{"blite"}
#'
#' @rdname bliteMethods
#' @export
vcov.blite <- function(object, ...) {
  if (!inherits(object, "blite")) {
    stop("use only with \"blite\" objects")
  }
  vc <- stats::cov(object)
  vc_names <- c("p[u]", "sigma[u]", "xi", "theta")
  dim(vc) <- c(length(vc_names), length(vc_names))
  dimnames(vc) <- list(vc_names, vc_names)
  return(vc)
}

# ================================ nobs.blite =============================== #

#' Extract the number of observations from a fit for class \code{"blite"}
#'
#' @rdname bliteMethods
#' @export
nobs.blite <- function(object,  ...) {
  if (!inherits(object, "blite")) {
    stop("use only with \"blite\" objects")
  }
  return(nobs(attr(object, "flite_object")))
}

# =============================== summary.blite ============================== #

#' Summarising times series extreme fits
#'
#' @param short A logical scalar that determines the form of the output. See
#'   \strong{Details}.
#' @param mean A logical scalar.  Determines the form of the output if
#'   \code{short = TRUE}. See \strong{Details}.
#' @param digits An integer. Passed to \code{\link[base:Round]{signif}} to
#'   round the values in the summary.
#' @rdname bliteMethods
#' @export
summary.blite <- function(object, short = TRUE, mean = TRUE,
                          digits = max(3, getOption("digits") - 3L), ...) {
  if (!inherits(object, "blite")) {
    stop("use only with \"blite\" objects")
  }
  res <- attributes(object)["call"]
  if (short) {
    if (mean) {
      posterior_means <- signif(as.vector(coef(object)), digits = digits)
      posterior_sds <- signif(as.vector(sqrt(diag(vcov(object)))),
                              digits = digits)
      res$matrix <- cbind(`Posterior mean` = posterior_means,
                          `Posterior SD` = posterior_sds)
    } else {
      posterior_medians <- signif(as.vector(coef(object, fun = stats::median)),
                                  digits = digits)
      uq <- as.vector(coef(object, fun = stats::quantile, probs = 0.75))
      lq <- as.vector(coef(object, fun = stats::quantile, probs = 0.25))
      posterior_iqrs <- signif(uq - lq, digits = digits)
      res$matrix <- cbind(`Posterior median` = posterior_medians,
                          `Posterior IQR` = posterior_iqrs)
    }
    rownames(res$matrix) <- colnames(coef(object))
  } else {
    temp <- object
    class(temp) <- "matrix"
    res$matrix <- summary(temp)
    colnames(res$matrix) <- colnames(coef(object))
  }
  class(res) <- "summary.blite"
  return(res)
}

# ============================ print.summary.blite =========================== #

#' Print method for objects of class \code{"summary.blite"}
#'
#' @rdname bliteMethods
#' @export
print.summary.blite <- function(x, ...) {
  if (!inherits(x, "summary.blite")) {
    stop("use only with \"summary.blite\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")
  print(x$matrix, ...)
  invisible(x)
}

# ================================ confint.blite ============================ #

#' Credible intervals for \code{"blite"} objects
#'
#' @param object An object of class \code{"blite"}, returned by
#'   \code{\link{blite}}.
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
#' @param level The credible level required.  A numeric scalar in (0, 1).
#' @rdname bliteMethods
#' @export
confint.blite <- function(object, parm = "all", level = 0.95, ...) {
  if (!inherits(object, "blite")) {
    stop("use only with \"blite\" objects")
  }
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
  # Set up a matrix to store the results
  ci_mat <- matrix(NA, nrow = length(parm), ncol = 2)
  # Estimate the equi-tailed interval
  low <- (1 - level) / 2
  # Use only the posterior samples for the required parameters
  parm_values <- c("pu", "sigmau", "xi", "theta")
  colnames(object) <- parm_values
  ci_mat <- t(apply(object[, which(is.element(parm_values, parm)),
                           drop = FALSE], 2, quantile,
                    probs = c(low, 1 - low)))
  return(ci_mat)
}
