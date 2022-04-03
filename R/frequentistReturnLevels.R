#' Frequentist threshold-based inference for return levels
#'
#' Calculates point estimates and confidence intervals for \code{m}-year
#' return levels for stationary time series fitted extreme value model objects
#' returned from \code{\link{flite}}.  Two types of interval may be returned:
#' (a) intervals based on approximate large-sample normality of the maximum
#' likelihood estimator for return level, which are symmetric about the point
#' estimate, and (b) profile likelihood-based intervals based on an (adjusted)
#' log-likelihood.
#'
#' @param x An object inheriting from class \code{"flite"} returned from
#'   \code{\link{flite}}.
#' @param m A numeric scalar.  The return period, in years.
#' @param level A numeric scalar in (0, 1).  The confidence level required for
#'   confidence interval for the \code{m}-year return level.
#' @param npy A numeric scalar.  The (mean) number of observations per year.
#'   \strong{Setting this appropriately is important}. See \strong{Details}.
#' @param prof A logical scalar.  Should we calculate intervals based on
#'   profile log-likelihood?
#' @param inc A numeric scalar. Only relevant if \code{prof = TRUE}. The
#'   increment in return level by which we move upwards and downwards from the
#'   MLE for the return level in the search for the lower and upper confidence
#'   limits.  If this is not supplied then \code{inc} is set to one hundredth
#'   of the length of the symmetric confidence interval for return level.
#' @param type A character scalar.  The argument \code{type} to the function
#'   returned by \code{\link[chandwich]{adjust_loglik}}, that is, the type of
#'   adjustment made to the independence log-likelihood function in creating
#'   an adjusted log-likelihood function.  See \strong{Details} and
#'   \strong{Value} in \code{\link[chandwich]{adjust_loglik}}.
#' @details \code{npy} provides information about the (intended) frequency of
#'   sampling in time, that is, the number of observations that would be
#'   observed in a year if there are no missing values.  If the number of
#'   observations may vary between years then \code{npy} should be set equal to
#'   the mean number of observations per year.
#'
#'   \strong{Supplying \code{npy}.}
#'   The value of \code{npy} may have been set in the call to
#'   \code{\link{flite}}.  If \code{npy} is supplied by the user in the call to
#'   \code{returnLevel} then this will be used in preference to the value
#'   stored in the fitted model object.  If these two values differ then no
#'   warning will be given.
#'
#'   For details of the definition and estimation of return levels see the
#'   Inference for return levels vignette.
#'
#'   The profile likelihood-based intervals are calculated by
#'   reparameterising in terms of the \code{m}-year return level and estimating
#'   the values at which the (adjusted) profile log-likelihood reaches
#'   the critical value \code{logLik(x) - 0.5 * stats::qchisq(level, 1)}.
#'   This is achieved by calculating the profile log-likelihood for a sequence
#'   of values of this return level as governed by \code{inc}. Once the profile
#'   log-likelihood drops below the critical value the lower and upper limits are
#'   estimated by interpolating linearly between the cases lying either side of
#'   the critical value. The smaller \code{inc} the more accurate (but slower)
#'   the calculation will be.
#' @return A object (a list) of class \code{"retlev", "lax"} with the
#'   components
#'   \item{rl_sym,rl_prof }{Named numeric vectors containing the respective
#'     lower 100\code{level}\% limit, the MLE and the upper
#'     100\code{level}\% limit for the return level.
#'     If \code{prof = FALSE} then \code{rl_prof} will be missing.}
#'   \item{rl_se }{Estimated standard error of the return level.}
#'   \item{max_loglik,crit,for_plot }{If \code{prof = TRUE} then
#'     these components will be present, containing respectively: the maximised
#'     log-likelihood; the critical value and a matrix with return levels in
#'     the first column (\code{ret_levs}) and the corresponding values of the
#'     (adjusted) profile log-likelihood (\code{prof_loglik}).}
#'   \item{m,level }{The input values of \code{m} and \code{level}.}
#'   \item{call }{The call to \code{return_level}.}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \doi{10.1007/978-1-4471-3675-0_3}
#' @seealso \code{\link{plot.retlev}} for plotting the profile log-likelihood
#'   for a return level.
#' @examples
#' ### Cheeseboro wind gusts
#'
#' # Make inferences
#' cdata <- exdex::cheeseboro
#' # Each column of the matrix cdata corresponds to data from a different year
#' # flite() sets cluster automatically to correspond to column (year)
#' cfit <- flite(cdata, u = 45, k = 3)
#'
#' # These data are hourly for one month (January) year so npy = 31 * 24
#' # Large inc set here for speed, sacrificing accuracy
#' rl <- returnLevel(cfit, inc = 2.5, npy = 31 * 24)
#' summary(rl)
#' rl
#' plot(rl)
#'
#' ###
#' @export
returnLevel <- function(x, m = 100, level = 0.95, npy, prof = TRUE,
                         inc = NULL,
                         type = c("vertical", "cholesky", "spectral", "none")) {
  if (!inherits(x, "flite")) {
    stop("use only with \"flite\" objects")
  }
  Call <- match.call(expand.dots = TRUE)
  type <- match.arg(type)
  # Check whether npy is supplied in the call to return_level
  npy_given <- ifelse(missing(npy), FALSE, TRUE)
  # Make inferences about return levels
  temp <- return_level_bingp(x, m, level, npy, prof, inc, type, npy_given)
  temp$m <- m
  temp$level <- level
  temp$call <- Call
  class(temp) <- c("retlev", "lax")
  return(temp)
}

# ------------------------------- plot.retlev ------------------------------- #

#' Plot diagnostics for a retlev object
#'
#' \code{plot} method for an objects of class \code{c("retlev", "lax")}.
#'
#' @param x an object of class \code{c("retlev", "lax")}, a result of
#'   a call to \code{\link{return_level}}, using \code{prof = TRUE}.
#' @param y Not used.
#' @param level A numeric scalar in (0, 1).  The confidence level required for
#'   the confidence interval for the \code{m}-year return level.
#'   If \code{level} is not supplied then \code{x$level} is used.
#'   \code{level} must be no larger than \code{x$level}.
#' @param legend A logical scalar.  Should we add a legend (in the top right
#'   of the plot) that gives the approximate values of the MLE and
#'   100\code{level}\% confidence limits?
#' @param digits An integer. Passed to \code{\link[base:Round]{signif}} to
#'   round the values in the legend.
#' @param plot A logical scalar.  If \code{TRUE} then the plot is produced.
#'   Otherwise, it is not, but the MLE and confidence limits are returned.
#' @param ... Further arguments to be passed to
#'   \code{\link[graphics:plot.default]{plot}}.
#' @details Plots the profile log-likelihood for a return level, provided that
#'   \code{x} returned by a call to \code{\link{return_level}} using
#'   \code{prof = TRUE}.  Horizontal lines indicate the values of the
#'   maximised log-likelihood and the critical level used to calculate
#'   the confidence limits.
#'   If \code{level} is smaller than \code{x$level} then approximate
#'   100\code{level}\% confidence limits are recalculated based on the
#'   information contained in \code{x$for_plot}.
#' @return A numeric vector of length 3 containing the lower
#'   100\code{level}\% confidence limit, the MLE and the upper
#'   100\code{level}\% confidence limit.
#' @seealso \code{\link{return_level}} to perform inferences about return
#'   levels.
#' @section Examples:
#' See the examples in \code{\link{return_level}}.
#' @export
plot.retlev <- function(x, y = NULL, level = NULL, legend = TRUE, digits = 3,
                        plot= TRUE, ...) {
  if (!inherits(x, "retlev")) {
    stop("use only with \"retlev\" objects")
  }
  if (!inherits(x, "lax")) {
    stop("use only with \"lax\" objects")
  }
  if (is.null(x$rl_prof)) {
    stop("No prof loglik info: call return_level() using prof = TRUE")
  }
  # If level is NULL then we use the intervals stored in x
  # Otherwise, we recalculate the confidence intervals
  if (is.null(level)) {
    level <- x$level
    crit_value <- x$crit
    low_lim <- x$rl_prof["lower"]
    up_lim <- x$rl_prof["upper"]
  } else if (level > x$level) {
    stop("level must be no larger than x$level")
  } else {
    crit_value <- x$max_loglik - 0.5 * stats::qchisq(level, 1)
    # Find where the curve crosses conf_line
    prof_loglik <- x$for_plot[, "prof_loglik"]
    ret_levs <- x$for_plot[, "ret_levs"]
    temp <- diff(prof_loglik - crit_value > 0)
    # Find the upper limit of the confidence interval
    loc <- which(temp == -1)
    x1 <- ret_levs[loc]
    x2 <- ret_levs[loc + 1]
    y1 <- prof_loglik[loc]
    y2 <- prof_loglik[loc + 1]
    up_lim <- x1 + (crit_value - y1) * (x2 - x1) / (y2 - y1)
    # Find the lower limit of the confidence interval
    loc <- which(temp == 1)
    x1 <- ret_levs[loc]
    x2 <- ret_levs[loc+1]
    y1 <- prof_loglik[loc]
    y2 <- prof_loglik[loc+1]
    low_lim <- x1 + (crit_value - y1) * (x2 - x1) / (y2 - y1)
  }
  my_xlab <- paste0(x$m, "-year return level")
  my_ylab <- "profile log-likelihood"
  my_plot <- function(x, y, ..., xlab = my_xlab, ylab = my_ylab, type = "l") {
    graphics::plot(x, y, ..., xlab = xlab, ylab = ylab, type = type)
  }
  hline <- function(x, ..., col = "blue", lty = 2) {
    graphics::abline(h = x, ..., col = col, lty = lty)
  }
  if (plot) {
    my_plot(x$for_plot[, "ret_levs"], x$for_plot[, "prof_loglik"], ...)
    hline(x$max_loglik, ...)
    hline(crit_value, ...)
    # Add a legend, if requested
    if (legend && length(level) == 1) {
      mle_leg <- paste0("     MLE ", signif(x$rl_prof["mle"], digits))
      conf_leg <- paste0(100 * x$level, "% CI (", signif(low_lim, digits), ",",
                         signif(up_lim, digits), ")")
      graphics::legend("topright", legend = c(mle_leg, conf_leg))
    }
  }
  res <- c(low_lim, x$rl_prof["mle"], up_lim)
  names(res) <- c("lower", "mle", "upper")
  return(res)
}

# ------------------------------ print.retlev ------------------------------- #

#' Print method for retlev object
#'
#' \code{print} method for an objects of class \code{c("retlev", "lax")}.
#'
#' @param x an object of class \code{c("retlev", "lax")}, a result of
#'   a call to \code{\link{return_level}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @param ... Additional arguments.  None are used in this function.
#' @details Prints the call to \code{\link{return_level}} and the estimates
#'   and 100\code{x$level}\% confidence limits for the \code{x$m}-year
#'   return level.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{return_level}}.
#' @section Examples:
#' See the examples in \code{\link{return_level}}.
#' @export
print.retlev <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "retlev")) {
    stop("use only with \"retlev\" objects")
  }
  if (!inherits(x, "lax")) {
    stop("use only with \"lax\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("MLE and ", 100 * x$level, "% confidence limits for the ", x$m,
      "-year return level\n\n", sep = "")
  cat("Normal interval:\n")
  print.default(format(x$rl_sym, digits = digits), print.gap = 2L,
                quote = FALSE)
  if (!is.null(x$rl_prof[1])) {
    cat("\n Profile likelihood-based interval:\n")
    print.default(format(x$rl_prof, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  return(invisible(x))
}

# ----------------------------- summary.retlev ------------------------------ #

#' Summary method for a \code{"retlev"} object
#'
#' \code{summary} method for an objects of class \code{c("retlev", "lax")}.
#'
#' @param object an object of class \code{c("retlev", "lax")}, a result of
#'   a call to \code{\link{return_level}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.  If \code{digits} is not specified
#'   (i.e. \code{\link{missing}}) then \code{signif()} will not be called
#'   (i.e. no rounding will be performed).
#' @param ... Additional arguments.  None are used in this function.
#' @return Returns a list containing the list element \code{object$call}
#'   and a numeric matrix \code{matrix} containing the MLE and estimated
#'   SE of the return level.
#' @seealso \code{\link{return_level}}.
#' @section Examples:
#' See the examples in \code{\link{return_level}}.
#' @export
summary.retlev <- function(object, digits, ...) {
  if (!inherits(object, "retlev")) {
    stop("use only with \"retlev\" objects")
  }
  if (!inherits(object, "lax")) {
    stop("use only with \"lax\" objects")
  }
  res <- object["call"]
  if (missing(digits)) {
    res$matrix <- cbind(`Estimate` = object$rl_sym["mle"],
                        `Std. Error` = object$rl_se)
  } else {
    res$matrix <- cbind(`Estimate` = signif(object$rl_sym["mle"],
                                            digits = digits),
                        `Std. Error` = signif(object$rl_se, digits = digits))
  }
  rownames(res$matrix) <- paste0("m = ", object$m)
  class(res) <- "summary.retlev"
  return(res)
}

# ---------------------------- print.summary.spm ---------------------------- #

#' Print method for objects of class \code{"summary.retlev"}
#'
#' \code{print} method for an object \code{x} of class \code{"summary.retlev"}.
#'
#' @param x An object of class "summary.retlev", a result of a call to
#'   \code{\link{summary.retlev}}.
#' @param ... Additional arguments passed on to \code{\link{print.default}}.
#' @details Prints the call and the numeric matrix \code{x$matrix} returned from
#'   \code{\link{summary.retlev}}.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{return_level}} to perform inferences about return
#'   levels.
#' @section Examples:
#' See the examples in \code{\link{return_level}}.
#' @export
print.summary.retlev <- function(x, ...) {
  if (!inherits(x, "summary.retlev")) {
    stop("use only with \"summary.retlev\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(x$matrix, ...)
  if (!is.null(x$warning)) {
    cat("\n")
    cat(x$warning)
  }
  invisible(x)
}
