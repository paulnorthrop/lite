#' Methods for objects of class \code{"returnLevel"}
#'
#' Methods for obects of class \code{"returnLevel"} returned from
#' \code{\link{returnLevel}}.
#' @param ... For \code{plot.returnLevel}: arguments passed to
#'   \code{\link[graphics:plot.default]{plot}}, such as graphical parameters.
#'
#'   For \code{print.summary.returnLevel}: additional arguments passed to
#'   \code{\link{print.default}}.
#' @return \code{plot.returnLevel}: a numeric vector of length 3 containing the
#'   lower 100\code{level}\% confidence limit, the MLE and the upper
#'   100\code{level}\% confidence limit is returned invisibly.
#'
#'   \code{print.returnLevel}: the argument \code{object}, invisibly.
#'
#'   \code{summary.returnLevel}: a list containing the list element
#'   \code{object$call} and a numeric matrix \code{matrix} containing the MLE
#'   and estimated SE of the return level.
#'
#'   \code{print.summary.returnLevel}: the argument \code{object}, invisibly.
#' @name returnLevelMethods
NULL
## NULL

# ------------------------------- plot.returnLevel ------------------------------- #

#' Plot diagnostics for a returnLevel object
#'
#' @param object an object of class \code{c("returnLevel", "lite")}, a result
#'   of a call to \code{\link{returnLevel}}, using \code{prof = TRUE}.
#' @param level A numeric scalar in (0, 1).  The confidence level required for
#'   the confidence interval for the \code{m}-year return level.
#'   If \code{level} is not supplied then \code{object$level} is used.
#'   \code{level} must be no larger than \code{object$level}.
#' @param legend A logical scalar.  Should we add a legend (in the top right
#'   of the plot) that gives the approximate values of the MLE and
#'   100\code{level}\% confidence limits?
#' @param digits An integer. Passed to \code{\link[base:Round]{signif}} to
#'   round the values in the legend.
#' @param plot A logical scalar.  If \code{TRUE} then the plot is produced.
#'   Otherwise, it is not, but the MLE and confidence limits are returned.
#' @details \code{plot.returnLevel} plots the profile log-likelihood for a
#'   return level, provided that \code{object} returned by a call to
#'   \code{\link{returnLevel}} using \code{prof = TRUE}.  Horizontal lines
#'   indicate the values of the maximised log-likelihood and the critical level
#'   used to calculate the confidence limits.
#'   If \code{level} is smaller than \code{object$level} then approximate
#'   100\code{level}\% confidence limits are recalculated based on the
#'   information contained in \code{object$for_plot}.
#' @seealso \code{\link{returnLevel}} to perform inferences about return
#'   levels.
#' @section Examples:
#' See \code{\link{returnLevel}}.
#' @rdname returnLevelMethods
#' @export
plot.returnLevel <- function(object, level = NULL, legend = TRUE, digits = 3,
                             plot = TRUE, ...) {
  if (!inherits(object, "returnLevel")) {
    stop("use only with \"returnLevel\" objects")
  }
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  if (is.null(object$rl_prof)) {
    stop("No prof loglik info: call returnLevel() using prof = TRUE")
  }
  # If level is NULL then we use the intervals stored in object
  # Otherwise, we recalculate the confidence intervals
  if (is.null(level)) {
    level <- object$level
    crit_value <- object$crit
    low_lim <- object$rl_prof["lower"]
    up_lim <- object$rl_prof["upper"]
  } else if (level > object$level) {
    stop("level must be no larger than object$level")
  } else {
    crit_value <- object$max_loglik - 0.5 * stats::qchisq(level, 1)
    # Find where the curve crosses conf_line
    prof_loglik <- object$for_plot[, "prof_loglik"]
    ret_levs <- object$for_plot[, "ret_levs"]
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
  my_xlab <- paste0(object$m, "-year return level")
  my_ylab <- "profile log-likelihood"
  my_plot <- function(x, y, ..., xlab = my_xlab, ylab = my_ylab,
                      type = "l") {
    graphics::plot(x, y, ..., xlab = xlab, ylab = ylab, type = type)
  }
  hline <- function(x, ..., col = "blue", lty = 2) {
    graphics::abline(h = x, ..., col = col, lty = lty)
  }
  if (plot) {
    my_plot(object$for_plot[, "ret_levs"], object$for_plot[, "prof_loglik"],
            ...)
    hline(object$max_loglik, ...)
    hline(crit_value, ...)
    # Add a legend, if requested
    if (legend && length(level) == 1) {
      mle_leg <- paste0("     MLE ", signif(object$rl_prof["mle"], digits))
      conf_leg <- paste0(100 * object$level, "% CI (", signif(low_lim, digits),
                         ",", signif(up_lim, digits), ")")
      graphics::legend("topright", legend = c(mle_leg, conf_leg))
    }
  }
  res <- c(low_lim, object$rl_prof["mle"], up_lim)
  names(res) <- c("lower", "mle", "upper")
  return(invisible(res))
}

# ------------------------------ print.returnLevel ------------------------------- #

#' Print method for returnLevel object
#'
#' @param object an object of class \code{c("returnLevel", "lite")}, a result of
#'   a call to \code{\link{returnLevel}}.
#' @param digits The argument \code{digits} to \code{\link{print.default}}.
#' @details \code{print.returnLevel} prints the call to
#'   \code{\link{returnLevel}} and the estimates and 100\code{object$level}\%
#'   confidence limits for the \code{object$m}-year return level.
#' @rdname returnLevelMethods
#' @export
print.returnLevel <- function(object, digits =
                                max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(object, "returnLevel")) {
    stop("use only with \"returnLevel\" objects")
  }
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
  }
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("MLE and ", 100 * object$level, "% confidence limits for the ", object$m,
      "-year return level\n\n", sep = "")
  cat("Normal interval:\n")
  print.default(format(object$rl_sym, digits = digits), print.gap = 2L,
                quote = FALSE)
  if (!is.null(object$rl_prof[1])) {
    cat("\n Profile likelihood-based interval:\n")
    print.default(format(object$rl_prof, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  return(invisible(object))
}

# ----------------------------- summary.returnLevel ------------------------------ #

#' Summary method for a \code{"returnLevel"} object
#'
#' @param object an object of class \code{c("returnLevel", "lite")}, a result of
#'   a call to \code{\link{returnLevel}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base:Round]{signif}}.  If \code{digits} is not specified
#'   (i.e. \code{\link{missing}}) then \code{signif()} will not be called
#'   (i.e. no rounding will be performed).
#' @rdname returnLevelMethods
#' @export
summary.returnLevel <- function(object, digits, ...) {
  if (!inherits(object, "returnLevel")) {
    stop("use only with \"returnLevel\" objects")
  }
  if (!inherits(object, "lite")) {
    stop("use only with \"lite\" objects")
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
  class(res) <- "summary.returnLevel"
  return(res)
}

# ---------------------------- print.summary.spm ---------------------------- #

#' Print method for objects of class \code{"summary.returnLevel"}
#'
#' @param object An object of class "summary.returnLevel", a result of a call
#'   to \code{\link{summary.returnLevel}}.
#' @rdname returnLevelMethods
#' @export
print.summary.returnLevel <- function(object, ...) {
  if (!inherits(object, "summary.returnLevel")) {
    stop("use only with \"summary.returnLevel\" objects")
  }
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  print(object$matrix, ...)
  if (!is.null(object$warning)) {
    cat("\n")
    cat(object$warning)
  }
  invisible(object)
}
