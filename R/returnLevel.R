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
#' @param ny A numeric scalar.  The (mean) number of observations per year.
#'   \strong{Setting this appropriately is important}. See \strong{Details}.
#' @param prof A logical scalar.  Should we calculate intervals based on
#'   profile log-likelihood?
#' @param inc A numeric scalar. Only relevant if \code{prof = TRUE}. The
#'   increment in return level by which we move upwards and downwards from the
#'   MLE for the return level in the search for the lower and upper confidence
#'   limits.  If this is not supplied then \code{inc} is set to one hundredth
#'   of the length of the symmetric confidence interval for return level.
#' @param type A character scalar.  The argument \code{type} to the function
#'   returned by the function \code{\link[chandwich]{adjust_loglik}}, that is, the type of
#'   adjustment made to the independence log-likelihood function in creating
#'   an adjusted log-likelihood function.  See \strong{Details} and
#'   \strong{Value} in \code{\link[chandwich]{adjust_loglik}}.
#' @details For information about return levels see the "Introducing lite"
#'   vignette.
#'
#'   \code{ny} provides information about the (intended) frequency of
#'   sampling in time, that is, the number of observations that would be
#'   observed in a year if there are no missing values.  If the number of
#'   observations may vary between years then \code{ny} should be set equal to
#'   the mean number of observations per year.
#'
#'   \strong{Supplying \code{ny}.}
#'   The value of \code{ny} may have been set in the call to
#'   \code{\link{flite}}.  If \code{ny} is supplied by the user in the call to
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
#' @return A object (a list) of class \code{"returnLevel", "lite"} with the
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
#'   \item{ny }{The value of \code{ny} used to infer the return level.}
#'   \item{call }{The call to \code{returnLevel}.}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \doi{10.1007/978-1-4471-3675-0_3}
#' @seealso \code{\link{returnLevelMethods}}, including plotting the (adjusted)
#'   profile log-likelihood for a return level.
#' @examples
#' ### Cheeseboro wind gusts
#'
#' # Make inferences
#' cdata <- exdex::cheeseboro
#' # Each column of the matrix cdata corresponds to data from a different year
#' # flite() sets cluster automatically to correspond to column (year)
#' cfit <- flite(cdata, u = 45, k = 3)
#'
#' # These data are hourly for one month (January) year so ny = 31 * 24
#' # Large inc set here for speed, sacrificing accuracy
#' # Default 95% confidence intervals
#' rl <- returnLevel(cfit, inc = 2.5, ny = 31 * 24)
#' summary(rl)
#' rl
#' oldrl <- plot(rl)
#' oldrl
#'
#' # Quickly recalculate/replot the intervals based on profile log-likelihood
#' # provided that level is smaller than that used to produce rl
#' newrl <- plot(rl, level = 0.9)
#' newrl
#' @export
returnLevel <- function(x, m = 100, level = 0.95, ny, prof = TRUE,
                         inc = NULL,
                         type = c("vertical", "cholesky", "spectral", "none")) {
  if (!inherits(x, "flite")) {
    stop("use only with \"flite\" objects")
  }
  Call <- match.call(expand.dots = TRUE)
  type <- match.arg(type)
  # Check whether ny is supplied in the call to returnLevel
  ny_given <- ifelse(missing(ny), FALSE, TRUE)
  # Make inferences about return levels
  temp <- return_level_bingp(x, m, level, ny, prof, inc, type, ny_given)
  temp$m <- m
  temp$level <- level
  temp$ny <- ny
  temp$call <- Call
  class(temp) <- c("returnLevel", "lite")
  return(temp)
}
