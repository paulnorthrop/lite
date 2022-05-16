# =========================== predict.blite ================================= #

#' Predictive inference for the largest value observed in \eqn{N} years.
#'
#' \code{predict} method for class "blite".  Performs predictive inference
#' about the largest value to be observed over a future time period of
#' \eqn{N} years.  Predictive inferences accounts for uncertainty in model
#' parameters and for uncertainty owing to the variability of future
#' observations.
#'
#' @inheritParams revdbayes::predict.evpost
#' @param object An object of class \code{"blite"} returned from
#'   \code{\link{blite}}.
#' @param x A numeric vector or a matrix with \code{n_years} columns.
#'   The meaning of \code{x} depends on \code{type}.
#'   \itemize{
#'     \item{\code{type = "p"} or \code{type = "d"}:} \code{x} contains
#'       quantiles at which to evaluate the distribution or density function.
#'       No element of \code{x} can be less than the threshold
#'       \code{attr(object, "inputs")$u}.
#'
#'       If \code{x} is not supplied then \code{n_year}-specific defaults are
#'       set: vectors of length \code{x_num} from the 0.1\% quantile to the
#'       99\% quantile, subject all values being greater than the threshold.
#'     \item{\code{type = "q"}:} \code{x} contains probabilities in (0,1)
#'       at which to evaluate the quantile function.  Any values outside
#'       (0, 1) will be removed without warning. No element of \code{p} can
#'       correspond to a predictive quantile that is below the threshold,
#'       \code{attr(object, "inputs")$u}.  That is, no element of \code{p} can
#'       be less than the value of \code{predict.evpost(object,}
#'       \code{type = "q", x = attr(object, "inputs")$u)}.
#'
#'       If \code{x} is not supplied then a default value of
#'       \code{c(0.025, 0.25, 0.5, 0.75, 0.975)} is used.
#'     \item{\code{type = "i"} or \code{type = "r"}:} \code{x} is not relevant.
#'   }
#' @param ny A numeric scalar.  The (mean) number of observations per year.
#'   \strong{Setting this appropriately is important}. See \strong{Details}.
#' @details The function \code{\link[revdbayes]{predict.evpost}} in the
#'   \code{\link[revdbayes]{revdbayes}} package is used to perform the
#'   predictive inferences.  The effect of adjusting for the values of the
#'   extremal index \eqn{\theta} in the posterior sample in
#'   \code{object$sim_vals[, "theta]} is to change the effective time horizon
#'   from \eqn{N} to \eqn{\theta N}.
#'
#'   \code{ny} provides information about the (intended) frequency of
#'   sampling in time, that is, the number of observations that would be
#'   observed in a year if there are no missing values.  If the number of
#'   observations may vary between years then \code{ny} should be set equal to
#'   the mean number of observations per year.
#'
#'   \strong{Supplying \code{ny}.}
#'   The value of \code{ny} may have been set in the call to
#'   \code{\link{blite}}.  If \code{ny} is supplied by the user in the call to
#'   \code{predict.blite} then this will be used in preference to the value
#'   stored in the fitted model object.  If these two values differ then no
#'   warning will be given.
#' @return An object of class "evpred", a list containing a subset of the
#'   following components:
#'     \item{type}{The argument \code{type} supplied to \code{predict.blite}.
#'     Which of the following components are present depends \code{type}.}
#'     \item{x}{A matrix containing the argument \code{x} supplied to
#'       \code{predict.blite}, or set within \code{predict.blite} if \code{x}
#'       was not supplied, replicated to have \code{n_years} columns
#'       if necessary.
#'       Only present if \code{type} is \code{"p", "d"} or \code{"q"}.}
#'     \item{y}{The content of \code{y} depends on \code{type}:
#'     \itemize{
#'       \item{\code{type = "p", "d", "q"}:}  A matrix with the same
#'       dimensions as \code{x}.  Contains distribution function values
#'       (\code{type = "p"}), predictive density (\code{type = "d"})
#'       or quantiles (\code{type = "q"}).
#'       \item{\code{type = "r"}:} A numeric matrix with \code{length(n_years)}
#'       columns and number of rows equal to the size of the posterior sample.
#'       \item{\code{type = "i"}:} \code{y} is not present.
#'       }}
#'       \item{long}{A \code{length(n_years)*length(level)} by 4 numeric
#'         matrix containing the equi-tailed limits with columns:
#'         lower limit, upper limit, n_years, level.
#'         Only present if \code{type = "i"}.  If an interval extends below
#'         the threshold then \code{NA} is returned.}
#'       \item{short}{A matrix with the same structure as \code{long}
#'         containing the HPD limits.  Only present if \code{type = "i"}.
#'         Columns 1 and 2 contain \code{NA}s if \code{hpd = FALSE}
#'         or if the corresponding equi-tailed interval extends below
#'         the threshold.}
#'   The arguments \code{n_years, level, hpd, lower_tail, log} supplied
#'   to \code{predict.blite} are also included, as is the value of \code{ny}
#'   and \code{model = "bingp"}.
#' @examples
#' ### Cheeseboro wind gusts
#'
#' cdata <- exdex::cheeseboro
#' # Each column of the matrix cdata corresponds to data from a different year
#' # blite() sets cluster automatically to correspond to column (year)
#' cpost <- blite(cdata, u = 45, k = 3, ny = 31 * 24)
#'
#' # Need revdbayes v1.5.9 for the following examples to be correct
#' # (At the moment they are based on theta = 1)
#'
#' # Interval estimation
#' predict(cpost)$long
#' predict(cpost, hpd = TRUE)$short
#'
#' # Density function
#' plot(predict(cpost, type = "d", n_years = c(100, 1000)))
#'
#' # Distribution function
#' plot(predict(cpost, type = "p", n_years = c(100, 1000)))
#'
#' # Quantiles
#' predict(cpost, type = "q", n_years = c(100, 1000))$y
#'
#' # Random generation
#' plot(predict(cpost, type = "r"))
#' @export
predict.blite <- function(object, type = c("i", "p", "d", "q", "r"), x = NULL,
                           x_num = 100, n_years = 100, ny = NULL, level = 95,
                           hpd = FALSE, lower_tail = TRUE, log = FALSE,
                           big_q = 1000, ...) {
  # Create an object that has the same structure, at least in terms of things
  # that matter for predictive inference, as an object of class "evpost"
  # returned from revdbayes::rpost(..., model = "bingp")
  temp <- list()
  # Threshold
  temp$thresh <- attr(object, "inputs")$u
  # Name of model (binomial-GP)
  temp$model <- "bingp"
  # A function to evalulate the posterior density of the posterior for the
  # exceedance probability p_u, the arguments to this function and the values
  # simulated from this posterior density
  temp$bin_logf <- attr(object, "Bernoulli")$bin_logf
  temp$bin_logf_args <- attr(object, "Bernoulli")$bin_logf_args
  temp$bin_sim_vals <- attr(object, "Bernoulli")$bin_sim_vals
  # Extract the simulated values of (sigma_u, xi, theta)
  temp$sim_vals <- object[, c("sigma[u]", "xi", "theta")]
  # Check whether ny is supplied in the call to predict.blite()
  ny_given <- ifelse(missing(ny), FALSE, TRUE)
  # If ny is not supplied then check whether a value was supplied to blite()
  if (!ny_given) {
    ny <- attr(object, "inputs")$ny
  }
  if (is.na(ny)) {
    stop("'ny' has not been supplied")
  }
  # In predict.evpost() ny is called npy
  npy <- ny
  # Call predict.evpost()
  class(temp) <- "evpost"
  res <- predict(temp, type = type, x = x, x_num = x_num, n_years = n_years,
                 npy = npy, level = level, hpd = hpd, lower_tail = lower_tail,
                 log = log, big_q = big_q)
  # Rename npy to ny
  res$ny <- res$npy
  res$npy <- NULL
  return(res)
}
