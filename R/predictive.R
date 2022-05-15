# =========================== predict.blite ================================= #

#' Predictive inference for the largest value observed in N years.
#'
#' \code{predict} method for class "blite".  Performs predictive inference
#' about the largest value to be observed over a future time period of
#' N years.  Predictive inferences accounts for uncertainty in model
#' parameters and for uncertainty owing to the variability of future
#' observations.
#'
#' @inheritParams revdbayes::predict.evpost
#' @examples
#' ### Cheeseboro wind gusts
#'
#' cdata <- exdex::cheeseboro
#' # Each column of the matrix cdata corresponds to data from a different year
#' # blite() sets cluster automatically to correspond to column (year)
#' cpost <- blite(cdata, u = 45, k = 3)
#' summary(cpost)
#' # Need revdbayes v1.5.9 for the following examples
#' #predict(cpost, ny = 31 * 24)$long
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
  return(res)
}
