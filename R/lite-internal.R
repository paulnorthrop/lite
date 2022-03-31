#' Internal lite functions
#'
#' Internal lite functions
#' @details
#' These functions are not intended to be called by the user.
#' @name lite-internal
#' @keywords internal
NULL

# ================================= gp_mle ================================== #

#' @keywords internal
#' @rdname lite-internal
gp_mle <- function(gp_data, u) {
  # Maximum likelihood estimation for the generalized Pareto distribution
  #
  # Performs maximum likelihood estimation for the generalized Pareto
  # distribution.  Uses the function \code{gpdmle} associated with
  # Grimshaw (1993), which returns MLEs of sigma and k = - \xi.
  #
  # Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
  #   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
  #   and Computing (1991) 1, 129-133. https://doi.org/10.1007/BF01889987.
  #
  # Args:
  #   gp_data : A numeric vector containing positive values, assumed to be a
  #             random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A list with components
  #     mle  : A numeric vector.  MLEs of GP parameters sigma and xi.
  #     nllh : A numeric scalar.  The negated log-likelihood at the MLE.
  #
  # Create excesses of threshold u
  excesses <- gp_data[gp_data > u] - u
  # Call Grimshaw (1993) function, note: k is -xi, a is sigma
  grimshaw_fit <- revdbayes::grimshaw_gp_mle(excesses)
  res <- list()
  # mle for (sigma, xi)
  res$mle <- c(sigma_u = grimshaw_fit$a, xi = -grimshaw_fit$k)
  # number of threshold excesses
  res$nexc <- length(excesses)
  sc <- rep_len(res$mle[1], res$nexc)
  xi <- res$mle[2]
  res$cov <- solve(gp_obs_info(res$mle, excesses))
  res$nllh <- sum(log(sc)) + sum(log(1 + xi * excesses / sc) * (1 / xi + 1))
  # Return the exceedances (values that lie above the threshold)
  res$data <- gp_data[gp_data > u]
  res$threshold <- u
  class(res) <- "GP"
  return(res)
}

# =========================== gp_obs_info ===========================

#' @keywords internal
#' @rdname lite-internal
gp_obs_info <- function(gp_pars, y) {
  # Observed information for the generalized Pareto distribution
  #
  # Calculates the observed information matrix for a random sample \code{y}
  # from the generalized Pareto distribution, i.e. the negated Hessian matrix
  # of the generalized Pareto log-likelihood, evaluated at \code{gp_pars}.
  #
  # Args:
  #   gp_pars : A numeric vector. Parameters sigma and xi of the
  #   generalized Pareto distribution.
  #   y       : A numeric vector. A sample of positive data values.
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  s <- gp_pars[1]
  x <- gp_pars[2]
  i <- matrix(NA, 2, 2)
  i[1,1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1,2] <- i[2,1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  i[2,2] <- sum(2 * log(1 + x * y / s) / x ^ 3 - 2 * y / (s + x * y) / x ^ 2 -
                  (1 + 1 / x) * y ^ 2 / (s + x * y) ^ 2)
  return(i)
}

# ====================== Log-likelihood adjustment function ================= #

#' @keywords internal
#' @rdname lite-internal
adjust_object <- function(x, cluster = NULL, ...) {
  #
  # Set logLikVector function
  #
  loglik_fn <- function(pars, fitted_object, ...) {
    return(logLikVector(fitted_object, pars = pars))
  }
  #
  # Set H ----------
  #
  H <- -solve(vcov(x))
  #
  # Set mle and nobs ----------
  #
  mle <- coef(x)
  n_obs <- nobs(x)
  #
  # Set V, using meat() or meatCL() from the sandwich package ----------
  #
  if (is.null(cluster)) {
    V <- sandwich::meat(x, fitted_object = x, loglik_fn = loglik_fn,
                        ...) * n_obs
  } else {
    V <- sandwich::meatCL(x, cluster = cluster, fitted_object = x,
                          loglik_fn = loglik_fn, ...) * n_obs
  }
  # We don't pass cluster because it would only be used in the estimation of
  # V: we have already estimated V using sandwich::meat() or sandwich::meatCL()
  res <- chandwich::adjust_loglik(loglik = loglik_fn,
                                  fitted_object = x,
                                  p = length(mle),
                                  par_names = names(mle),
                                  name = paste(class(x), collapse = "_"),
                                  mle = mle, H = H, V = V)
  # If cluster was supplied then overwrite the default (1,2, ...) returned by
  # chandwich::adjust_loglik()
  if (!is.null(cluster)) {
    attr(res, "cluster") <- cluster
  }
  # Add the original fitted model object as an attribute
  attr(res, "original_fit") <- x
  class(res) <- c("lite", "chandwich", class(x))
  return(res)
}

# =============================== kgaps_loglik ============================== #
# Included because this is not exported from exdex
# The argument n_kgaps is not used here but it is included because it is
# included in the list returned by exdex::kgaps_stat()

#' @keywords internal
#' @rdname lite-internal
kgaps_loglik <- function(theta, N0, N1, sum_qs, n_kgaps){
  if (theta < 0 || theta > 1) {
    return(-Inf)
  }
  loglik <- 0
  if (N1 > 0) {
    loglik <- loglik + 2 * N1 * log(theta) - sum_qs * theta
  }
  if (N0 > 0) {
    loglik <- loglik + N0 * log(1 - theta)
  }
  return(loglik)
}
