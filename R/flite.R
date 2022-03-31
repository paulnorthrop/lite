#' Frequentist inference for time series extremes
#'
#' Explain
#'
#' @param data A numeric vector or numeric matrix of raw data.  If \code{data}
#'   is a matrix then the log-likelihood is constructed as the sum of
#'   (independent) contributions from different columns. A common situation is
#'   where each column relates to a different year.
#'
#'   If \code{data} contains missing values then \code{\link{split_by_NAs}} is
#'   used to divide the data further into sequences of non-missing values,
#'   stored in different columns in a matrix.  Again, the log-likelihood
#'   is constructed as a sum of contributions from different columns.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param cluster Explain ... vector if data is a vector, matrix if data is a matrix
#' @param ... Arguments to be passed to \code{\link[exdex]{kgaps}}, namely
#'   \code{k} and/or \code{inc_cens}.
#' @details Add details
#' @examples
#' cdata <- exdex::cheeseboro
#' u <- quantile(cdata, probs = 0.9, na.rm = TRUE)
#' cloglik <- flite(cdata, u)
#' summary(cloglik)
#' @export
flite <- function(data, u, cluster, ...) {
  #
  # 1. Check and manipulate data and cluster
  #
  # If cluster is missing and data is a matrix then set cluster so that data
  # in different columns are in different clusters.
  # If cluster is missing and data is a vector then set cluster so that each
  # observation forms its own cluster.
  # If cluster is supplied then check that it is consistent with data
  if (missing(cluster)) {
    if (is.matrix(data)) {
      cluster <- rep(1:ncol(data), each = nrow(data))
    } else {
      cluster <- 1:length(data)
    }
  } else {
    if (is.matrix(data)) {
      if (dim(cluster) != dim(data)) {
        stop("cluster is a matrix: dim(cluster) must equal dim(data)")
      }
      # Make cluster a vector
      cluster <- as.vector(cluster)
    } else {
      if (length(cluster) != length(data)) {
        stop("cluster is a vector: length(cluster) must equal length(data)")
      }
    }
  }
  # For GP and Bernoulli inferences we make data a vector and omit missings
  data_vector <- as.vector(data)
  not_missing <- !is.na(data_vector)
  data_vector <- data_vector[not_missing]
  cluster <- cluster[not_missing]
  # Subset cluster to contain only values for threshold exceedances
  gp_cluster <- cluster[data_vector > u]
  #
  # 2. Fit GP(sigma_u, xi) to excesses of u using an independence
  #    log-likelihood and adjust the inferences using the chandwich package
  #
  # The returned object has class "GP"
  gp_fit <- fitGP(data_vector, u)
  aloglik_gp <- adjust_object(gp_fit, cluster = gp_cluster, ...)
  #
  # 3. Fit Bernoulli(p_u) to indicators of threshold exceedance and adjust
  #    the inferences using chandwich::adjust_loglik()
  #
  # The returned object has class "Bernoulli"
  bernoulli_fit <- fitBernoulli(data_vector > u)
  aloglik_bernoulli <- adjust_object(bernoulli_fit, cluster = cluster, ...)
  #
  # 4. Make inferences about the extremal index theta using the K-gaps model
  #    via the exdex::kgaps()
  #
  # Check k
#  if (!is.wholenumber(k) || k < 0 || length(k) != 1) {
#    stop("k must be a positive integer")
#  }
  # exdex::kgaps() accept the original data matrix and omit missings itself
  theta_fit <- exdex::kgaps(data = data, u = u, k = 1, inc_cens = TRUE)
  loglik_theta <- function(tval) {
    theta_list <- c(list(theta = tval), theta_fit$ss)
    return(do.call(kgaps_loglik, theta_list))
  }
  #
  # 5. Put GP, binomial and K-gaps inferences together to return an adjusted
  #    log-likelihood for the parameters (sigma_u, xi, p_u, theta)
  #
  bernoulli_gp_theta_loglik <- function(pars, type) {
    bernoulli_loglik <- aloglik_bernoulli(pars[1], type = type)
    gp_loglik <- aloglik_gp(pars[2:3], type = type)
    theta_loglik <- loglik_theta(pars[4])
    return(bin_loglik + gp_negloglik + theta_negloglik)
  }
  class(bernoulli_gp_theta_loglik) <- c("lite", "chandwich",
                                        "bernoulli_gp_theta")
  attr(bernoulli_gp_theta_loglik, "bernoulli") <- aloglik_bernoulli
  attr(bernoulli_gp_theta_loglik, "gp") <- aloglik_gp
  attr(bernoulli_gp_theta_loglik, "kgaps") <- theta_fit
  attr(bernoulli_gp_theta_loglik, "call") <- match.call(expand.dots = TRUE)
  return(bernoulli_gp_theta_loglik)
}
