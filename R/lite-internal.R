#' Internal lite functions
#'
#' Internal lite functions
#' @details
#' These functions are not intended to be called by the user.
#' @name lite-internal
#' @keywords internal
NULL

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
  class(res) <- c("chandwich", class(x))
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

# ============================ check_logLik_flite =========================== #
# Included to provide a check of logLik.flite()

#' @keywords internal
#' @rdname lite-internal
check_logLik_flite <- function(object, ...) {
  if (!inherits(object, "flite")) {
    stop("use only with \"flite\" objects")
  }
  bfit <- attr(object, "bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "kgaps")
  bloglik <- attr(bfit, "max_loglik")
  gloglik <- attr(gfit, "max_loglik")
  kloglik <- kfit$max_loglik
  val <- bloglik + gloglik + kloglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- 4
  class(val) <- "logLik"
  return(val)
}
