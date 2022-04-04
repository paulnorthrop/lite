# Create an estfun methods for the "Bernoulli" and GP" classes.
# Leave the default method here for the moment, even though it is not used.
# It could be useful for testing.

# Default

#' @export
estfun.default <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}

# Bernoulli

#' @export
estfun.Bernoulli <- function(x, ...) {
  U <- x$data / x$mle - (1 - x$data) / (1 - x$mle)
  dim(U)<- c(length(U), 1)
  colnames(U) <- names(coef(x))
  return(U)
}

# GP

#' @export
estfun.GP <- function(x, eps = 1e-5, m = 3, ...) {
  if (eps <= 0) {
    stop("'eps' must be positive")
  }
  if (m < 0) {
    stop("'m' must be non-negative")
  }
  pars <- coef(x)
  sigma <- pars[1]
  xi <- pars[2]
  z <- xi / sigma
  # Threshold excesses
  y <- x$exceedances - x$threshold
  zy <- z * y
  t0 <- 1 + zy
  U <- matrix(NA, nrow = length(y), ncol = 2)
  U[, 1] <- -1 / sigma + (xi + 1) * y / (t0 * sigma ^ 2)
  if (abs(xi) < eps) {
    i <- 0:m
    sum_fn <- function(zy) {
      return(sum((-1) ^ i * (i + 1) * zy ^ i / (i + 2)))
    }
    tsum <- vapply(zy, sum_fn, 0.0)
    U[, 2] <- y ^ 2 * tsum / sigma ^ 2 - y / (sigma * t0)
  } else {
    t1 <- log1p(zy) / z ^ 2
    t2 <- y / (z * t0)
    U[, 2] <- (t1 - t2) / sigma ^ 2 - y / (sigma * t0)
  }
  colnames(U) <- names(coef(x))
  return(U)
}
