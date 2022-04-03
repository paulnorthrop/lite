# Create an estfun method for the "GP" class default estfun method and, for safety, create individual methods
# for all the classes currently involved in lax.  At the moment they all
# use numDeriv::jacobian() to do the calculation

#' @export
estfun.default <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
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
  t0 <- 1 + z * y
  U <- matrix(NA, nrow = length(y), ncol = 2)
  U[, 1] <- -1 / sigma + (xi + 1) * y / (t0 * sigma ^ 2)
  if (abs(xi) < eps) {
    i <- 0:m
    zy <- z * y
    sum_fn <- function(zy) {
      return(sum((-1) ^ i * (i + 1) * zy ^ i / (i + 2)))
    }
    tsum <- vapply(zy, sum_fn, 0.0)
    U[, 2] <- y ^ 2 * tsum / sigma ^ 2 - y / (sigma * t0)
  } else {
    t1 <- log(t0) / z ^ 2
    t2 <- y / (z * t0)
    U[, 2] <- (t1 - t2) / sigma ^ 2 - y / (sigma * t0)
  }
  colnames(U) <- names(coef(x))
  return(U)
}
