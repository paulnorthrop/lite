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
estfun.GP <- function(x, eps = 1e-5, ...) {
  print("Using my code")
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
    s1 <- 1 / 2
    s2 <- 2 * z * y / 3
    s3 <- 3 * z ^ 2 * y ^ 2 / 4
    U[, 2] <- y ^ 2 * (s1 - s2 + s3) / sigma ^ 2 - y / (sigma * t0)
  } else {
    t1 <- log(t0) / z ^ 2
    t2 <- y / (z * t0)
    U[, 2] <- (t1 - t2) / sigma ^ 2 - y / (sigma * t0)
  }
  colnames(U) <- names(coef(x))
  return(U)
}
