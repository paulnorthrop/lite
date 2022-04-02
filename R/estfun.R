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
estfun.GP <- function(x, loglik_fn, ...) {
  U <- numDeriv::jacobian(loglik_fn, x = coef(x), ...)
  colnames(U) <- names(coef(x))
  return(U)
}
