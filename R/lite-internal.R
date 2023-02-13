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
  # Set V, using meatCL() from the sandwich package ----------
  #
  V <- sandwich::meatCL(x, cluster = cluster, fitted_object = x,
                        loglik_fn = loglik_fn, ...) * n_obs
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
  bfit <- attr(object, "Bernoulli")
  gfit <- attr(object, "gp")
  kfit <- attr(object, "theta")
  bloglik <- attr(bfit, "max_loglik")
  gloglik <- attr(gfit, "max_loglik")
  kloglik <- kfit$max_loglik
  val <- bloglik + gloglik + kloglik
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- 4
  class(val) <- "logLik"
  return(val)
}

# ==================== Binomial-GP return levels functions ================== #

#' @keywords internal
#' @rdname lite-internal
return_level_bingp <- function(x, m, level, ny, prof, inc, type, ny_given) {
  # Extract the threshold and ny from the original fitted model object
  u <- attr(x, "inputs")$u
  if (!ny_given) {
    ny <- attr(x, "inputs")$ny
  }
  if (is.na(ny)) {
    stop("'ny' has not been supplied")
  }
  # MLE and symmetric conf% CI for the return level
  rl_sym <- bingp_rl_CI(x, m, level, ny, type, u)
  # Extract SE
  rl_se <- rl_sym["se"]
  # Remove SE
  rl_sym <- rl_sym[c("lower", "mle", "upper")]
  if (!prof) {
    return(list(rl_sym = rl_sym, rl_prof = NULL, rl_se = rl_se))
  }
  # CIs based on profile log-likelihood
  # If k = 0 then use a function for which theta is set equal to 1
  if (attr(x, "theta")$k > 0) {
    temp <- bingp_rl_prof(x, m, level, ny, inc, type, rl_sym, u)
  } else {
    temp <- bingp_rl_prof0(x, m, level, ny, inc, type, rl_sym, u)
  }
  return(list(rl_sym = rl_sym, rl_prof = temp$rl_prof, rl_se = rl_se,
              max_loglik = logLik(x), crit = temp$crit,
              for_plot = temp$for_plot))
}

#' @keywords internal
#' @rdname lite-internal
bingp_rl_CI <- function (x, m, level, ny, type, u) {
  # Extract the MLEs
  mles <- coef(x)
  pu <- mles["p[u]"]
  sigmau <- mles["sigma[u]"]
  xi <- mles["xi"]
  theta <- mles["theta"]
  # Create the covariance matrix for all 4 parameters: (pu, sigmau, xi, theta)
  # Get the correct matrix by supplying adjust = FALSE for type = "none" and
  # adjust = TRUE otherwise
  adjust <- ifelse(type == "none", FALSE, TRUE)
  mat <- vcov(x, adjust = adjust)
  # pmny is approximately equal to 1 / (m * ny * theta)
  con <- 1 - 1 / m
  nt <- ny * theta
  pmny <- 1 - con ^ (1 / nt)
  p <- pmny / pu
  rp <- 1 / p
  if (p > 1) {
    stop("The return level required is below the threshold")
  }
  rl_mle <- revdbayes::qgp(p, loc = u, scale = sigmau, shape = xi,
                           lower.tail = FALSE)
  delta <- matrix(0, 4, 1)
  delta[1, ] <- sigmau * pu ^ (xi - 1) / pmny ^ xi
  delta[2, ] <- revdbayes::qgp(p, loc = 0, scale = 1, shape = xi,
                               lower.tail = FALSE)
  delta[3, ] <- sigmau * box_cox_deriv(rp, lambda = xi)
  # Add information about theta
  delta[4, ] <- -sigmau * rp ^ (xi - 1) * pu * con ^ (1 / nt) * log(con) /
    (ny * theta ^ 2 * pmny ^ 2)
  rl_var <- t(delta) %*% mat %*% delta
  rl_se <- sqrt(rl_var)
  z_val <- stats::qnorm(1 - (1 - level) / 2)
  rl_lower <- rl_mle - z_val * rl_se
  rl_upper <- rl_mle + z_val * rl_se
  res <- c(lower = rl_lower, mle = rl_mle, upper = rl_upper, se = rl_se)
  return(res)
}

#' @keywords internal
#' @rdname lite-internal
bingp_rl_prof <- function(x, m, level, ny, inc, type, rl_sym, u) {
  # Convert inc to a length relative to the length of the symmetric CI
  inc <- (rl_sym["upper"] - rl_sym["lower"]) * inc
  bingptheta_mle <- coef(x)
  # Function to calculate the negated adjusted binomial-GP log-likelihood
  bingp_negloglik <- function(pars, type) {
    return(-x(pars, type = type))
  }
  # Calculates the negated profile log-likelihood of the m-year return level
  bingp_neg_prof_loglik <- function(a, xp) {
    # a[1] is pu, a[2] is xi, a[3] is theta
    # Check that pu is in (0, 1) and theta is in (0, 1]
    if (a[1] <= 0 || a[1] >= 1 || a[3] <= 0 || a[3] > 1) {
      return(10 ^ 10)
    }
    pmny <- 1 - (1 - 1 / m) ^ (1 / (ny * a[3]))
    p <- pmny / a[1]
    # Check that p is in [0, 1]
    if (p > 1 || p < 0) {
      return(10 ^ 10)
    }
    sigmau <- (xp - u) / revdbayes::qgp(p, loc = 0, scale = 1, shape = a[2],
                                        lower.tail = FALSE)
    # Check that sigmau is positive
    if (sigmau <= 0) {
      return(10 ^ 10)
    }
    bingp_pars <- c(a[1], sigmau, a[2], a[3])
    return(bingp_negloglik(bingp_pars, type = type))
  }
  max_loglik <- logLik(x)
  rl_mle <- rl_sym["mle"]
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- rl_mle
  v2[1] <- v1[1] <- max_loglik
  #
  # Starting from the MLE, we search upwards and downwards until we pass the
  # cutoff for the 100level% confidence interval
  #
  ### Upper tail ...
  xp <- rl_mle
  my_val <- max_loglik
  ii <- 1
  sol <- bingptheta_mle[-2]
  while (my_val > conf_line){
    xp <- xp + inc
    opt <- stats::optim(sol, bingp_neg_prof_loglik, method = "BFGS", xp = xp)
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- xp
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }
  sol_up <- sol
  ### Lower tail ...
  xp <- rl_mle
  my_val <- max_loglik
  ii <- 1
  sol <- bingptheta_mle[-2]
  while (my_val > conf_line){
    xp <- xp - inc
    opt <- stats::optim(sol, bingp_neg_prof_loglik, method = "BFGS", xp = xp)
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- xp
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  sol_low <- sol
  #
  # Find the limits of the confidence interval
  #
  prof_lik <- c(rev(v1), v2)
  ret_levs <- c(rev(x1), x2)
  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc <- which(temp == -1)
  x1 <- ret_levs[loc]
  x2 <- ret_levs[loc + 1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc + 1]
  up_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  # Find the lower limit of the confidence interval
  loc <- which(temp == 1)
  x1 <- ret_levs[loc]
  x2 <- ret_levs[loc+1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc+1]
  low_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  rl_prof <- c(lower = low_lim, rl_mle, upper = up_lim)
  return(list(rl_prof = rl_prof, crit = conf_line,
              for_plot = cbind(ret_levs = ret_levs, prof_loglik = prof_lik)))
}

# Profile log-likelihood-based CIs in the K = 0 case, where theta is 1

#' @keywords internal
#' @rdname lite-internal
bingp_rl_prof0 <- function(x, m, level, ny, inc, type, rl_sym, u) {
  if (is.null(inc)) {
    inc <- (rl_sym["upper"] - rl_sym["lower"]) / 100
  }
  bingptheta_mle <- coef(x)
  # Function to calculate the negated adjusted binomial-GP log-likelihood
  bingp_negloglik <- function(pars, type) {
    return(-x(pars, type = type))
  }
  # Calculates the negated profile log-likelihood of the m-year return level
  bingp_neg_prof_loglik <- function(a, xp) {
    # a[1] is pu, a[2] is xi
    # Check that pu is in (0, 1)
    if (a[1] <= 0 || a[1] >= 1) {
      return(10 ^ 10)
    }
    pmny <- 1 - (1 - 1 / m) ^ (1 / ny)
    p <- pmny / a[1]
    # Check that p is in [0, 1]
    if (p > 1 || p < 0) {
      return(10 ^ 10)
    }
    sigmau <- (xp - u) / revdbayes::qgp(p, loc = 0, scale = 1, shape = a[2],
                                        lower.tail = FALSE)
    # Check that sigmau is positive
    if (sigmau <= 0) {
      return(10 ^ 10)
    }
    bingp_pars <- c(a[1], sigmau, a[2], 1)
    return(bingp_negloglik(bingp_pars, type = type))
  }
  max_loglik <- logLik(x)
  rl_mle <- rl_sym["mle"]
  conf_line <- max_loglik - 0.5 * stats::qchisq(level, 1)
  v1 <- v2 <- x1 <- x2 <- NULL
  x2[1] <- x1[1] <- rl_mle
  v2[1] <- v1[1] <- max_loglik
  #
  # Starting from the MLE, we search upwards and downwards until we pass the
  # cutoff for the 100level% confidence interval
  #
  ### Upper tail ...
  xp <- rl_mle
  my_val <- max_loglik
  ii <- 1
  sol <- bingptheta_mle[-c(2, 4)]
  while (my_val > conf_line){
    xp <- xp + inc
    opt <- stats::optim(sol, bingp_neg_prof_loglik, method = "BFGS", xp = xp)
    sol <- opt$par
    ii <- ii + 1
    x2[ii] <- xp
    v2[ii] <- -opt$value
    my_val <- v2[ii]
  }
  sol_up <- sol
  ### Lower tail ...
  xp <- rl_mle
  my_val <- max_loglik
  ii <- 1
  sol <- bingptheta_mle[-c(2, 4)]
  while (my_val > conf_line){
    xp <- xp - inc
    opt <- stats::optim(sol, bingp_neg_prof_loglik, method = "BFGS", xp = xp)
    sol <- opt$par
    ii <- ii + 1
    x1[ii] <- xp
    v1[ii] <- -opt$value
    my_val <- v1[ii]
  }
  sol_low <- sol
  #
  # Find the limits of the confidence interval
  #
  prof_lik <- c(rev(v1), v2)
  ret_levs <- c(rev(x1), x2)
  # Find where the curve crosses conf_line
  temp <- diff(prof_lik - conf_line > 0)
  # Find the upper limit of the confidence interval
  loc <- which(temp == -1)
  x1 <- ret_levs[loc]
  x2 <- ret_levs[loc + 1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc + 1]
  up_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  # Find the lower limit of the confidence interval
  loc <- which(temp == 1)
  x1 <- ret_levs[loc]
  x2 <- ret_levs[loc+1]
  y1 <- prof_lik[loc]
  y2 <- prof_lik[loc+1]
  low_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  rl_prof <- c(lower = low_lim, rl_mle, upper = up_lim)
  return(list(rl_prof = rl_prof, crit = conf_line,
              for_plot = cbind(ret_levs = ret_levs, prof_loglik = prof_lik)))
}

# ============================== box_cox_deriv ============================== #

#' @keywords internal
#' @rdname lite-internal
box_cox_deriv <- function(x, lambda = 1, lambda_tol = 1 / 50,
                          poly_order = 3) {
  #
  # Computes the derivative with respect to lambda the Box-Cox
  # transformation.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   poly_order : order of Taylor series polynomial in lambda used as
  #                an approximation if abs(lambda) < lambda_tol
  #
  # Returns:
  #   A numeric vector.  The derivative with respect to lambda of
  #     (x^lambda - 1) / lambda
  #
  lnx <- log(x)
  if (abs(lambda) > lambda_tol) {
    retval <- (lambda * x ^ lambda * lnx - x ^ lambda + 1) / lambda ^ 2
  } else {
    i <- 0:poly_order
    retval <- sum(lnx ^ (i + 2) * lambda ^ i / ((i + 2) * factorial(i)))
  }
  return(retval)
}

# ========================== create_ru_list ================================= #

#' @keywords internal
#' @rdname lite-internal
make_ru_list <- function(model, trans, rotate, min_xi, max_xi) {
  #
  # Creates a list of arguments to pass to the functions ru() or ru_rcpp()
  # in the rust package to perform ratio-of-uniforms sampling from a
  # posterior density.
  #
  # Args:
  #   model     : character string specifying the extreme value model.
  #   trans     : "none", no transformation.
  #               "BC", marginal Box-Cox transformation
  #   rotate    : if TRUE rotate posterior using Cholesky decomposition of
  #               Hessian of negated log-posterior.
  #   min_xi    : the smallest xi with a non-zero posterior density
  #   max_xi    : the largest xi with a non-zero posterior density
  #
  # Returns: a list containing the inputs model, trans, rotate and
  #   d         : the dimension of the density (number of model parameters)
  #   lower     : vector of lower bounds on the arguments of logf.
  #   upper     : vector of upper bounds on the arguments of logf.
  #   var_names : the names of the variables (posterior parameters)
  #
  if (model == "gp") {
    d <- 2L
    if (trans == "none") {
      lower <- c(0, min_xi)
      upper <- c(Inf, max_xi)
    } else if (trans == "BC") {
      lower <- c(0, 0)
      upper <- c(Inf, Inf)
    } else {
      lower <- rep(-Inf, 2)
      upper <- rep(Inf, 2)
    }
    var_names <- c("sigma[u]", "xi")
  }
  if (model == "pp") {
    d <- 3L
    if (trans == "none") {
      lower <- c(-Inf, 0, min_xi)
      upper <- c(Inf, Inf, max_xi)
    } else if (trans == "BC") {
      lower <- c(-Inf, 0, 0)
      upper <- c(Inf, Inf, Inf)
    } else {
      lower <- rep(-Inf, 3)
      upper <- rep(Inf, 3)
    }
    var_names = c("mu","sigma", "xi")
  }
  return(list(d = d, lower = lower, upper = upper, var_names = var_names))
}
