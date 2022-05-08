#' Bayesian threshold-based inference for time series extremes
#'
#' Performs threshold-based Bayesian inference for 3 aspects of stationary
#' time series extremes: the probability that the threshold is exceeded, the
#' marginal distribution of threshold excesses and the extent of clustering of
#' extremes, as summarised by the extremal index.
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
#' @param u A numeric scalar.  The extreme value threshold applied to the data.
#'   See \strong{Details} for information about choosing \code{u}.
#' @param cluster This argument is used to set the argument \code{cluster} to
#'   \code{\link[sandwich:vcovCL]{meatCL}}, which calculates the matrix \eqn{V}
#'   passed as the argument \code{V} to \code{\link[chandwich]{adjust_loglik}}.
#'   If \code{data} is a matrix and \code{cluster} is missing then
#'   \code{cluster} is set so that data in different columns are in different
#'   clusters.  If \code{data} is a vector and \code{cluster} is missing then
#'   cluster is set so that each observation forms its own cluster.
#'
#'   If \code{cluster} is supplied then it must have the same structure as
#'   \code{data}: if \code{data} is a matrix then \code{cluster} must be a
#'   matrix with the same dimensions as \code{data} and if  \code{data} is a
#'   vector then \code{cluster} must be a vector of the same length as
#'   \code{data}.  Each entry in \code{cluster} sets the cluster of the
#'   corresponding component of \code{data}.
#' @param k,inc_cens Arguments passed to \code{\link[exdex]{kgaps}}.
#'   \code{k} sets the value of the run parameter \eqn{K} in the \eqn{K}-gaps
#'   model for the extremal index.
#'   \code{inc_cens} determines whether contributions from right-censored
#'   inter-exceedance times are used. See \strong{Details} for information
#'   about choosing \code{k}.
#' @param ny A numeric scalar.  The (mean) number of observations per year.
#'   Setting this appropriately is important when making inferences about
#'   return levels, using \code{\link{returnLevel}}, but \code{ny} is not
#'   used by \code{flite} so it need not be supplied now.  If \code{ny} is
#'   supplied to \code{flite} then it is stored for use by
#'   \code{\link{returnLevel}}.  Alternatively, \code{ny} can be supplied in
#'   a later call to \code{\link{returnLevel}}.  If \code{ny} is supplied to
#'   both \code{flite} and \code{\link{returnLevel}} then the value supplied to
#'   \code{\link{returnLevel}} will take precedence.
#' @param gp_prior A list to specify a prior distribution for the GP parameters
#'   (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'   \eqn{\xi}), set using \code{\link[revdbayes]{set_prior}}.
#' @param b_prior A list to specify a prior distribution for the Bernoulli
#'   parameter \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{\sigma_u}}, set
#'   using \code{\link[revdbayes]{set_bin_prior}}.
#' @param theta_prior_pars A numerical vector of length 2 containing the
#'   respective values of the parameters \eqn{\alpha} and \eqn{\beta} of a
#'   Beta(\eqn{\alpha}, \eqn{\beta}) prior for the extremal index \eqn{\theta}.
#' @param n An integer scalar.  The size of posterior sample required.
#' @param type A character scalar.  Either \code{"vertical"} to adjust the
#'   independence log-likelihood vertically, or \code{"none"} for no
#'   adjustment.  Horizontal adjustment is not offered because it does not
#'   preserve the correct support of the posterior distribution.
#' @param ... Further arguments to be passed to the function
#'   \code{\link[sandwich:vcovCL]{meatCL}} in the sandwich package.
#'   In particular, the clustering adjustment argument \code{cadjust}
#'   may make a difference if the number of clusters is not large.
#' @details There are 3 independent parts to the inference, all performed using
#'   maximum likelihood estimation.
#'     \enumerate{
#'       \item{A Bernoulli(\eqn{p_u}) model for whether a given observation
#'         exceeds the threshold \eqn{u}.}
#'       \item{A generalised Pareto, GP(\eqn{\sigma_u}, \eqn{\xi}), model for
#'         the marginal distribution of threshold excesses.}
#'       \item{The \eqn{K}-gaps model for the extremal index \eqn{\theta}.}
#'     }
#'   The general approach follows Fawcett and Walshaw (2012).
#'
#'   For parts 1 and 2, inferences based on a mis-specified independence
#'   log-likelihood are adjusted to account for clustering in the data. Here,
#'   we follow Chandler and Bate (2007) to estimate adjusted log-likelihood
#'   functions for \eqn{p_u} and for (\eqn{\sigma_u}, \eqn{\xi}), with the
#'   argument \code{cluster} defining the clusters. This aspect of the
#'   calculations is performed using the \code{\link[chandwich]{adjust_loglik}}
#'   in the \code{\link[chandwich]{chandwich}} package (Northrop and Chandler,
#'   2021). The GP distribution initial fit of the GP distribution to threshold
#'   excesses is performed using the \code{\link[revdbayes]{grimshaw_gp_mle}}
#'   function in the \code{\link[revdbayes]{revdbayes}} package
#'   (Northrop, 2020).
#'
#'   In part 3, the methodology described in Suveges and Davison (2010) is
#'   implemented using the \code{\link[exdex]{exdex}} package
#'   (Northrop and Christodoulides, 2022).
#'
#'   Two tuning parameters need to be chosen: a threshold \eqn{u} and the
#'   \eqn{K}-gaps run parameter \eqn{K}.  The \code{\link[exdex]{exdex}}
#'   package has a function \code{\link[exdex]{choose_uk}} to inform this
#'   choice.
#'
#'   Each of part of the inference produces a log-likelihood function (adjusted
#'   for parts 1 and 2).  These log-likelihoods are combined (summed) to form
#'   a log-likelihood function for the parameter vector
#'   \eqn{(p_u, \sigma_u, \xi, \theta)}.  Return levels are a function of these
#'   parameters and therefore inferences for return levels can be based on
#'   this log-likelihood.
#' @return An object of class \code{c("flite", "lite", "chandwich")}.
#'   This object is a function with 2 arguments:
#'     \itemize{
#'       \item{\code{pars}, a numeric vector of length 4 to supply the value of
#'         the parameter vector \eqn{(p_u, \sigma_u, \xi, \theta)},}
#'       \item{\code{type}, a character scalar specifying the type of
#'         adjustment made to the independence log-likelihood in parts
#'         1 and 2, one of \code{"vertical"}, \code{"none"}, \code{"cholesky"},
#'         or \code{"spectral"}.  For details see Chandler and Bate (2007).
#'         The default is \code{"vertical"} for the reason given in
#'         the description of the argument \code{adj_type} in
#'         \code{\link{plot.flite}}.}
#'     }
#'  The object also has the attributes \code{"Bernoulli"}, \code{"gp"},
#'  \code{"theta"}, which provide the fitted model objects returned from
#'  \code{\link[chandwich]{adjust_loglik}} (for \code{"Bernoulli"} and
#'  \code{"gp"}) and \code{\link[exdex]{kgaps}} (for \code{"theta"}).
#'
#'  The named input arguments are returned in a list as the attribute
#'  \code{inputs}.  If \code{ny} was not supplied then its value is \code{NA}.
#'
#'   Objects inheriting from class \code{"flite"} have \code{coef},
#'   \code{logLik}, \code{nobs}, \code{plot}, \code{summary} and \code{vcov}
#'   methods.  See \code{\link{fliteMethods}}.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered.
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \doi{10.1093/biomet/asm015}
#' @references Fawcett, L. and Walshaw, D. (2012), Estimating return levels
#'   from serially dependent extremes. \emph{Environmetrics}, \strong{23},
#'   272-283. \doi{10.1002/env.2133}
#' @references Northrop, P. J. and Chandler, R. E. (2021).
#'   chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment. R package
#'   version 1.1.5. \url{https://CRAN.R-project.org/package=chandwich}.
#' @references Northrop, P. J. and Christodoulides, C. (2022). exdex:
#' Estimation of the Extremal Index. R package version 1.1.1.
#' \url{https://CRAN.R-project.org/package=exdex/}.
#' @references  Northrop, P. J. (2020). revdbayes: Ratio-of-Uniforms Sampling
#' for Bayesian Extreme Value Analysis. R package version 1.3.9.
#' \url{https://paulnorthrop.github.io/revdbayes/}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \doi{10.1214/09-AOAS292}
#' @seealso \code{\link{fliteMethods}}, including plotting (adjusted)
#'   log-likelihoods for \eqn{(p_u, \sigma_u, \xi, \theta)}.
#' @seealso \code{\link{Bernoulli}} for maximum likelihood inference for the
#'   Bernoulli distribution.
#' @seealso \code{\link{generalisedPareto}} for maximum likelihood inference
#'   for the generalised Pareto distribution.
#' @seealso \code{\link[exdex]{kgaps}} for maximum likelihood inference from
#'   the \eqn{K}-gaps model for the extremal index.
#' @seealso \code{\link[exdex]{choose_uk}} to inform the choice of the
#'   threshold \eqn{u} and run parameter \eqn{K}.
#' @examples
#' ### Cheeseboro wind gusts
#'
#' # Make inferences
#' cdata <- exdex::cheeseboro
#' # Each column of the matrix cdata corresponds to data from a different year
#' # blite() sets cluster automatically to correspond to column (year)
#' cpost <- blite(cdata, u = 45, k = 3)
#' summary(cpost)
#' @export
blite <- function(data, u, cluster, k = 1, inc_cens = TRUE, ny,
                  gp_prior = revdbayes::set_prior(prior = "mdi", model = "gp"),
                  b_prior = revdbayes::set_bin_prior(prior = "jeffreys"),
                  theta_prior_pars = c(1, 1), n = 1000,
                  type = c("vertical", "none"), ...) {
  type <- match.arg(type)
  #
  # 1. Call flite() to create an object that contains the functions and
  #    information needed to sample from the posterior distribution of
  #    (p[u], sigma[u], xi, theta)
  #
  x <- flite(data = data, u = u, cluster = cluster, k = k, inc_cens = inc_cens,
             ny = ny, ...)
  # Save cluster and ny so that they can be returned in the "inputs" attribute
  cluster <- attr(x, "inputs")$cluster
  ny <- attr(x, "inputs")$ny
  #
  # 2. Sample from a posterior for the GP parameters (sigma[u], xi) based on
  #    the adjusted log-likelihood in attr(x, "gp") and the prior in gp_prior
  #
  # Extract the "chandwich" object
  adj_gp_loglik <- attr(x, "gp")
  # Extract min_xi and max_xi from prior (if supplied)
  min_xi <- ifelse(is.null(gp_prior$min_xi), -Inf, gp_prior$min_xi)
  max_xi <- ifelse(is.null(gp_prior$max_xi), +Inf, gp_prior$max_xi)
  # Create a function to evaluate the (adjusted) GP posterior density
  gp_logpost <- function(pars) {
    loglik <- do.call(adj_gp_loglik, list(x = pars, type = type))
    logprior <- do.call(gp_prior$prior, c(list(pars), gp_prior[-1]))
    return(loglik + logprior)
  }
  # Use the MLEs as initial estimates
  gp_init <- coef(x)[c("sigma[u]", "xi")]
  # Create list of objects to send to function rust::ru()
  fr <- make_ru_list(model = "gp", trans = "none", rotate = TRUE,
                     min_xi = min_xi, max_xi = max_xi)
  for_ru <- c(list(logf = gp_logpost), fr, list(init = gp_init, n = n))
  gp_posterior <- do.call(rust::ru, for_ru)
  #
  # 3. Sample from a posterior for the Bernoulli parameter p[u] based on the
  #    adjusted log-likelihood in attr(x, "Bernoulli") and the prior in b_prior
  #
  # Extract the "chandwich" object
  adj_b_loglik <- attr(x, "Bernoulli")
  # We can use revdbayes::binpost() to sample from a posterior distribution
  # based on the (vertically) adjusted log-likelihood.  The effect of the
  # vertical adjustment is to multiply the numbers of successes n_s (threshold
  # exceedances) and the number of failures n_f (threshold nn-exceedances) by
  # the ratio of the estimate of HA and HI.
  b_fit <- attr(adj_b_loglik, "original_fit")
  n_s <- b_fit$n1
  n_f <- b_fit$n0
  if (type == "none") {
    HA_over_HI <- 1
  } else {
    HA_over_HI <- attr(adj_b_loglik, "HA") / attr(adj_b_loglik, "HI")
  }
  # Adjust n_s and n_f
  n_s <- HA_over_HI * n_s
  n_f <- HA_over_HI * n_f
  # Create list of arguments for revdbayes::binpost
  # If the prior is user-supplied then we the default param = "logit", so
  # we sample on the logit scale and back-transform to the p[u]-scale
  ds_bin <- list(n_raw = n_s + n_f, m = n_s)
  b_posterior <- revdbayes::binpost(n = n, prior = b_prior, ds_bin = ds_bin)
  # Make the simulated values a 1-column matrix
  dim(b_posterior$bin_sim_vals) <- c(length(b_posterior$bin_sim_vals), 1)
  colnames(b_posterior$bin_sim_vals) <- "p[u]"
  #
  # 4. Sample from a posterior for the K-gaps parameter theta using
  #    revdbayes::kgaps_post()
  #
  alpha <- theta_prior_pars[1]
  beta <- theta_prior_pars[2]
  theta_posterior <- kgaps_post(data = data, thresh = u, k = k, n = n,
                                inc_cens = inc_cens, alpha = alpha,
                                beta = beta)
  # Make the simulated values a 1-column matrix
  dim(theta_posterior$sim_vals) <- c(length(theta_posterior$sim_vals), 1)
  colnames(theta_posterior$sim_vals) <- "theta"
  #
  # 5. Combine 2,3,4 to form a posterior sample for (p[u], sigma[u], xi, theta)
  #
  posterior_sample <- matrix(NA, nrow = n, ncol = 4)
  posterior_sample[, 1] <- b_posterior$bin_sim_vals
  posterior_sample[, 2:3] <- gp_posterior$sim_vals
  posterior_sample[, 4] <- theta_posterior$sim_vals
  names(posterior_sample) <- c("p[u]", "sigma[u]", "xi", "theta")
  # Set the class of the return object and add attributes
  class(posterior_sample) <- c("blite", "lite", "chandwich")
  attr(posterior_sample, "flite_object") <- x
  attr(posterior_sample, "Bernoulli") <- b_posterior
  attr(posterior_sample, "gp") <- gp_posterior
  attr(posterior_sample, "theta") <- theta_posterior
  attr(posterior_sample, "call") <- match.call(expand.dots = TRUE)
  inputs <- list(data = data, u = u, cluster = cluster, k = k,
                 inc_cens = inc_cens, ny = ny, gp_prior = gp_prior,
                 b_prior = b_prior, theta_prior_pars = theta_prior_pars,
                 n = n, type = type)
  attr(posterior_sample, "inputs") <- inputs
  return(posterior_sample)
}
