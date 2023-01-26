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
#'   Setting this appropriately is important when making predictive inferences
#'   using \code{\link{predict.blite}}, but \code{ny} is not used by
#'   \code{blite} so it need not be supplied now.  If \code{ny} is supplied to
#'   \code{blite} then it is stored for use by \code{\link{predict.blite}}.
#'   Alternatively, \code{ny} can be supplied in a later call to
#'   \code{\link{predict.blite}}.  If \code{ny} is supplied to
#'   both \code{blite} and \code{\link{predict.blite}} then the value supplied
#'   to \code{\link{predict.blite}} will take precedence, with no warning
#'   given.
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
#' @details See \code{\link{flite}} for details of the (adjusted) likelihoods
#'   on which these Bayesian inferences are based.
#'
#'   The likelihood is based on a model for 3 independent aspects.
#'     \enumerate{
#'       \item{A Bernoulli(\ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}}) model
#'         for whether a given observation exceeds the threshold \eqn{u}.}
#'       \item{A generalised Pareto,
#'         GP(\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'         \eqn{\xi}), model for the marginal distribution of threshold
#'         excesses.}
#'       \item{The \eqn{K}-gaps model for the extremal index \eqn{\theta}.}
#'     }
#'   The general approach follows Fawcett and Walshaw (2012).
#'
#'   The contributions to the likelihood for
#'   \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}} and
#'   (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}}, \eqn{\xi})
#'   are based on the vertically-adjusted likelihoods described in
#'   \code{\link{flite}}.  This is an example of Bayesian inference using a
#'   composite likelihood Ribatet et al (2012). Priors for
#'   \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}}
#'   (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}}, \eqn{\xi})
#'   and \eqn{\theta} are set using the arguments \code{gp_prior},
#'   \code{b_prior} and \code{theta_prior_pars}.
#'   Currently, only priors where
#'   \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}}
#'   (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}}, \eqn{\xi})
#'   and \eqn{\theta} are independent a priori are allowed.
#'
#'   Two tuning parameters need to be chosen: a threshold \eqn{u} and the
#'   \eqn{K}-gaps run parameter \eqn{K}.  The \code{\link[exdex]{exdex}}
#'   package has a function \code{\link[exdex]{choose_uk}} to inform this
#'   choice.
#'
#'   Random samples are simulated from the posteriors for
#'   \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}} and
#'   (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}}, \eqn{\xi})
#'   (using \code{\link[rust]{ru}}) and \eqn{\theta} (using
#'   \code{\link[revdbayes]{kgaps_post}}).
#' @return An object of class \code{c("blite", "lite", "chandwich")}.
#'   This object is an \code{n} \eqn{\times 4}{x 4} matrix containing the
#'   posterior samples, with column names
#'   \code{c("p[u]", "sigma[u]", "xi", "theta")}.
#'
#'  The object also has the attributes \code{"Bernoulli"}, \code{"gp"},
#'  \code{"theta"}, which provide the fitted model objects returned from
#'  \code{\link[chandwich]{adjust_loglik}} (for \code{"Bernoulli"} and
#'  \code{"gp"}) and \code{\link[exdex]{kgaps}} (for \code{"theta"}).
#'  The named input arguments are returned in a list as the attribute
#'  \code{inputs}.  If \code{ny} was not supplied then its value is \code{NA}.
#'  The call to \code{blite} is provided in the attribute \code{"call"}.
#'  A call to \code{\link{flite}} is used to create adjusted log-likelihoods
#'  for \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}} and
#'  (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}}, \eqn{\xi}).
#'  The object returned from the call is provided as the attribute
#'  \code{"flite_object"}.
#'
#'  Objects inheriting from class \code{"blite"} have \code{coef},
#'  \code{nobs}, \code{plot}, \code{summary} and \code{vcov}
#'  methods.  See \code{\link{bliteMethods}}.
#' @references Fawcett, L. and Walshaw, D. (2012), Estimating return levels
#'   from serially dependent extremes. \emph{Environmetrics}, \strong{23},
#'   272-283. \doi{10.1002/env.2133}
#' @references Ribatet, M., Cooley, D., & Davison, A. C. (2012). Bayesian
#'   inference from composite likelihoods, with an application to spatial
#'   extremes. \emph{Statistica Sinica}, \strong{22}(2), 813-845.
#' @seealso \code{\link{bliteMethods}}, including plotting the posterior
#'   samples.
#' @seealso \code{\link{flite}} for frequentist threshold-based inference
#'   for time series extremes.
#' @seealso \code{\link[exdex]{choose_uk}} to inform the choice of the
#'   threshold \eqn{u} and run parameter \eqn{K}.
#' @examples
#' ### Cheeseboro wind gusts
#'
#' cdata <- exdex::cheeseboro
#' # Each column of the matrix cdata corresponds to data from a different year
#' # blite() sets cluster automatically to correspond to column (year)
#' cpost <- blite(cdata, u = 45, k = 3)
#' summary(cpost)
#' plot(cpost)
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
  colnames(posterior_sample) <- c("p[u]", "sigma[u]", "xi", "theta")
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
