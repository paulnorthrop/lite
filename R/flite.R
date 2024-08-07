#' Frequentist threshold-based inference for time series extremes
#'
#' Performs threshold-based frequentist inference for 3 aspects of stationary
#' time series extremes: the probability that the threshold is exceeded, the
#' marginal distribution of threshold excesses and the extent of clustering of
#' extremes, as summarised by the extremal index.
#'
#' @param data A numeric vector or numeric matrix of raw data.  If \code{data}
#'   is a matrix then the log-likelihood is constructed as the sum of
#'   (independent) contributions from different columns. A common situation is
#'   where each column relates to a different year.
#'
#'   If \code{data} contains missing values then
#'   \code{\link[exdex]{split_by_NAs}} is used to divide the data further into
#'   sequences of non-missing values, stored in different columns in a matrix.
#'   Again, the log-likelihood is constructed as a sum of contributions from
#'   different columns.
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
#'   \code{\link{returnLevel}} will take precedence, with no warning given.
#' @param ... Further arguments to be passed to the function
#'   \code{\link[sandwich:vcovCL]{meatCL}} in the sandwich package.
#'   In particular, the clustering adjustment argument \code{cadjust}
#'   may make a difference if the number of clusters is not large.
#' @details There are 3 independent parts to the inference, all performed using
#'   maximum likelihood estimation.
#'     \enumerate{
#'       \item{A Bernoulli(\ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}})
#'         model for whether a given observation exceeds the threshold
#'         \eqn{u}.}
#'       \item{A generalised Pareto,
#'         GP(\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'         \eqn{\xi}), model for the marginal distribution of threshold
#'         excesses.}
#'       \item{The \eqn{K}-gaps model for the extremal index \eqn{\theta}.}
#'     }
#'   The general approach follows Fawcett and Walshaw (2012).
#'
#'   For parts 1 and 2, inferences based on a mis-specified independence
#'   log-likelihood are adjusted to account for clustering in the data. Here,
#'   we follow Chandler and Bate (2007) to estimate adjusted log-likelihood
#'   functions for \ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}}
#'   and for (\ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'   \eqn{\xi}), with the
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
#'   Each part of the inference produces a log-likelihood function (adjusted
#'   for parts 1 and 2).  These log-likelihoods are combined (summed) to form
#'   a log-likelihood function for the parameter vector
#'   (\ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'   \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'   \eqn{\xi}, \eqn{\theta}).  Return levels are a function of these
#'   parameters and therefore inferences for return levels can be based on
#'   this log-likelihood.
#' @return An object of class \code{c("flite", "lite", "chandwich")}.
#'   This object is a function with 2 arguments:
#'     \itemize{
#'       \item{\code{pars}, a numeric vector of length 4 to supply the value of
#'         the parameter vector
#'         (\ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'         \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'         \eqn{\xi}, \eqn{\theta}),}
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
#'  The named input arguments are returned in a list as the attribute
#'  \code{inputs}.  If \code{ny} was not supplied then its value is \code{NA}.
#'  The call to \code{flite} is provided in the attribute \code{"call"}.
#'
#'  Objects inheriting from class \code{"flite"} have \code{coef},
#'  \code{logLik}, \code{nobs}, \code{plot}, \code{summary}, \code{vcov}
#'  and \code{confint} methods.  See \code{\link{fliteMethods}}.
#'
#'  \code{\link{returnLevel}} can be used to make frequentist inferences about
#'   return levels.
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
#'   log-likelihoods for
#'   (\ifelse{html}{\eqn{p}\out{<sub>u</sub>}}{\eqn{p_u}},
#'   \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}},
#'   \eqn{\xi}, \eqn{\theta}).
#' @seealso \code{\link{returnLevel}} to make frequentist inferences about
#'   return levels.
#' @seealso \code{\link{blite}} for Bayesian threshold-based inference
#'   for time series extremes.
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
#' # flite() sets cluster automatically to correspond to column (year)
#' cfit <- flite(cdata, u = 45, k = 3)
#' summary(cfit)
#'
#' # 2 ways to find the maximised log-likelihood value
#' cfit(coef(cfit))
#' logLik(cfit)
#'
#' # Plots of (adjusted) log-likelihoods
#' plot(cfit)
#' plot(cfit, which = "gp")
#'
#' ## Confidence intervals
#' # Based on an adjusted profile log-likelihood
#' confint(cfit)
#' # Symmetric intervals based on large sample normality
#' confint(cfit, profile = FALSE)
#' @export
flite <- function(data, u, cluster, k = 1, inc_cens = TRUE, ny, ...) {
  # Check that the threshold does not lie above all the data
  if (u >= max(data, na.rm = TRUE)) {
    stop("'u' must be less than 'max(data, na.rm = TRUE)'")
  }
  # Check that the threshold does not lie below all the data
  if (u < min(data, na.rm = TRUE)) {
    stop("'u' must be greater than 'min(data, na.rm = TRUE)'")
  }
  # If ny is not supplied then store NA
  if (missing(ny)) {
    ny <- NA
  }
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
      if (!all(dim(cluster) == dim(data))) {
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
  # 2. Fit GP(sigma[u], xi) to excesses of u using an independence
  #    log-likelihood and adjust the inferences using the chandwich package
  #
  # The returned object has class "GP"
  gp_fit <- fitGP(data_vector, u)
  aloglik_gp <- adjust_object(gp_fit, cluster = gp_cluster, ...)
  #
  # 3. Fit Bernoulli(p[u]) to indicators of threshold exceedance and adjust
  #    the inferences using chandwich::adjust_loglik()
  #
  # The returned object has class "Bernoulli"
  bernoulli_fit <- fitBernoulli(data_vector > u)
  aloglik_bernoulli <- adjust_object(bernoulli_fit, cluster = cluster, ...)
  #
  # 4. Make inferences about the extremal index theta using the K-gaps model
  #    via the exdex::kgaps()
  #
  # exdex::kgaps() accept the original data matrix and omit missings itself
  theta_fit <- exdex::kgaps(data = data, u = u, k = k, inc_cens = inc_cens)
  loglik_theta <- function(tval) {
    theta_list <- c(list(theta = tval), theta_fit$ss)
    return(do.call(kgaps_loglik, theta_list))
  }
  #
  # 5. Put GP, binomial and K-gaps inferences together to return an adjusted
  #    log-likelihood for the parameters (p[u], sigma[u], xi, theta)
  #
  bernoulli_gp_theta_loglik <- function(pars, type = c("vertical", "cholesky",
                                                       "spectral", "none")) {
    type <- match.arg(type)
    bernoulli_loglik <- aloglik_bernoulli(pars[1], type = type)
    gp_loglik <- aloglik_gp(pars[2:3], type = type)
    theta_loglik <- loglik_theta(pars[4])
    val <- bernoulli_loglik + gp_loglik + theta_loglik
    names(val) <- NULL
    return(val)
  }
  # Set the class of the return object and add attributes
  class(bernoulli_gp_theta_loglik) <- c("flite", "lite", "chandwich")
  attr(bernoulli_gp_theta_loglik, "Bernoulli") <- aloglik_bernoulli
  attr(bernoulli_gp_theta_loglik, "gp") <- aloglik_gp
  attr(bernoulli_gp_theta_loglik, "theta") <- theta_fit
  attr(bernoulli_gp_theta_loglik, "call") <- match.call(expand.dots = TRUE)
  inputs <- list(data = data, u = u, cluster = cluster, k = k,
                 inc_cens = inc_cens, ny = ny)
  attr(bernoulli_gp_theta_loglik, "inputs") <- inputs
  return(bernoulli_gp_theta_loglik)
}
