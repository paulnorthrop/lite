
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lite

[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/lite?branch=main&svg=true)](https://ci.appveyor.com/project/paulnorthrop/lite)
[![R-CMD-check](https://github.com/paulnorthrop/lite/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paulnorthrop/lite/actions/workflows/R-CMD-check.yaml)
[![Coverage
Status](https://codecov.io/github/paulnorthrop/lite/coverage.svg?branch=main)](https://codecov.io/github/paulnorthrop/lite?branch=main)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/lite)](https://cran.r-project.org/package=lite)
[![Downloads
(monthly)](https://cranlogs.r-pkg.org/badges/lite?color=brightgreen)](https://cran.r-project.org/package=lite)
[![Downloads
(total)](https://cranlogs.r-pkg.org/badges/grand-total/lite?color=brightgreen)](https://cran.r-project.org/package=lite)

## Likelihood-Based Inference for Time Series Extremes

The **lite** package performs likelihood-based inference for stationary
time series extremes. The general approach follows [Fawcett and Walshaw
(2012)](https://doi.org/10.1002/env.2133). There are 3 independent parts
to the inference.

1.  A Bernoulli (*p*<sub>*u*</sub>) model for whether a given
    observation exceeds the threshold $u$.
2.  A generalised Pareto, GP (*σ*<sub>*u*</sub>,*ξ*), model for the
    marginal distribution of threshold excesses.
3.  The $K$-gaps model for the extremal index $\theta$, based on
    inter-exceedance times.

For parts 1 and 2 it is necessary to adjust the inferences because we
expect that the data will exhibit cluster dependence. This is achieved
using the methodology developed in [Chandler and Bate
(2007)](https://doi.org/10.1093/biomet/asm015) to produce a
log-likelihood that is adjusted for this dependence. This is achieved
using the [chandwich
package](https://cran.r-project.org/package=chandwich). For part 3, the
methodology described in [Süveges and Davison
(2010)](https://doi.org/10.1214/09-AOAS292) is used, implemented by the
function `kgaps` in the [exdex
package](https://cran.r-project.org/package=exdex). The (adjusted)
log-likelihoods from parts 1, 2 and 3 are combined to make inferences
about return levels.

We illustrate the main functions in `lite` using the `cheeseboro` wind
gusts data from the [exdex
package](https://cran.r-project.org/package=exdex), which contains
hourly wind gust data from each January over the 10-year period
2000-2009.

### Frequentist inference

The function `flite` makes frequentist inferences about
$(p_u, \sigma_u, \xi, \theta)$ using maximum likelihood estimation.
First, we make inferences about the model parameters.

``` r
library(lite)
cdata <- exdex::cheeseboro
# Each column of the matrix cdata corresponds to data from a different year
# flite() sets cluster automatically to correspond to column (year)
cfit <- flite(cdata, u = 45, k = 3)
summary(cfit)
#> 
#> Call:
#> flite(data = cdata, u = 45, k = 3)
#> 
#>          Estimate Std. Error
#> p[u]      0.02771   0.005988
#> sigma[u]  9.27400   2.071000
#> xi       -0.09368   0.084250
#> theta     0.24050   0.023360
```

Then, we make inferences about the 100-year return level, including 95%
confidence intervals. The argument `ny` sets the number of observations
per year, which is $31 \times 24 = 744$ for these data.

``` r
rl <- returnLevel(cfit, m = 100, level = 0.95, ny = 31 * 24)
rl
#> 
#> Call:
#> returnLevel(x = cfit, m = 100, level = 0.95, ny = 31 * 24)
#> 
#> MLE and 95% confidence limits for the 100-year return level
#> 
#> Normal interval:
#>  lower     mle   upper  
#>  69.80   88.62  107.44  
#> Profile likelihood-based interval:
#>  lower     mle   upper  
#>  75.89   88.62  125.43
```

### Bayesian inference

The function `blite` performs Bayesian inferences about
$(p_u, \sigma_u, \xi, \theta)$, based on a likelihood constructed from
the (adjusted) log-likelihoods detailed above. First, we sample from the
posterior distribution of the parameters, using the default priors.

``` r
cpost <- blite(cdata, u = 45, k = 3, ny = 31 * 24)
summary(cpost)
#> 
#> Call:
#> blite(data = cdata, u = 45, k = 3, ny = 31 * 24)
#> 
#>          Posterior mean Posterior SD
#> p[u]            0.02820     0.005942
#> sigma[u]       10.16000     2.505000
#> xi             -0.07531     0.097800
#> theta           0.24150     0.023170
```

Then, we estimate a 95% highest predictive density (HPD) interval for
the largest value $M_{100}$ to be observed over a future time period of
length $100$ years.

``` r
predict(cpost, hpd = TRUE, n_years = 100)$short
#>         lower    upper n_years level
#> [1,] 73.45208 139.7568     100    95
```

Objects returned from `flite` and `blite` have `plot` methods to
summarise graphically, respectively, log-likelihoods and posterior
distributions.

### Installation

To get the current released version from CRAN:

``` r
install.packages("lite")
```

### Vignettes

See `vignette("lite-1-frequentist", package = "lite")` and
`vignette("lite-2-bayesian", package = "lite")`.
