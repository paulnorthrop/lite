# lite 1.1.1

## Bug fixes and minor improvements

* Fixed 2 \link{} targets with missing (exdex) package anchors in Rd files for `blite()` and `flite()`

# lite 1.1.0

## New features

* The new function `blite()` performs Bayesian threshold-based inference for time series extremes.  It is a Bayesian version of the existing function `flite()`, which performs frequentist inference.  

* Objects returned from `blite()` have a predict S3 method `predict.blite()` based on the `predict.evpost()` method from the `revdbayes`package.  It provides predictive inferences for the largest value observed in *N* years.

* Objects of class `flite` returned from `flite()` now have a `confint` method.

## Bug fixes and minor improvements

* In `flite()`, the arguments `k` and `inc_cens` were not passed to `exdex::kgaps()`.  This has been corrected.

* In the (unexported, internal) function `bingp_rl_CI()` an error is triggered if the return level requested is lower than the threshold used to fit the model.

* In the (unexported, internal) function `bingp_rl_prof()`, which calculates a confidence interval for a return level based on a profile log-likelihood, a check is made on the value `p` to be passed to `revdbayes::qgp()` to check that it is in [0, 1].
