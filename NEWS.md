# lite 1.0.0.9000

## New features

* Objects of class `flite` returned from `flite()` now have a `confint` method.

* The new function `blite()` performs Bayesian threshold-based inference for time series extremes.  It is a Bayesian version of the existing function `flite()`, which performs frequentist inference.  

* Objects returned from `blite()` have a predict S3 method `predict.blite()` based on the `predict.evpost()` method from the `revdbayes`package.  It provides predictive inferences for the largest value observed in *N* years.

## Bug fixes and minor improvements
