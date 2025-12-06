# CoxBoost 1.5.1

* CRAN release.

# CoxBoost 1.5

* Latest GitHub release since the package was archived on CRAN on November 11th 2020.

# CoxBoost 1.4

* Added a formula interface through `iCoxBoost`
* Added generic function `coef` for extracting estimated coefficients
* Added a plot routine that provides coefficient paths
* Added support for package `parallel` (removing support for `multicore` and older R versions)
* Convergence problems for unpenalized covariates now are caught

# CoxBoost 1.3

* Added option `criterion` to allow for selection according to unpenalized scores
* Added `criterion="hpscore"` and `criterion="hscore"` for heuristic evaluation
  of only a subset of covariates in each boosting step
* Fixed a bug where results from `predict()` without `"newdata"` and `"linear.predictor"`
  in CoxBoost objects would have the wrong order (introduced in 1.2-1)
* Added missing value check for covariate matrix
* Implemented observation weights

# CoxBoost 1.2-2

* Fixed a bug in the predict function occurred when all coefficients
  were equal to zero
* Fixed bug where `estimPVal` with using only one boosting step
* `estimPVal` now also works for zero boosting steps

# CoxBoost 1.2-1

* Improved speed of the core selection routine
* Added faster code for the special case of binary covariate data
* Added an option for not returning the matrix with the score statistics
  for saving memory in applications with a huge number of covariates
* Optimized memory usage for a large number of covariates
* Covariates with standard deviation equal to zero now only are centered
* A matrix of the employed penalties know is only stored if the penalties,
  changed. Otherwise the 'element' penalty is just a vector
* Added support for `multicore` package for cross-validation and p-value estimation
* Added an option for fitting on subsets of observations
* The coefficient matrix is now stored as a sparse matrix, employing package `Matrix`
* Fixed the implementation of the p-value estimation

# CoxBoost 1.2

* Added function `estimPVal()` for permutation-based p-value estimation
* Improved the speed of the penalty updating code in PathBoost

# CoxBoost 1.1-1

* fixed bug in print method (introduced in 1.0-1) where the number of
  non-zero coefficients would be taken from a wrong boosting step

# CoxBoost 1.1

* Implemented penalty modification factors and penalty change distribution
  via a connection matrix
* Implemented estimation of models for competing risks

# CoxBoost 1.0-1

* Implemented data adaptive rule for default penalty value
* Fixed bug where output of the selected covariate would print the wrong name in 
  presence of unpenalized covariates
* Boosting now starts a step 0, i.e., also the model before updating
  any of the coefficients of the penalized covariates is considered.
  However, the unpenalized covariates will already have non-zero
  values in boosting step 0.
  This change breaks code that relies on the size of elements
  `"coefficients"`, `"linear.predictors"`, or `"Lambda"` of CoxBoost objects
* Implemented parallel evaluation of cross-validation folds, via package `snowfall`
* Speed improvements by replacing 'apply' and 'rbind', most noticeably
  for a large number of observations with a small number of covariates

# CoxBoost 1.0

* Initial public release
