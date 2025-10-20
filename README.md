
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CoxBoost

<img src="man/figures/warning.png" width = "42" /> This package is
**under maintenance so that it is restored in CRAN**

[Package website](TOADD)

<!-- badges: start -->

[![CRAN
Status](https://www.r-pkg.org/badges/version-ago/CoxBoost)](https://cran.r-project.org/package=CoxBoost)
<!-- badges: end -->

Cox-Likelihood Based Boosting for right-censored single-event survival
tasks and competing risks.

## Installation

You can install the development version of CoxBoost from
[GitHub](https://github.com/binderh/CoxBoost/) with:

``` r
# install.packages("pak")
pak::pak("binderh/CoxBoost")
```

## Example

This is a basic example which shows you how to fit a `CoxBoost` model
and predict on new data.

First, we generate some survival data with 10 informative covariates and
define the train and test sets:

``` r
library(CoxBoost)
set.seed(42)

# generate survival data
n = 200; p = 100
beta = c(rep(1,10),rep(0,p-10))
x = matrix(rnorm(n*p),n,p)
real.time = -(log(runif(n)))/(10*exp(drop(x %*% beta)))
cens.time = rexp(n,rate=1/10)
status = ifelse(real.time <= cens.time,1,0)
obs.time = ifelse(real.time <= cens.time,real.time,cens.time)

# define training and test set
train.index = 1:150
test.index = 151:200
```

We fit CoxBoost to the training data and see the model’s `summary()`:

``` r
cbfit = CoxBoost(
  time = obs.time[train.index],
  status = status[train.index],
  x = x[train.index,],
  stepno = 300, penalty = 100
)

summary(cbfit)
```

    # 300 boosting steps resulting in 67 non-zero coefficients  
    # partial log-likelihood: -353.4452 
    # 
    # parameter estimates > 0:
    #  V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V15, V20, V21, V23, V30, V31, V33, V34, V37, V38, V41, V47, V48, V49, V60, V63, V71, V73, V74, V77, V78, V80, V82, V87, V90, V91, V96, V97 
    # parameter estimates < 0:
    #  V13, V18, V22, V24, V27, V28, V32, V36, V40, V42, V44, V51, V52, V53, V56, V59, V67, V69, V70, V72, V76, V83, V84, V86, V89, V93, V99

Plot the mean partial log-likelihood for test set in every boosting
step:

``` r
step_pll = predict(
  cbfit,
  newdata = x[test.index,],
  newtime = obs.time[test.index],
  newstatus = status[test.index],
  at.step = 0:300, type = "logplik"
)

plot(step_pll)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

Names of covariates with non-zero coefficients at boosting step with
maximal test set partial log-likelihood:

``` r
print(cbfit$xnames[cbfit$coefficients[which.max(step_pll),] != 0])
```

    #  [1] "V1"  "V2"  "V3"  "V4"  "V5"  "V6"  "V7"  "V8"  "V9"  "V10" "V12" "V18"
    # [13] "V20" "V24" "V27" "V30" "V31" "V32" "V33" "V37" "V38" "V41" "V42" "V44"
    # [25] "V53" "V56" "V59" "V63" "V69" "V71" "V72" "V73" "V76" "V78" "V80" "V82"
    # [37] "V84" "V87" "V99"

We refit the `CoxBoost` model but with covariates 1 and 2 as mandatory:

``` r
cbfit_mand = CoxBoost(
  time = obs.time,
  status = status,
  x = x, 
  unpen.index = c(1,2), stepno = 100, penalty = 100
)

summary(cbfit_mand)
```

    # 100 boosting steps resulting in 45 non-zero coefficients (with 2 being mandatory) 
    # partial log-likelihood: -567.6112 
    # 
    # parameter estimates > 0:
    #  V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V20, V30, V31, V37, V38, V41, V49, V61, V63, V68, V71, V73, V78, V80, V85, V90, V97 
    # parameter estimates < 0:
    #  V24, V27, V32, V36, V45, V50, V52, V53, V59, V72, V76, V84, V88, V93, V99, V100 
    # Parameter estimates for mandatory covariates:
    #    Estimate
    # V1   0.9673
    # V2   0.9816
