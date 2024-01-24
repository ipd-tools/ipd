
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IPD: Inference on Predicted Data

<!-- badges: start -->

[![R-CMD-check](https://github.com/awanafiaz/IPD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/awanafiaz/IPD/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## <img src="man/figures/IPD_LOGO.png" align="right" height="200" style="float:right; height:200px;">

With the rapid advancement of artificial intelligence and machine
learning (AI/ML) algorithms, and owing to financial and domain-specific
constraints, researchers from a wide range of disciplines are now
increasingly using predictions from pre-trained algorithms as outcome
variables in statistical analyses. However, reifying
algorithmically-derived values as measured outcomes may lead to
potentially biased estimates and anti-conservative inference ([Wang et
al., 2020](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117)). In
particular, the statistical challenges encountered when drawing
*inference on predicted data (IPD)* include: (1) understanding the
relationship between predicted outcomes and their true unobserved
counterparts, (2) quantifying the robustness of the AI/ML models to
resampling or uncertainty about the training data, and (3) appropriately
propagating both bias and uncertainty from predictions into downstream
inferential tasks.

Recently, several works have proposed methods to address this general
problem of IPD. These include *post-prediction inference (PostPI)* by
[Wang et al.,
2020](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117),
*prediction-powered inference (PPI)* and *PPI++* by [Angelopoulos et
al., 2023a](https://www.science.org/doi/10.1126/science.adi6000) and
[Angelopoulos et al., 2023b](https://arxiv.org/abs/2311.01453),
respectively, and *assumption-lean and data-adaptive post-prediction
inference (POP-Inf)* by [Miao et al.,
2023](https://arxiv.org/abs/2311.14220). These methods have been
developed in quick succession in response to the ever-growing practice
of using predicted data directly to conduct statistical inference. To
allow researchers and practitioners the ability to fully utilize these
state-of-the-art methods, we have developed `IPD`, a comprehensive and
open-source software package which implement these existing methods
under the umbrella of IPD.

To make the utilization of the package convenient for users, we provide
guidance on installation and use of the package and its functions in the
following:

## Installation

You can install the development version of `IPD` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")   ## If devtools is not already installed

devtools::install_github("awanafiaz/IPD")
```

## Usage

We provide a simple example to demonstrate the basic use of the
functions included in `IPD`. We build the premise in the following
manner to build an unifying example to be used across all available
methods. To start, we assume the user must have a dataset,
$D = L \cup U$, consisting of

- $L = \{(Y_i, f_i, \boldsymbol{X}_i, R_i);\ i = 1, \ldots, n\}$
  *labeled* samples of the observed outcome, $Y$, predicted outcome,
  $f$, and features of interest, $\boldsymbol{X}$

- $U = \{(f_i, \boldsymbol{X}_i, R_i);\ i = n + 1, \ldots, n + N\}$
  *unlabeled* samples, where the outcome is not observed.

- Here, $R$ is a variable indicating whether the $i$th observation is
  considered *labeled*, where $R_i = 1$ for the $n$ *labeled*
  observations and $R_i = 0$ for the $N$ *unlabeled* observations.

As an illustrative example, we consider simulated data, where

- $X = (X_1, X_2, X_3, X_4)$

- $Y = \beta_1X_1 + \beta_2 X_2 + \beta_3 \ g(X_3) + \beta_4 \ g(X_4) + \epsilon,$
  where, $\epsilon = N(0, 1)$ and $g(\cdot)$ refers to some smoother
  function

- We specify, $(\beta_1, \beta_2, \beta_3, \beta_4) = (1, 0.5, 3, 4)$.

Our interest is in performing inference on $H_0: \beta_1^* = 0$ vs
$H_1: \beta_1^* \ne 0$

``` r

#- Load the Library

library(IPD)

#- Generate Example Data

set.seed(12345)

dat <- simdat(n = c(300, 300, 300), beta1 = 1)
```

#### 1.1 Post-Prediction Inference - Analytic Correction (PostPI Analytic; Wang et al., 2020)

``` r

# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 

rel_form <- Y ~ Yhat  ## we consider a basic linear function

## 2. inference model

inf_form <- Yhat ~ X1

IPD::postpi_analytic_ols(rel_form, inf_form, dat = dat)
```

#### 1.2 Post-Prediction Inference - Bootstrap (PostPI Bootstrap; Wang et al., 2020)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1
## 3. we also need to specify the number of bootstraps 
nboot <- 200

IPD::postpi_boot_ols(rel_form, inf_form, dat = dat, nboot)
# the function returns both parametric (par) and 
# non-parametric (npar) estimate of std.error (se)
```

#### 2. Prediction-Powered Inference (PPI; Angelopoulos et al., 2023)

``` r
form <- Y - Yhat ~ X1 ## formula

# Labeled data set
X_l <- model.matrix(form, data = dat[dat$set == "labeled",]) 
Y_l <- dat[dat$set == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
f_l <- dat[dat$set == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)

# Unlabeled data set
X_u <- model.matrix(form, data = dat[dat$set == "unlabeled",])
f_u <- dat[dat$set == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)

n <- nrow(X_l)
p <- ncol(X_l)
N <- nrow(X_u)

IPD::ppi_ols(X_l, Y_l, f_l, X_u, f_u)
```

#### 3. PPI++ (Angelopoulos et al., 2023)

``` r
form <- Y - Yhat ~ X1 ## formula

# Labeled data set
X_l <- model.matrix(form, data = dat[dat$set == "labeled",])
Y_l <- dat[dat$set == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
f_l <- dat[dat$set == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)

# Unlabeled data set
X_u <- model.matrix(form, data = dat[dat$set == "unlabeled",])
f_u <- dat[dat$set == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)

n <- nrow(X_l)
p <- ncol(X_l)
N <- nrow(X_u)

IPD::ppi_plusplus_ols(X_l, Y_l, f_l, X_u, f_u)
```

#### 4. Assumption-Lean and Data-Adaptive Post-Prediction Inference (POP-Inf; Miao et al., 2023)

``` r
## rectifier formula
rec_form <- Y - Yhat ~ X1

## inference model 
inf_form <- Yhat ~ X1

IPD::popinf_ols(rec_form, inf_form, dat = dat)
```

#### 5. Multiple-Imputation (PostPI-MI; Leek et al., 2023)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1
m <- 100

IPD::postpi_mi_ols(rel_form, inf_form, dat = dat, m)
```

## Vignette

For more advanced users and researchers, we provide more use cases and
examples in the package `vignettes`.

``` r
vignette("IPD")
```

This provides an extensive tutorial on `IPD` and discusses
method-specific usage in details.

## Feedback

For questions and comments or any other feedback, please contact the
developers.
