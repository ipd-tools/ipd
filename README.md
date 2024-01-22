
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
particular, the statistical challenges encountered when drawing include:
(1) understanding the relationship between predicted outcomes and their
true unobserved counterparts, (2) quantifying the robustness of the
AI/ML models to resampling or uncertainty about the training data, and
(3) appropriately propagating both bias and uncertainty from predictions
into downstream inferential tasks.

Recently, several works have proposed methods to address this general
problem of IPD. These include by [Wang et al.,
2020](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117), and by
[Angelopoulos et al.,
2023a](https://www.science.org/doi/10.1126/science.adi6000) and
[Angelopoulos et al., 2023b](https://arxiv.org/abs/2311.01453),
respectively, and by [Miao et al.,
2023](https://arxiv.org/abs/2311.14220). These methods have been
developed in quick succession in response to the ever-growing practice
of using predicted data directly to conduct statistical inference. To
allow researchers and practitioners the ability to fully utilize these
state-of-the-art methods, we have developed , a comprehensive and
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
methods.

(i). Assume that we have access to a well-performing and fairly accurate
AI/ML/DL algorithm $f_{\text{x}}(\cdot)$ that can predict our outcome of
interest $Y$.

(ii). Next, consider that we have 2 data sets, a labeled data set (which
we will call the **test set** $(X_{te}, Y_{te})$), and an unlabeled data
set (which we call the **validation set** $(X_{val)}$). Typically the
the labeled/test set is considerably smaller in size compared to the
unlabeled/validation set. Here we will consider them to be equal for
brevity.

- We consider the regressors, $X = (X_1, X_2, X_3, X_4)$ and let $Y$ be
  a scalar.
- The true data generating mechanism is
  $Y = \beta_1X_1 + \beta_2 X_2 + \beta_3 \ g(X_3) + \beta_4 \ g(X_4) + \epsilon,$
  where, $\epsilon = N(0, 1)$ and $g(\cdot)$ refers to some smoother
  function.
- We specify, $(\beta_1, \beta_2, \beta_3, \beta_4) = (1, 0.5, 3, 4)$.

(iii). Our interest is in performing inference on $H_0: \beta_1^* = 0$
vs $H_1: \beta_1^* \ne 0$. That is, our inference model is,

$$
Y_{val} = \beta_0^* + \beta_1^* X_{val} + \epsilon^*,
$$

where $\epsilon^* = N(0, 1)$.

(iv). However, we do not observe $Y_{val}$. We instead only have access
to the predicted $\hat Y_{val} = f_{\text{x}}(X_{val})$.

We will now obtain the estimates from each method below:

``` r
## Load the library
library(IPD)

## generate an example data set consisting of training, labeled, unlabeled data
set.seed(2023)
dat <- simdat(n = c(300, 300, 300), beta1 = 1)
```

#### 1.1) Analytic correction method from Wang et al. (2020)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1

IPD::postpi_analytic_ols(rel_form, inf_form, dat = dat)
#> $est
#> [1] 11.561038  0.955248
#> 
#> $se
#> [1] 0.2401493 0.1706713
```

#### 1.2) Bootstrap method from Wang et al. (2020)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1
## 3. we also need to specify the number of bootstraps 
nboot <- 200

IPD::postpi_boot_ols(rel_form, inf_form, dat = dat, nboot)
#> $est
#> [1] 11.5462734  0.9444002
#> 
#> $se_par
#> [1] 0.2394565 0.1706623
#> 
#> $se_npar
#> [1] 0.2364742 0.1807781
# the function returns both parametric (par) and 
# non-parametric (npar) estimate of std.error (se)
```

#### 2. Prediction-powered inference method (Angelopoulos et al., 2023)

``` r
form <- Y - Yhat ~ X1 ## formula

# Labeled data set
X_l <- model.matrix(form, data = dat[dat$set == "tst",]) 
Y_l <- dat[dat$set == "tst", all.vars(form)[1]] |> matrix(ncol = 1)
f_l <- dat[dat$set == "tst", all.vars(form)[2]] |> matrix(ncol = 1)

# Unlabeled data set
X_u <- model.matrix(form, data = dat[dat$set == "val",])
f_u <- dat[dat$set == "val", all.vars(form)[2]] |> matrix(ncol = 1)

n <- nrow(X_l)
p <- ncol(X_l)
N <- nrow(X_u)
IPD::ppi_ols(X_l, Y_l, f_l, X_u, f_u, n, p, N)
#> $est
#> [1] 11.697627  0.817898
#> 
#> $se
#> [1] 0.2357630 0.1790149
```

#### 3. PPI++ (Angelopoulos et al., 2023)

``` r
form <- Y - Yhat ~ X1 ## formula

# Labeled data set
X_l <- model.matrix(form, data = dat[dat$set == "tst",])
Y_l <- dat[dat$set == "tst", all.vars(form)[1]] |> matrix(ncol = 1)
f_l <- dat[dat$set == "tst", all.vars(form)[2]] |> matrix(ncol = 1)

# Unlabeled data set
X_u <- model.matrix(form, data = dat[dat$set == "val",])
f_u <- dat[dat$set == "val", all.vars(form)[2]] |> matrix(ncol = 1)

n <- nrow(X_l)
p <- ncol(X_l)
N <- nrow(X_u)

IPD::ppi_plusplus_ols(X_l, Y_l, f_l, X_u, f_u, n, p, N)
#> $est
#> X(Intercept)          XX1 
#>   11.7737782    0.8460529 
#> 
#> $se
#> [1] 0.2200044 0.1615084
```

#### 4. Assumption-lean and data-adaptive Post-Prediction Inference (POP-Inf) (Miao et al., 2023)

``` r
## rectifier formula
rec_form <- Y - Yhat ~ X1

## inference model 
inf_form <- Yhat ~ X1

IPD::popinf_ols(rec_form, inf_form, dat = dat)
#> j: 6.474096e-09j: 1
#> $est
#>                   [,1]
#> (Intercept) 11.8378384
#> X1           0.8697373
#> 
#> $se
#> (Intercept)          X1 
#>   0.2465003   0.1804958 
#> 
#> $w
#> [1] 6.474096e-09 1.000000e+00
```

#### 5. Multiple-imputation method from Leek et al., (2023)

``` r
# Requires the specification of 
## 1. relationship model between observed y and predicted Y-hat 
rel_form <- Y ~ Yhat  ## we consider a basic linear function
## 2. inference model
inf_form <- Yhat ~ X1
m <- 100

IPD::postpi_mi_ols(rel_form, inf_form, dat = dat, m)
#> $est
#> [1] 11.5855979  0.9351817
#> 
#> $se
#> [1] 0.3228812 0.2336300
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
