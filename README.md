
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IPD: Inference on Predicted Data

<!-- badges: start -->

[![R-CMD-check](https://github.com/awanafiaz/IPD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/awanafiaz/IPD/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## <img src="man/figures/IPD_LOGO.png" align="right" height="200" style="float:right; height:200px;">

With the rapid advancement of artificial intelligence and machine
learning (AI/ML), researchers from a wide range of disciplines
increasingly use predictions from pre-trained algorithms as outcome
variables in statistical analyses. However, reifying
algorithmically-derived values as measured outcomes may lead to biased
estimates and anti-conservative inference ([Wang et al.,
2020](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117)). The
statistical challenges encountered when drawing inference on predicted
data (IPD) include:

1.  Understanding the relationship between predicted outcomes and their
    true unobserved counterparts
2.  Quantifying the robustness of the AI/ML models to resampling or
    uncertainty about the training data
3.  Appropriately propagating both bias and uncertainty from predictions
    into downstream inferential tasks

Several works have proposed methods for IPD, including post-prediction
inference (PostPI) by [Wang et al.,
2020](https://www.pnas.org/doi/suppl/10.1073/pnas.2001238117),
prediction-powered inference (PPI) and PPI++ by [Angelopoulos et al.,
2023a](https://www.science.org/doi/10.1126/science.adi6000) and
[Angelopoulos et al., 2023b](https://arxiv.org/abs/2311.01453), and
assumption-lean and data-adaptive post-prediction inference (POP-Inf) by
[Miao et al., 2023](https://arxiv.org/abs/2311.14220). To enable
researchers and practitioners interested in these state-of-the-art
methods, we have developed `IPD`, a open-source `R` package that
implements methods under the umbrella of IPD.

## Installation

You can install the development version of `IPD` from
[GitHub](https://github.com/) with:

``` r
#-- Install devtools if it is not already installed

# install.packages("devtools")   

#-- Install the IPD package

devtools::install_github("awanafiaz/IPD")
```

## Usage

We provide a simple example to demonstrate the basic use of the
functions included in `IPD`.

### Example Setup

1.  Assume we have access to an algorithm $f(X)$ that can predict our
    outcome of interest $Y$.
2.  Consider we have two datasets: a labeled dataset,
    $(Y_L, X_L, f(X))$, and an unlabeled dataset, $(X_u, f(X_U))$. The
    labeled set is typically smaller in size compared to the unlabeled
    set.
3.  Our interest is in performing inference on $H_0: \beta_1^* = 0$ vs
    $H_1: \beta_1^* \ne 0$. Our inference model is:

$$Y_{val} = \beta_0^* + \beta_1^* X_{val} + \epsilon^*,$$

where $\epsilon^*$ is normally distributed with mean 0 and variance 1.
We do not observe $Y_U$, but only the predicted $\hat{Y}_U = f(X_U)$.

### Data Generation and Loading

``` r
#-- Load the IPD Library

library(IPD)

#-- Generate Example Data

set.seed(2023)

dat <- simdat(n = c(300, 300, 300), beta1 = 1)

head(dat)
#>            X1         X2           X3        X4         Y Yhat      set
#> 1  0.91621564  2.4201431  0.838470912 1.7995475  8.021947   NA training
#> 2  0.01705625 -0.3033343 -0.027365235 0.3132418  7.089503   NA training
#> 3 -0.87506732  1.0594201 -0.096441058 1.9303638  6.822057   NA training
#> 4  0.81385534 -0.3065133  0.702090044 3.4198123  7.792005   NA training
#> 5  0.36651430  1.3240447 -0.007524424 1.5385638  9.156789   NA training
#> 6  2.09079746  1.1153281  0.370369207 1.5954569 10.621224   NA training
```

### 1.1) PostPI Analytic Correction (Wang et al., 2020)

``` r
#-- Specify the Formula

formula <- Y - Yhat ~ X1

#-- Fit the PostPI Analytic Correction

IPD::ipd(formula, 
         
  method = "postpi_analytic", model = "ols", data = dat, label_index = "set")
#> $est
#> [1] 11.561038  0.955248
#> 
#> $se
#> [1] 0.2401493 0.1706713
#> 
#> $ci
#> [1] 11.0903542  0.6207384 12.0317220  1.2897577
#> 
#> attr(,"class")
#> [1] "ipd"
```

### 1.2) PostPI Bootstrap Correction (Wang et al., 2020)

``` r
#-- Fit the PostPI Bootstrap Correction

nboot <- 200

IPD::ipd(formula, 
         
  method = "postpi_boot", model = "ols", data = dat, label_index = "set", 
  
  nboot = nboot)
#> $est
#> [1] 12.70224908 -0.03039477
#> 
#> $se
#> [1] 0.2502446 0.1781699
#> 
#> $ci
#> [1] 12.2117787 -0.3796013 13.1927195  0.3188118
#> 
#> $nboot
#> [1] 200
#> 
#> attr(,"class")
#> [1] "ipd"
```

### 2. Prediction-Powered Inference (PPI; Angelopoulos et al., 2023)

``` r
#-- Fit the PPI Correction

IPD::ipd(formula, 
         
  method = "ppi", model = "ols", data = dat, label_index = "set")
#> $est
#> [1] 11.697627  0.817898
#> 
#> $se
#> [1] 0.2357630 0.1790149
#> 
#> $ci
#> [1] 11.2355396  0.4670352 12.1597136  1.1687608
#> 
#> attr(,"class")
#> [1] "ipd"
```

### 3. PPI++ (Angelopoulos et al., 2023)

``` r
#-- Fit the PPI++ Correction

IPD::ipd(formula, 
         
  method = "ppi_plusplus", model = "ols", data = dat, label_index = "set")
#> $est
#> X(Intercept)          XX1 
#>   11.7737782    0.8460529 
#> 
#> $se
#> [1] 0.2200044 0.1615084
#> 
#> $ci
#> X(Intercept)          XX1 X(Intercept)          XX1 
#>   11.3425774    0.5295022   12.2049789    1.1626035 
#> 
#> attr(,"class")
#> [1] "ipd"
```

### 4. Assumption-Lean and Data-Adaptive Post-Prediction Inference (POP-Inf; Miao et al., 2023)

``` r
#-- Fit the POP-Inf Correction

IPD::ipd(formula, 
         
  method = "popinf", model = "ols", data = dat, label_index = "set")
#> $est
#> [1] 11.7620945  0.8441666
#> 
#> $se
#> [1] 0.2307502 0.1721795
#> 
#> $ci
#> [1] 11.309833  0.506701 12.214357  1.181632
#> 
#> attr(,"class")
#> [1] "ipd"
```

## Vignette

For additional details, we provide more use cases and examples in the
package vignette:

``` r
vignette("IPD")
```

## Feedback

For questions, comments, or any other feedback, please contact the
developers.
