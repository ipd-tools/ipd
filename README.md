
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
estimates and anti-conservative inference ([Hoffman et al.,
2023](https://arxiv.org/abs/2401.08702)). The statistical challenges
encountered when drawing inference on predicted data (IPD) include:

1.  Understanding the relationship between predicted outcomes and their
    true, unobserved counterparts
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
implements these methods under the umbrella of IPD.

This README provides an overview of the package, including installation
instructions, basic usage examples, and links to further documentation.
The examples show how to generate data, fit models, and use custom
methods provided by the package.

## Installation

To install the development version of `IPD` from
[GitHub](https://github.com/awanafiaz/IPD), you can use the `devtools`
package:

``` r
#-- Install devtools if it is not already installed

install.packages("devtools")   

#-- Install the IPD package from GitHub

devtools::install_github("awanafiaz/IPD")
```

## Usage

We provide a simple example to demonstrate the basic use of the
functions included in the `IPD` package.

### Example Setup

1.  We have two datasets: a labeled dataset,
    $\mathcal{L} = \left\{Y^\mathcal{L}, X^\mathcal{L}, f\left(X^\mathcal{L}\right)\right\}$,
    and an unlabeled dataset,
    $\left\{X^\mathcal{U}, f\left(X^\mathcal{U}\right)\right\}$. The
    labeled set is typically smaller in size compared to the unlabeled
    set.
2.  We have access to an algorithm $f(X)$ that can predict our outcome
    of interest $Y$.
3.  Our interest is in performing inference on a quantity such as the
    outcome mean or quantile, or to recover a downstream inferential
    (mean) model:

$$\mathbb{E}\left[Y^{\mathcal{U}} \mid \boldsymbol{X}^{\mathcal{U}}\right] = g^{-1}\left(\boldsymbol{X}^{\mathcal{U}'}\beta\right),$$

where $\beta$ is a vector of regression coefficients and $g(\cdot)$ is a
given link function, such as the identity link for linear regression,
the logistic link for logistic regression, or the log link for Poisson
regression. However, we do not observe $Y^\mathcal{U}$, only the
predicted $f(X^\mathcal{U})$. We can use methods for IPD to obtain
corrected estimates and standard errors when we replace these unobserved
$Y^\mathcal{U}$ by $f(X^\mathcal{U})$.

### Data Generation

You can generate synthetic datasets for different types of regression
models using the provided `simdat` function. The function currently
supports “mean”, “quantile”, “ols”, “logistic”, and “poisson” models.

``` r
#-- Load the IPD Library

library(IPD)

#-- Generate Example Data for Linear Regression

set.seed(12345)

dat <- simdat(c(10000, 500, 1000), effect = 1, sigma_Y = 1, model = "ols")

head(dat[dat$set == "training",])
#>           X1         X2          X3         X4         Y  f      set
#> 1  0.5855288  0.7487320 -1.08153963  1.7380841 10.084643 NA training
#> 2  0.7094660  0.8386136  0.04123505 -0.1996555  1.679300 NA training
#> 3 -0.1093033 -1.8377414  0.49054926  1.1245687  7.163712 NA training
#> 4 -0.4534972  0.2044950 -1.33044649  0.5389589 -7.564027 NA training
#> 5  0.6058875  0.8086635 -0.39129058 -2.4365934 24.318713 NA training
#> 6 -1.8179560 -0.3886214  1.16798808  0.2939337  2.763427 NA training

head(dat[dat$set == "labeled",])
#>               X1         X2         X3         X4          Y         f     set
#> 10001  0.1612775 -1.3730215  0.8816500  0.3013700  3.0791726 12.421518 labeled
#> 10002  0.5319214 -0.4466659 -1.0221340  0.7586302 -1.4589775 -3.676138 labeled
#> 10003 -1.2979034 -0.8122093  1.2188757 -0.1630360  4.2082109 13.639386 labeled
#> 10004  1.2455783 -0.1915237 -1.4300525  1.2729837 -0.1437806 -6.324083 labeled
#> 10005 -0.9916905  1.8768000  1.1086366 -0.2568603  4.8118202 13.043107 labeled
#> 10006 -0.3181931 -0.9006203  0.6209725  0.5786088  1.8822368  9.673776 labeled

head(dat[dat$set == "unlabeled",])
#>               X1         X2          X3         X4         Y         f
#> 10501  1.5058747 -1.4321016  0.72924860  1.5121492 12.875098 12.856373
#> 10502 -0.9105015  1.5949058 -1.25467772  0.3954227 -5.286333 -7.370062
#> 10503  0.4059670 -1.6423932  0.41524420 -0.9161881  4.675317  8.358913
#> 10504 -0.8801105 -1.5265113  0.03174338 -0.5398788  1.919659  3.649073
#> 10505  0.1844883  0.2574885  0.37794643 -0.8719697  2.878582  7.824472
#> 10506  0.3271606 -0.5748017 -0.12104284  0.6207108  1.407498  3.931931
#>             set
#> 10501 unlabeled
#> 10502 unlabeled
#> 10503 unlabeled
#> 10504 unlabeled
#> 10505 unlabeled
#> 10506 unlabeled
```

``` r

library(tidyverse)
#> Warning: package 'ggplot2' was built under R version 4.2.3
#> Warning: package 'tibble' was built under R version 4.2.3
#> Warning: package 'tidyr' was built under R version 4.2.3
#> Warning: package 'readr' was built under R version 4.2.3
#> Warning: package 'purrr' was built under R version 4.2.3
#> Warning: package 'dplyr' was built under R version 4.2.3
#> Warning: package 'stringr' was built under R version 4.2.3
#> Warning: package 'lubridate' was built under R version 4.2.3
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

library(patchwork)
#> Warning: package 'patchwork' was built under R version 4.2.3

dat_labeled <- dat[dat$set == "labeled",]

p1 <- dat_labeled |> ggplot(aes(x = X1, y = Y)) + geom_point()

p2 <- dat_labeled |> ggplot(aes(x = X1, y = f)) + geom_point()

p3 <- dat_labeled |> ggplot(aes(x = f, y = Y)) + geom_point()

p1 + p2 + p3
```

<img src="man/figures/README-plot-1.png" width="100%" />

### Model Fitting

You can fit the various IPD methods to your data and obtain summaries
using the provided wrapper function, `ipd()`:

#### 1.1) PostPI Bootstrap Correction (Wang et al., 2020)

``` r
#-- Specify the Formula

formula <- Y - f ~ X1

#-- Fit the PostPI Bootstrap Correction

nboot <- 200

IPD::ipd(formula, 
         
  method = "postpi_boot", model = "ols", data = dat, label = "set", 
  
  nboot = nboot)
#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows
#> (Intercept) X1
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#> 4.853974127 0.002707768
```

#### 1.2) PostPI Analytic Correction (Wang et al., 2020)

``` r
#-- Fit the PostPI Analytic Correction

IPD::ipd(formula, 
         
  method = "postpi_analytic", model = "ols", data = dat, label = "set")
#> (Intercept) X1
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#>   5.2443801   0.9654006
```

#### 2. Prediction-Powered Inference (PPI; Angelopoulos et al., 2023)

``` r
#-- Fit the PPI Correction

IPD::ipd(formula, 
         
  method = "ppi", model = "ols", data = dat, label = "set")
#> (Intercept) X1
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#>    5.238435    1.870321
```

#### 3. PPI++ (Angelopoulos et al., 2023)

``` r
#-- Fit the PPI++ Correction

IPD::ipd(formula, 
         
  method = "ppi_plusplus", model = "ols", data = dat, label = "set")
#> (Intercept) X1
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#>    5.130324    1.841966
```

#### 4. Assumption-Lean and Data-Adaptive Post-Prediction Inference (POP-Inf; Miao et al., 2023)

``` r
#-- Fit the POP-Inf Correction

IPD::ipd(formula, 
         
  method = "popinf", model = "ols", data = dat, label = "set")
#> (Intercept) X1
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#>    5.110967    1.851909
```

### Printing and Tidying

The package also provides custom print, summary, tidy, glance, and
augment methods to facilitate easy model inspection:

``` r

#-- Fit the PostPI Bootstrap Correction

nboot <- 200

fit_postpi <- IPD::ipd(formula, 
         
  method = "postpi_boot", model = "ols", data = dat, label = "set", 
  
  nboot = nboot)
#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows

#> Warning: 'newdata' had 1000 rows but variables found have 500 rows
#> (Intercept) X1
  
#-- Print the Model

print(fit_postpi)
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#> 4.853974127 0.002707768

#-- Summarize the Model

summ_fit_postpi <- summary(fit_postpi)
  
#-- Print the Model Summary

print(summ_fit_postpi)
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: postpi_boot 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>               Estimate  Std.Error   Lower.CI Upper.CI
#> (Intercept)  4.8539741  0.5965569  3.6847441   6.0232
#> X1           0.0027078  0.6216147 -1.2156347   1.2211

#-- Tidy the Model Output

tidy(fit_postpi)
#>                    term    estimate std.error  conf.low conf.high
#> (Intercept) (Intercept) 4.853974127 0.5965569  3.684744  6.023204
#> X1                   X1 0.002707768 0.6216147 -1.215635  1.221050

#-- Get a One-Row Summary of the Model

glance(fit_postpi)
#>        method model intercept  nobs       call
#> 1 postpi_boot   ols      TRUE 11500 Y - f ~ X1

#-- Augment the Original Data with Fitted Values and Residuals

df <- dat[which(dat$set != "training"),]

augmented_df <- augment(fit_postpi, data = df)

head(augmented_df)
#>               X1         X2         X3         X4          Y         f     set
#> 10001  0.1612775 -1.3730215  0.8816500  0.3013700  3.0791726 12.421518 labeled
#> 10002  0.5319214 -0.4466659 -1.0221340  0.7586302 -1.4589775 -3.676138 labeled
#> 10003 -1.2979034 -0.8122093  1.2188757 -0.1630360  4.2082109 13.639386 labeled
#> 10004  1.2455783 -0.1915237 -1.4300525  1.2729837 -0.1437806 -6.324083 labeled
#> 10005 -0.9916905  1.8768000  1.1086366 -0.2568603  4.8118202 13.043107 labeled
#> 10006 -0.3181931 -0.9006203  0.6209725  0.5786088  1.8822368  9.673776 labeled
#>        .fitted      .resid
#> 10001 4.854411 -1.77523821
#> 10002 4.855414 -6.31439198
#> 10003 4.850460 -0.64224885
#> 10004 4.857347 -5.00112744
#> 10005 4.851289 -0.03946866
#> 10006 4.853113 -2.97087578
```

## Vignette

For additional details, we provide more use cases and examples in the
package vignette:

``` r
vignette("ipd")
```

## Feedback

For questions, comments, or any other feedback, please contact the
developers at [ssalerno@fredhutch.org](ssalerno@fredhutch.org).

## Contributing

Contributions are welcome! Please open an issue or submit a pull request
on GitHub.

## License

This package is licensed under the MIT License.
