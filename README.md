
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
models using the provided `simdat` function by specifying the sizes of
the datasets, the effect size, residual variance, and the type of model.
The function currently supports “mean”, “quantile”, “ols”, “logistic”,
and “poisson” models.

``` r
#-- Load the IPD Library

library(IPD)

#-- Generate Example Data for Linear Regression

set.seed(12345)

dat <- simdat(c(10000, 500, 1000), effect = 1, sigma_Y = 1, model = "ols")

#-- Print First 6 Rows of Training, Labeled, and Unlabeled Subsets

options(digits=2)

head(dat[dat$set == "training",])
#>      X1    X2     X3    X4    Y  f      set
#> 1  0.59  0.75 -1.082  1.74  2.1 NA training
#> 2  0.71  0.84  0.041 -0.20  1.5 NA training
#> 3 -0.11 -1.84  0.491  1.12  2.1 NA training
#> 4 -0.45  0.20 -1.330  0.54 -2.4 NA training
#> 5  0.61  0.81 -0.391 -2.44  2.2 NA training
#> 6 -1.82 -0.39  1.168  0.29 -1.8 NA training

head(dat[dat$set == "labeled",])
#>          X1    X2    X3    X4      Y    f     set
#> 10001  0.16 -1.37  0.88  0.30  0.911 1.70 labeled
#> 10002  0.53 -0.45 -1.02  0.76 -0.769 0.28 labeled
#> 10003 -1.30 -0.81  1.22 -0.16 -0.720 0.57 labeled
#> 10004  1.25 -0.19 -1.43  1.27  1.578 0.61 labeled
#> 10005 -0.99  1.88  1.11 -0.26  0.931 0.85 labeled
#> 10006 -0.32 -0.90  0.62  0.58 -0.012 0.98 labeled

head(dat[dat$set == "unlabeled",])
#>          X1    X2     X3    X4      Y     f       set
#> 10501  1.51 -1.43  0.729  1.51  3.266  2.91 unlabeled
#> 10502 -0.91  1.59 -1.255  0.40 -0.606 -1.34 unlabeled
#> 10503  0.41 -1.64  0.415 -0.92  1.337  1.48 unlabeled
#> 10504 -0.88 -1.53  0.032 -0.54  0.827 -0.17 unlabeled
#> 10505  0.18  0.26  0.378 -0.87 -0.117  1.28 unlabeled
#> 10506  0.33 -0.57 -0.121  0.62 -0.033  0.93 unlabeled
```

The `simdat` function provides observed and unobserved outcomes for both
the labeled and unlabeled datasets, though in practice the observed
outcomes are not in the unlabeled set. We can visualize the
relationships between these variables:

<img src="man/figures/README-plot-1.png" width="100%" />

We can see that:

- The relationship between the true outcome and the covariate (plot A)
  is less variable than the relationship between the predicted outcome
  and the covariate (plot B)
- There is uncertainty in predicting the outcomes that needs to be
  accounted for (plot C)

### Model Fitting

We compare two non-`IPD` approaches to analyzin the data to methods
included in the `IPD` package.

#### 0.1) ‘Naive’ Regression Using the Predicted Outcomes

``` r
#--- Fit the Naive Regression

lm(f ~ X1, data = dat[dat$set == "unlabeled",]) |> 
  
  summary()
#> 
#> Call:
#> lm(formula = f ~ X1, data = dat[dat$set == "unlabeled", ])
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.9426 -0.6501 -0.0101  0.6543  3.0193 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   0.7418     0.0303    24.5   <2e-16 ***
#> X1            0.9805     0.0314    31.2   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.96 on 998 degrees of freedom
#> Multiple R-squared:  0.494,  Adjusted R-squared:  0.493 
#> F-statistic:  972 on 1 and 998 DF,  p-value: <2e-16
```

#### 0.2) ‘Classic’ Regression Using only the Labeled Data

``` r
#--- Fit the Classic Regression

lm(Y ~ X1, data = dat[dat$set == "labeled",]) |> 
  
  summary()
#> 
#> Call:
#> lm(formula = Y ~ X1, data = dat[dat$set == "labeled", ])
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -23.491  -0.885  -0.120   0.831  12.736 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   0.7658     0.0916    8.36  6.5e-16 ***
#> X1            1.0405     0.0936   11.12  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 2 on 498 degrees of freedom
#> Multiple R-squared:  0.199,  Adjusted R-squared:  0.197 
#> F-statistic:  124 on 1 and 498 DF,  p-value: <2e-16
```

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
  
  nboot = nboot) |> 
  
  summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: postpi_boot 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)   0.8215    0.0859   0.6532     0.99
#> X1            1.1019    0.0893   0.9269     1.28
```

#### 1.2) PostPI Analytic Correction (Wang et al., 2020)

``` r
#-- Fit the PostPI Analytic Correction

IPD::ipd(formula, 
         
  method = "postpi_analytic", model = "ols", data = dat, label = "set") |> 
  
  summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: postpi_analytic 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)   0.8146    0.0858   0.6464     0.98
#> X1            1.0993    0.0890   0.9249     1.27
```

#### 2. Prediction-Powered Inference (PPI; Angelopoulos et al., 2023)

``` r
#-- Fit the PPI Correction

IPD::ipd(formula, 
         
  method = "ppi", model = "ols", data = dat, label = "set") |> 
  
  summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: ppi 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)    0.808     0.084    0.644     0.97
#> X1             1.051     0.101    0.854     1.25
```

#### 3. PPI++ (Angelopoulos et al., 2023)

``` r
#-- Fit the PPI++ Correction

IPD::ipd(formula, 
         
  method = "ppi_plusplus", model = "ols", data = dat, label = "set") |> 
  
  summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: ppi_plusplus 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)   0.8006    0.0833   0.6373     0.96
#> X1            1.0492    0.1024   0.8485     1.25
```

#### 4. Assumption-Lean and Data-Adaptive Post-Prediction Inference (POP-Inf; Miao et al., 2023)

``` r
#-- Fit the POP-Inf Correction

IPD::ipd(formula, 
         
  method = "popinf", model = "ols", data = dat, label = "set") |> 
  
  summary()
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Method: popinf 
#> Model: ols 
#> Intercept: Yes 
#> 
#> Coefficients:
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)   0.7971    0.0831   0.6342     0.96
#> X1            1.0507    0.1002   0.8543     1.25
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
  
#-- Print the Model

print(fit_postpi)
#> 
#> Call:
#>  Y - f ~ X1 
#> 
#> Coefficients:
#> (Intercept)          X1 
#>        0.82        1.10

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
#>             Estimate Std.Error Lower.CI Upper.CI
#> (Intercept)   0.8215    0.0859   0.6532     0.99
#> X1            1.1019    0.0893   0.9269     1.28

#-- Tidy the Model Output

tidy(fit_postpi)
#>                    term estimate std.error conf.low conf.high
#> (Intercept) (Intercept)     0.82     0.086     0.65      0.99
#> X1                   X1     1.10     0.089     0.93      1.28

#-- Get a One-Row Summary of the Model

glance(fit_postpi)
#>        method model intercept  nobs       call
#> 1 postpi_boot   ols      TRUE 11500 Y - f ~ X1

#-- Augment the Original Data with Fitted Values and Residuals

df <- dat[which(dat$set != "training"),]

augmented_df <- augment(fit_postpi, data = df)

head(augmented_df)
#>          X1    X2    X3    X4      Y    f     set .fitted .resid
#> 10001  0.16 -1.37  0.88  0.30  0.911 1.70 labeled    1.00 -0.088
#> 10002  0.53 -0.45 -1.02  0.76 -0.769 0.28 labeled    1.41 -2.177
#> 10003 -1.30 -0.81  1.22 -0.16 -0.720 0.57 labeled   -0.61 -0.112
#> 10004  1.25 -0.19 -1.43  1.27  1.578 0.61 labeled    2.19 -0.616
#> 10005 -0.99  1.88  1.11 -0.26  0.931 0.85 labeled   -0.27  1.202
#> 10006 -0.32 -0.90  0.62  0.58 -0.012 0.98 labeled    0.47 -0.483
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
