
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IPD

<!-- badges: start -->

[![R-CMD-check](https://github.com/awanafiaz/IPD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/awanafiaz/IPD/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of IPD is to …

## Installation

You can install the development version of IPD like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(IPD)
## basic example code
head(simdat())
#>            X1         X2         X3       X4        Y Yhat set
#> 1  0.04050977  1.6491773  0.4520383 3.923058 18.13042   NA trn
#> 2  1.58266222  0.4323382  1.6251682 2.301778 15.92563   NA trn
#> 3  0.11305996 -0.2203519  0.1929704 3.428823 12.24805   NA trn
#> 4  0.83598899  0.0802341  1.8760026 1.422858 14.24691   NA trn
#> 5  1.33509244  2.3725053 -1.1384842 2.868554 13.35183   NA trn
#> 6 -0.08802149 -0.4022710  0.6951442 2.441598 12.17421   NA trn
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
