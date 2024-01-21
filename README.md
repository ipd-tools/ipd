
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IPD

<!-- badges: start -->
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
#>          X1          X2         X3         X4         Y Yhat set
#> 1 0.1950123  0.71320956  1.7586172  1.7698342  8.914110   NA trn
#> 2 0.3582049 -0.31353315  0.1441178  0.6194818  5.987283   NA trn
#> 3 0.2590394  2.11566439 -0.1788052  2.2846895  5.786989   NA trn
#> 4 1.6731867  1.12790807  1.8621499  1.0567431  6.527794   NA trn
#> 5 0.2131707 -0.46024201  1.1124455 -0.5719875  8.871602   NA trn
#> 6 0.7334329 -0.07483258  1.0380555  4.3763716 12.422600   NA trn
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
