# Glance at an ipd fit

Glance at an ipd fit

## Usage

``` r
# S3 method for class 'ipd'
glance(x, ...)
```

## Arguments

- x:

  An object of class `ipd`.

- ...:

  Ignored.

## Value

A one-row [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
summarizing the fit.

## Examples

``` r
dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

glance(fit)
#> # A tibble: 1 × 6
#>   method model intercept nobs_labeled nobs_unlabeled call      
#>   <chr>  <chr> <lgl>            <int>          <int> <chr>     
#> 1 pspa   ols   TRUE               300            300 Y - f ~ X1
```
