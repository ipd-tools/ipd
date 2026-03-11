# Tidy an ipd fit

Tidy an ipd fit

## Usage

``` r
# S3 method for class 'ipd'
tidy(x, ...)
```

## Arguments

- x:

  An object of class `ipd`.

- ...:

  Ignored.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
columns `term, estimate, std.error, conf.low, conf.high`.

## Examples

``` r
dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

tidy(fit)
#> # A tibble: 2 × 5
#>   term        estimate std.error conf.low conf.high
#>   <chr>          <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.740     0.107    0.529     0.950
#> 2 X1             1.13      0.130    0.880     1.39 
```
