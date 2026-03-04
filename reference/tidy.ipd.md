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
#> 1 (Intercept)    0.773    0.0843    0.608     0.938
#> 2 X1             0.947    0.0874    0.776     1.12 
```
