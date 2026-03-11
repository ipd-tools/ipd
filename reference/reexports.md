# tidy re-exported from generics packages

These objects are imported from other packages. Follow the links below
to see their documentation.

- generics:

  [`augment`](https://generics.r-lib.org/reference/augment.html),
  [`glance`](https://generics.r-lib.org/reference/glance.html),
  [`tidy`](https://generics.r-lib.org/reference/tidy.html)

## Value

A wrapper for the `tidy` generic. See
[`tidy`](https://generics.r-lib.org/reference/tidy.html) for details.

A wrapper for the `glance` generic. See
[`glance`](https://generics.r-lib.org/reference/glance.html) for
details.

A wrapper for the `augment` generic. See
[`augment`](https://generics.r-lib.org/reference/augment.html) for
details.

## See also

[`tidy`](https://generics.r-lib.org/reference/tidy.html)

[`glance`](https://generics.r-lib.org/reference/glance.html)

[`augment`](https://generics.r-lib.org/reference/augment.html)

## Examples

``` r
dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

tidy(fit)
#> # A tibble: 2 × 5
#>   term        estimate std.error conf.low conf.high
#>   <chr>          <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.705    0.0804    0.547     0.862
#> 2 X1             0.979    0.0839    0.814     1.14 


dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

glance(fit)
#> # A tibble: 1 × 6
#>   method model intercept nobs_labeled nobs_unlabeled call      
#>   <chr>  <chr> <lgl>            <int>          <int> <chr>     
#> 1 pspa   ols   TRUE               300            300 Y - f ~ X1


dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

augmented_df <- augment(fit)

head(augmented_df)
#>             X1          X2         X3          X4            Y         f
#> 601  0.4620383  1.21499611  0.5345855 -0.40777087  0.329666273 1.6930041
#> 602  0.3899389 -0.64306635 -0.4182949 -0.79746344 -0.148048554 0.5528934
#> 603 -0.7429560  0.31666331  0.5612329  0.61575629 -1.894593153 0.5041894
#> 604  0.6816739  0.07468319 -1.0140786 -0.22577715  0.007768632 0.2278700
#> 605  0.6658204 -0.92026387 -0.1138219 -1.13904740  2.813899904 1.1473102
#> 606  1.3897301 -0.40171531  0.6574756  0.09777209 -0.114173350 2.7292120
#>     set_label     .fitted     .resid
#> 601 unlabeled  1.12809139 -0.7984251
#> 602 unlabeled  1.05637443 -1.2044230
#> 603 unlabeled -0.07051029 -1.8240829
#> 604 unlabeled  1.34656179 -1.3387932
#> 605 unlabeled  1.33079238  1.4831075
#> 606 unlabeled  2.05086155 -2.1650349
```
