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
#> 1 (Intercept)    0.882    0.0931    0.700      1.06
#> 2 X1             0.859    0.0930    0.676      1.04


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
#>              X1         X2        X3          X4          Y         f set_label
#> 601  0.08481881 -0.2421913 0.1462140  0.34435404  0.3316460 0.9978946 unlabeled
#> 602  1.17282021 -0.6453282 0.2635753  2.23303739  2.5853066 1.9413377 unlabeled
#> 603 -0.41763512 -1.7609660 0.5193892  0.07973693  1.4059495 0.9457906 unlabeled
#> 604 -0.69857173 -0.8754849 0.2377680  0.83858637 -0.1150827 0.3497647 unlabeled
#> 605  0.55615268 -0.1253110 1.4646278 -0.70362920  0.1359178 2.9478723 unlabeled
#> 606  0.60065447 -0.7087863 0.2795710  0.59996179 -0.3521365 1.5726748 unlabeled
#>       .fitted     .resid
#> 601 0.8930284 -0.5613824
#> 602 1.9940361  0.5912705
#> 603 0.3845679  1.0213816
#> 604 0.1002729 -0.2153556
#> 605 1.3699967 -1.2340790
#> 606 1.4150305 -1.7671670
```
