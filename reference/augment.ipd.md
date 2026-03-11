# Augment data from an ipd fit

Augment data from an ipd fit

## Usage

``` r
# S3 method for class 'ipd'
augment(x, data = x@data_u, ...)
```

## Arguments

- x:

  An object of class `ipd`.

- data:

  A `data.frame` to augment; defaults to `x@data_u`.

- ...:

  Ignored.

## Value

The `data.frame` with columns `.fitted` and `.resid`.

## Examples

``` r
dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

augmented_df <- augment(fit)

head(augmented_df)
#>             X1         X2         X3          X4          Y          f
#> 601  0.2360958 -0.5785249  1.1436235  0.04909540 -1.1141895  2.0340122
#> 602  0.6289534 -0.1691090 -1.3817136 -0.03280005 -1.8322586  0.2564554
#> 603  0.4179257 -1.9192325 -0.8440696 -0.51092478  2.9149733  0.7119315
#> 604  1.9767585 -1.5342664 -1.8365835  0.35643054  1.0645962  1.2354387
#> 605 -0.5062863 -1.1147612  1.0441708  0.41794614  3.3228702  1.3354883
#> 606 -1.1099689  1.5978116  0.2152101  0.57920526  0.6491361 -0.1933974
#>     set_label    .fitted    .resid
#> 601 unlabeled  0.9681065 -2.082296
#> 602 unlabeled  1.3734770 -3.205736
#> 603 unlabeled  1.1557278  1.759246
#> 604 unlabeled  2.7642111 -1.699615
#> 605 unlabeled  0.2020786  3.120792
#> 606 unlabeled -0.4208319  1.069968
```
