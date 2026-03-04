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
#> 601  0.6289534 -0.1691090 -1.3817136 -0.03280005 -1.8322586  0.2566973
#> 602  0.4179257 -1.9192325 -0.8440696 -0.51092478  2.9149733  0.6961904
#> 603  1.9767585 -1.5342664 -1.8365835  0.35643054  1.0645962  1.2211797
#> 604 -0.5062863 -1.1147612  1.0441708  0.41794614  3.3228702  1.3169226
#> 605 -1.1099689  1.5978116  0.2152101  0.57920526  0.6491361 -0.1833841
#> 606 -0.9487057 -0.6398051  0.8003496 -1.47515865  2.6688457  0.7158423
#>     set_label    .fitted    .resid
#> 601 unlabeled  1.3682125 -3.200471
#> 602 unlabeled  1.1509001  1.764073
#> 603 unlabeled  2.7561567 -1.691560
#> 604 unlabeled  0.1991641  3.123706
#> 605 unlabeled -0.4224968  1.071633
#> 606 unlabeled -0.2564310  2.925277
```
