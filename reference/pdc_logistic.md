# PDC Logistic

Helper function for PDC logistic regression estimation.

## Usage

``` r
pdc_logistic(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE)
```

## Arguments

- X_l:

  (matrix): n x p matrix of covariates in the labeled data.

- Y_l:

  (vector): n-vector of labeled outcomes.

- f_l:

  (vector): n-vector of predictions in the labeled data.

- X_u:

  (matrix): N x p matrix of covariates in the unlabeled data.

- f_u:

  (vector): N-vector of predictions in the unlabeled data.

- intercept:

  (Logical): Do the design matrices include intercept columns? Default
  is `TRUE`.

## Value

(list): A list containing the following:

- est:

  (vector): vector of PDC logistic regression coefficient estimates.

- se:

  (vector): vector of standard errors of the coefficients.

## Details

Prediction de-correlated inference: A safe approach for post-prediction
inference (Gan et al., 2024)
[doi:10.1111/anzs.12429](https://doi.org/10.1111/anzs.12429)

## Examples

``` r
dat <- simdat(model = "logistic")

form <- Y - f ~ X1

X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>
  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

pdc_logistic(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE)
#> $est
#> [1] 0.4910017 0.7564819
#> 
#> $se
#> (Intercept)          X1 
#>    0.128485    0.136995 
#> 
```
