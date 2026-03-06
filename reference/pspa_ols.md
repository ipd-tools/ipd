# PSPA OLS Estimation

Helper function for PSPA OLS for linear regression

## Usage

``` r
pspa_ols(X_l, Y_l, f_l, X_u, f_u, weights = NA, alpha = 0.05)
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

- weights:

  (array): p-dimensional array of weights vector for variance reduction.
  PSPA will estimate the weights if not specified.

- alpha:

  (scalar): type I error rate for hypothesis testing - values in (0, 1);
  defaults to 0.05.

## Value

A list of outputs: estimate of inference model parameters and
corresponding standard error.

## Details

Post-prediction adaptive inference (Miao et al., 2023)
[doi:10.48550/arXiv.2311.14220](https://doi.org/10.48550/arXiv.2311.14220)

## Examples

``` r
dat <- simdat(model = "ols")

form <- Y - f ~ X1

X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>
  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

pspa_ols(X_l, Y_l, f_l, X_u, f_u)
#> $est
#> [1] 0.7294793 0.9438626
#> 
#> $se
#> [1] 0.0976508 0.0885703
#> 
```
