# PostPI OLS (Bootstrap Correction)

Helper function for PostPI OLS estimation (bootstrap correction)

## Usage

``` r
postpi_boot_ols(
  X_l,
  Y_l,
  f_l,
  X_u,
  f_u,
  nboot = 100,
  se_type = "par",
  rel_func = "lm",
  scale_se = TRUE,
  n_t = Inf
)
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

- nboot:

  (integer): Number of bootstrap samples. Defaults to 100.

- se_type:

  (string): Which method to calculate the standard errors. Options
  include "par" (parametric) or "npar" (nonparametric). Defaults to
  "par".

- rel_func:

  (string): Method for fitting the relationship model. Options include
  "lm" (linear model), "rf" (random forest), and "gam" (generalized
  additive model). Defaults to "lm".

- scale_se:

  (boolean): Logical argument to scale relationship model error
  variance. Defaults to TRUE; FALSE option is retained for posterity.

- n_t:

  (integer, optional) Size of the dataset used to train the prediction
  function (necessary if `n_t` \< `nrow(X_l)`. Defaults to `Inf`.

## Value

A list of outputs: estimate of inference model parameters and
corresponding standard error based on both parametric and non-parametric
bootstrap methods.

## Details

Methods for correcting inference based on outcomes predicted by machine
learning (Wang et al., 2020)
[doi:10.1073/pnas.2001238117](https://doi.org/10.1073/pnas.2001238117)

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

postpi_boot_ols(X_l, Y_l, f_l, X_u, f_u, nboot = 200)
#> $est
#> [1] 0.8403114 1.0299061
#> 
#> $se
#> [1] 0.1150041 0.1175506
#> 
```
