# One-step update for obtaining estimator

`optim_est` function for One-step update for obtaining estimator

## Usage

``` r
optim_est(
  X_l,
  X_u,
  Y_l,
  f_l,
  f_u,
  w,
  theta,
  quant = NA,
  method = c("ols", "quantile", "mean", "logistic", "poisson")
)
```

## Arguments

- X_l:

  Array or data.frame containing observed covariates in labeled data.

- X_u:

  Array or data.frame containing observed or predicted covariates in
  unlabeled data.

- Y_l:

  Array or data.frame of observed outcomes in labeled data.

- f_l:

  Array or data.frame of predicted outcomes in labeled data.

- f_u:

  Array or data.frame of predicted outcomes in unlabeled data.

- w:

  weights vector PSPA linear regression (d-dimensional, where d equals
  the number of covariates).

- theta:

  parameter theta

- quant:

  quantile for quantile estimation

- method:

  indicates the method to be used for M-estimation. Options include
  "mean", "quantile", "ols", "logistic", and "poisson".

## Value

estimator
