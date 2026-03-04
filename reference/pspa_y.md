# PSPA M-Estimation for ML-predicted labels

`pspa_y` function conducts post-prediction M-Estimation.

## Usage

``` r
pspa_y(
  X_l = NA,
  X_u = NA,
  Y_l,
  f_l,
  f_u,
  alpha = 0.05,
  weights = NA,
  quant = NA,
  intercept = FALSE,
  method
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

- alpha:

  Specifies the confidence level as 1 - alpha for confidence intervals.

- weights:

  weights vector PSPA linear regression (d-dimensional, where d equals
  the number of covariates).

- quant:

  quantile for quantile estimation

- intercept:

  Boolean indicating if the input covariates' data contains the
  intercept (TRUE if the input data contains)

- method:

  indicates the method to be used for M-estimation. Options include
  "mean", "quantile", "ols", "logistic", and "poisson".

## Value

A summary table presenting point estimates, standard error, confidence
intervals (1 - alpha), P-values, and weights.
