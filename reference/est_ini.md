# Initial estimation

`est_ini` function for initial estimation

## Usage

``` r
est_ini(
  X,
  Y,
  quant = NA,
  method = c("ols", "quantile", "mean", "logistic", "poisson")
)
```

## Arguments

- X:

  Array or data.frame containing covariates

- Y:

  Array or data.frame of outcomes

- quant:

  quantile for quantile estimation

- method:

  indicates the method to be used for M-estimation. Options include
  "mean", "quantile", "ols", "logistic", and "poisson".

## Value

initial estimator
