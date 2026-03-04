# Estimating equation

`psi` function for estimating equation

## Usage

``` r
psi(
  X,
  Y,
  theta,
  quant = NA,
  method = c("ols", "quantile", "mean", "logistic", "poisson")
)
```

## Arguments

- X:

  Array or data.frame containing covariates

- Y:

  Array or data.frame of outcomes

- theta:

  parameter theta

- quant:

  quantile for quantile estimation

- method:

  indicates the method to be used for M-estimation. Options include
  "mean", "quantile", "ols", "logistic", and "poisson".

## Value

estimating equation
