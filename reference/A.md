# Calculation of the matrix A based on single dataset

`A` function for the calculation of the matrix A based on single dataset

## Usage

``` r
A(
  X,
  Y,
  quant = NA,
  theta,
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

- theta:

  parameter theta

- method:

  indicates the method to be used for M-estimation. Options include
  "mean", "quantile", "ols", "logistic", and "poisson".

## Value

matrix A based on single dataset
