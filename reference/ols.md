# Ordinary Least Squares

Computes the ordinary least squares coefficients.

## Usage

``` r
ols(X, Y, return_se = FALSE)
```

## Arguments

- X:

  (matrix): n x p matrix of covariates.

- Y:

  (vector): p-vector of outcome values.

- return_se:

  (bool, optional): Whether to return the standard errors of the
  coefficients.

## Value

(list): A list containing the following:

- theta:

  (vector): p-vector of ordinary least squares estimates of the
  coefficients.

- se:

  (vector): If return_se == TRUE, return the p-vector of standard errors
  of the coefficients.
