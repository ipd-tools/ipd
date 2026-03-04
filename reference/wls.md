# Weighted Least Squares

Computes the weighted least squares estimate of the coefficients.

## Usage

``` r
wls(X, Y, w = NULL, return_se = FALSE)
```

## Arguments

- X:

  (matrix): n x p matrix of covariates.

- Y:

  (vector): p-vector of outcome values.

- w:

  (vector, optional): n-vector of sample weights.

- return_se:

  (bool, optional): Whether to return the standard errors of the
  coefficients.

## Value

(list): A list containing the following:

- theta:

  (vector): p-vector of weighted least squares estimates of the
  coefficients.

- se:

  (vector): If return_se == TRUE, return the p-vector of standard errors
  of the coefficients.
