# OLS Gradient and Hessian

Computes the statistics needed for the OLS-based prediction-powered
inference.

## Usage

``` r
ols_get_stats(
  est,
  X_l,
  Y_l,
  f_l,
  X_u,
  f_u,
  w_l = NULL,
  w_u = NULL,
  use_u = TRUE
)
```

## Arguments

- est:

  (vector): Point estimates of the coefficients.

- X_l:

  (matrix): Covariates for the labeled data set.

- Y_l:

  (vector): Labels for the labeled data set.

- f_l:

  (vector): Predictions for the labeled data set.

- X_u:

  (matrix): Covariates for the unlabeled data set.

- f_u:

  (vector): Predictions for the unlabeled data set.

- w_l:

  (vector, optional): Sample weights for the labeled data set.

- w_u:

  (vector, optional): Sample weights for the unlabeled data set.

- use_u:

  (boolean, optional): Whether to use the unlabeled data set.

## Value

(list): A list containing the following:

- grads:

  (matrix): n x p matrix gradient of the loss function with respect to
  the coefficients.

- grads_hat:

  (matrix): n x p matrix gradient of the loss function with respect to
  the coefficients, evaluated using the labeled predictions.

- grads_hat_unlabeled:

  (matrix): N x p matrix gradient of the loss function with respect to
  the coefficients, evaluated using the unlabeled predictions.

- inv_hessian:

  (matrix): p x p matrix inverse Hessian of the loss function with
  respect to the coefficients.
