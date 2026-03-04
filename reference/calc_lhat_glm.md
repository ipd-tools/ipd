# Estimate PPI++ Power Tuning Parameter

Calculates the optimal value of lhat for the prediction-powered
confidence interval for GLMs.

## Usage

``` r
calc_lhat_glm(
  grads,
  grads_hat,
  grads_hat_unlabeled,
  inv_hessian,
  coord = NULL,
  clip = FALSE
)
```

## Arguments

- grads:

  (matrix): n x p matrix gradient of the loss function with respect to
  the parameter evaluated at the labeled data.

- grads_hat:

  (matrix): n x p matrix gradient of the loss function with respect to
  the model parameter evaluated using predictions on the labeled data.

- grads_hat_unlabeled:

  (matrix): N x p matrix gradient of the loss function with respect to
  the parameter evaluated using predictions on the unlabeled data.

- inv_hessian:

  (matrix): p x p matrix inverse of the Hessian of the loss function
  with respect to the parameter.

- coord:

  (int, optional): Coordinate for which to optimize `lhat`. If `None`,
  it optimizes the total variance over all coordinates. Must be in (1,
  ..., d) where d is the shape of the estimand.

- clip:

  (boolean, optional): Whether to clip the value of lhat to be
  non-negative. Defaults to `False`.

## Value

(float): Optimal value of `lhat` in \[0,1\].
