# PPI++ Logistic Regression (Point Estimate)

Helper function for PPI++ logistic regression (point estimate)

## Usage

``` r
ppi_plusplus_logistic_est(
  X_l,
  Y_l,
  f_l,
  X_u,
  f_u,
  lhat = NULL,
  coord = NULL,
  opts = NULL,
  w_l = NULL,
  w_u = NULL
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

- lhat:

  (float, optional): Power-tuning parameter (see
  [doi:10.48550/arXiv.2311.01453](https://doi.org/10.48550/arXiv.2311.01453)
  ). The default value, `NULL`, will estimate the optimal value from the
  data. Setting `lhat = 1` recovers PPI with no power tuning, and
  setting `lhat = 0` recovers the classical point estimate.

- coord:

  (int, optional): Coordinate for which to optimize `lhat = 1`. If
  `NULL`, it optimizes the total variance over all coordinates. Must be
  in (1, ..., d) where d is the dimension of the estimand.

- opts:

  (list, optional): Options to pass to the optimizer. See ?optim for
  details.

- w_l:

  (ndarray, optional): Sample weights for the labeled data set. Defaults
  to a vector of ones.

- w_u:

  (ndarray, optional): Sample weights for the unlabeled data set.
  Defaults to a vector of ones.

## Value

(vector): vector of prediction-powered point estimates of the logistic
regression coefficients.

## Details

PPI++: Efficient Prediction Powered Inference (Angelopoulos et al.,
2023)
[doi:10.48550/arXiv.2311.01453](https://doi.org/10.48550/arXiv.2311.01453)

## Examples

``` r
dat <- simdat(model = "logistic")

form <- Y - f ~ X1

X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>

  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u)
#>           [,1]
#> [1,] 0.4458863
#> [2,] 0.6481211
```
