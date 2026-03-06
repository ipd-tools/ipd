# PPI++ Mean Estimation

Helper function for PPI++ mean estimation

## Usage

``` r
ppi_plusplus_mean(
  Y_l,
  f_l,
  f_u,
  alpha = 0.05,
  alternative = "two-sided",
  lhat = NULL,
  coord = NULL,
  w_l = NULL,
  w_u = NULL
)
```

## Arguments

- Y_l:

  (vector): n-vector of labeled outcomes.

- f_l:

  (vector): n-vector of predictions in the labeled data.

- f_u:

  (vector): N-vector of predictions in the unlabeled data.

- alpha:

  (scalar): type I error rate for hypothesis testing - values in (0, 1);
  defaults to 0.05.

- alternative:

  (string): Alternative hypothesis. Must be one of `"two-sided"`,
  `"less"`, or `"greater"`.

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

- w_l:

  (ndarray, optional): Sample weights for the labeled data set. Defaults
  to a vector of ones.

- w_u:

  (ndarray, optional): Sample weights for the unlabeled data set.
  Defaults to a vector of ones.

## Value

tuple: Lower and upper bounds of the prediction-powered confidence
interval for the mean.

## Details

PPI++: Efficient Prediction Powered Inference (Angelopoulos et al.,
2023)
[doi:10.48550/arXiv.2311.01453](https://doi.org/10.48550/arXiv.2311.01453)

## Examples

``` r
dat <- simdat(model = "mean")

form <- Y - f ~ 1

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>

  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

ppi_plusplus_mean(Y_l, f_l, f_u)
#>          lower    upper
#> [1,] 0.9910999 1.202698
```
