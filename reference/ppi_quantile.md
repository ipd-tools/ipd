# PPI Quantile Estimation

Helper function for PPI quantile estimation

## Usage

``` r
ppi_quantile(Y_l, f_l, f_u, q, alpha = 0.05, exact_grid = FALSE)
```

## Arguments

- Y_l:

  (vector): n-vector of labeled outcomes.

- f_l:

  (vector): n-vector of predictions in the labeled data.

- f_u:

  (vector): N-vector of predictions in the unlabeled data.

- q:

  (float): Quantile to estimate. Must be in the range (0, 1).

- alpha:

  (scalar): type I error rate for hypothesis testing - values in (0, 1);
  defaults to 0.05.

- exact_grid:

  (bool, optional): Whether to compute the exact solution (TRUE) or an
  approximate solution based on a linearly spaced grid of 5000 values
  (FALSE).

## Value

tuple: Lower and upper bounds of the prediction-powered confidence
interval for the quantile.

## Details

Prediction Powered Inference (Angelopoulos et al., 2023)
[doi:10.1126/science.adi6000](https://doi.org/10.1126/science.adi6000)

## Examples

``` r
dat <- simdat(model = "quantile")

form <- Y - f ~ X1

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>

  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

ppi_quantile(Y_l, f_l, f_u, q = 0.5)
#> [1] 0.8377598 1.4303978
```
