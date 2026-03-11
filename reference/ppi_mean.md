# PPI Mean Estimation

Helper function for PPI mean estimation

## Usage

``` r
ppi_mean(Y_l, f_l, f_u, alpha = 0.05, alternative = "two-sided")
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

## Value

tuple: Lower and upper bounds of the prediction-powered confidence
interval for the mean.

## Details

Prediction Powered Inference (Angelopoulos et al., 2023)
[doi:10.1126/science.adi6000](https://doi.org/10.1126/science.adi6000)

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

ppi_mean(Y_l, f_l, f_u)
#>          lower    upper
#> [1,] 0.7840516 1.177333
```
