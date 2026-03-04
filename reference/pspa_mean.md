# PSPA Mean Estimation

Helper function for PSPA mean estimation

## Usage

``` r
pspa_mean(Y_l, f_l, f_u, weights = NA, alpha = 0.05)
```

## Arguments

- Y_l:

  (vector): n-vector of labeled outcomes.

- f_l:

  (vector): n-vector of predictions in the labeled data.

- f_u:

  (vector): N-vector of predictions in the unlabeled data.

- weights:

  (array): 1-dimensional array of weights vector for variance reduction.
  PSPA will estimate the weights if not specified.

- alpha:

  (scalar): type I error rate for hypothesis testing - values in (0, 1);
  defaults to 0.05.

## Value

A list of outputs: estimate of inference model parameters and
corresponding standard error.

## Details

Post-prediction adaptive inference (Miao et al., 2023)
<https://arxiv.org/abs/2311.14220>

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

pspa_mean(Y_l = Y_l, f_l = f_l, f_u = f_u)
#> $est
#> [1] 1.004022
#> 
#> $se
#> [1] 0.06009168
#> 
```
