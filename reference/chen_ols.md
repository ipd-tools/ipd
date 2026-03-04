# Chen & Chen OLS

Helper function for Chen & Chen OLS estimation

## Usage

``` r
chen_ols(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE)
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

- intercept:

  (Logical): Do the design matrices include intercept columns? Default
  is `TRUE`.

## Value

(list): A list containing the following:

- est:

  (vector): vector of Chen & Chen OLS regression coefficient estimates.

- se:

  (vector): vector of standard errors of the coefficients.

## Details

Another look at inference after prediction (Gronsbell et al., 2025)
<https://arxiv.org/pdf/2411.19908>

## Examples

``` r
dat <- simdat(model = "ols")

form <- Y - f ~ X1

X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>
  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

chen_ols(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE)
#> $est
#> [1] 0.6637123 1.0015489
#> 
#> $se
#> (Intercept)          X1 
#>  0.09240834  0.09422155 
#> 
```
