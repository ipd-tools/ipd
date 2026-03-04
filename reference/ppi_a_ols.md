# PPI "All" OLS

Helper function for PPI "All" for OLS estimation

## Usage

``` r
ppi_a_ols(X_l, Y_l, f_l, X_u, f_u, w_l = NULL, w_u = NULL)
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

- w_l:

  (ndarray, optional): Sample weights for the labeled data set. Defaults
  to a vector of ones.

- w_u:

  (ndarray, optional): Sample weights for the unlabeled data set.
  Defaults to a vector of ones.

## Value

(list): A list containing the following:

- est:

  (vector): vector of PPI OLS regression coefficient estimates.

- se:

  (vector): vector of standard errors of the coefficients.

- rectifier_est:

  (vector): vector of the rectifier OLS regression coefficient
  estimates.

## Details

Another look at inference after prediction (Gronsbell et al., 2025)
<https://arxiv.org/pdf/2411.19908>

## Examples

``` r
dat <- simdat()

form <- Y - f ~ X1

X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>

  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>

  matrix(ncol = 1)

ppi_a_ols(X_l, Y_l, f_l, X_u, f_u)
#> $est
#>                  [,1]
#> (Intercept) 0.7524312
#> X1          1.2743293
#> 
#> $se
#> (Intercept)          X1 
#>  0.09985817  0.11738334 
#> 
#> $rectifier_est
#>                    [,1]
#> (Intercept)  0.16640753
#> X1          -0.05876525
#> 
```
