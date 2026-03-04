# Build design matrices and outcome vectors

Build design matrices and outcome vectors

## Usage

``` r
.build_design(formula, data_l, data_u, intercept = TRUE, na_action = "na.fail")
```

## Arguments

- formula:

  two-sided formula `Y - f ~ X...`

- data_l:

  labeled data.frame

- data_u:

  unlabeled data.frame

- intercept:

  include intercept?

- na_action:

  "na.fail" or "na.omit"

## Value

list(X_l, Y_l, f_l, X_u, f_u)
