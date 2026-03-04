# Warn on differing factor levels between labeled and unlabeled data

Warn on differing factor levels between labeled and unlabeled data

## Usage

``` r
.warn_differing_levels(data_l, data_u, factor_vars)
```

## Arguments

- data_l:

  labeled data.frame

- data_u:

  unlabeled data.frame

- factor_vars:

  character vector of factor column names

## Value

Invisibly returns `NULL`. Messages are printed for each variable whose
levels differ.
