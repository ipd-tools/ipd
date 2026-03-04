# Validate and split input data

Validate and split input data

## Usage

``` r
.parse_inputs(data, label = NULL, unlabeled_data = NULL, na_action = "na.fail")
```

## Arguments

- data:

  data.frame or coercible object

- label:

  optional column flagging labeled vs. unlabeled rows

- unlabeled_data:

  optional data.frame of unlabeled observations

- na_action:

  how to handle missing data: "na.fail" or "na.omit"

## Value

A list with components data_l and data_u (data.frames)
