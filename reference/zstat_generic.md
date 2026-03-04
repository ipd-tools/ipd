# Compute Z-Statistic and P-Value

Computes the z-statistic and the corresponding p-value for a given test.

## Usage

``` r
zstat_generic(value1, value2, std_diff, alternative, diff = 0)
```

## Arguments

- value1:

  (numeric): The first value or sample mean.

- value2:

  (numeric): The second value or sample mean.

- std_diff:

  (numeric): The standard error of the difference between the two
  values.

- alternative:

  (character): The alternative hypothesis. Can be one of "two-sided" (or
  "2-sided", "2s"), "larger" (or "l"), or "smaller" (or "s").

- diff:

  (numeric, optional): The hypothesized difference between the two
  values. Default is 0.

## Value

(list): A list containing the following:

- zstat:

  (numeric): The computed z-statistic.

- pvalue:

  (numeric): The corresponding p-value for the test.
