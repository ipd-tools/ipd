# Normal Confidence Intervals

Calculates normal confidence intervals for a given alternative at a
given significance level.

## Usage

``` r
zconfint_generic(mean, std_mean, alpha, alternative)
```

## Arguments

- mean:

  (float): Estimated normal mean.

- std_mean:

  (float): Estimated standard error of the mean.

- alpha:

  (float): Significance level in \[0,1\]

- alternative:

  (string): Alternative hypothesis, either 'two-sided', 'larger' or
  'smaller'.

## Value

(vector): Lower and upper (1 - alpha) \* 100% confidence limits.
