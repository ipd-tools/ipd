# Rectified P-Value

Computes a rectified p-value.

## Usage

``` r
rectified_p_value(
  rectifier,
  rectifier_std,
  imputed_mean,
  imputed_std,
  null = 0,
  alternative = "two-sided"
)
```

## Arguments

- rectifier:

  (float or vector): Rectifier value.

- rectifier_std:

  (float or vector): Rectifier standard deviation.

- imputed_mean:

  (float or vector): Imputed mean.

- imputed_std:

  (float or vector): Imputed standard deviation.

- null:

  (float, optional): Value of the null hypothesis to be tested. Defaults
  to `0`.

- alternative:

  (str, optional): Alternative hypothesis, either 'two-sided', 'larger'
  or 'smaller'.

## Value

(float or vector): The rectified p-value.
