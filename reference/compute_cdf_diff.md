# Empirical CDF Difference

Computes the difference between the empirical CDFs of the data and the
predictions.

## Usage

``` r
compute_cdf_diff(Y, f, grid, w = NULL)
```

## Arguments

- Y:

  (matrix): n x 1 matrix of observed data.

- f:

  (matrix): n x 1 matrix of predictions.

- grid:

  (matrix): Grid of values to compute the CDF at.

- w:

  (vector, optional): n-vector of sample weights.

## Value

(list): Difference between the empirical CDFs of the data and the
predictions and its standard deviation at the specified grid points.
