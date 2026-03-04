# Empirical CDF of the Data

Computes the empirical CDF of the data.

## Usage

``` r
compute_cdf(Y, grid, w = NULL)
```

## Arguments

- Y:

  (matrix): n x 1 matrix of observed data.

- grid:

  (matrix): Grid of values to compute the CDF at.

- w:

  (vector, optional): n-vector of sample weights.

## Value

(list): Empirical CDF and its standard deviation at the specified grid
points.
