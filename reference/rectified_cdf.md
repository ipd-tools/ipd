# Rectified CDF

Computes the rectified CDF of the data.

## Usage

``` r
rectified_cdf(Y_l, f_l, f_u, grid, w_l = NULL, w_u = NULL)
```

## Arguments

- Y_l:

  (vector): Gold-standard labels.

- f_l:

  (vector): Predictions corresponding to the gold-standard labels.

- f_u:

  (vector): Predictions corresponding to the unlabeled data.

- grid:

  (vector): Grid of values to compute the CDF at.

- w_l:

  (vector, optional): Sample weights for the labeled data set.

- w_u:

  (vector, optional): Sample weights for the unlabeled data set.

## Value

(vector): Rectified CDF of the data at the specified grid points.
