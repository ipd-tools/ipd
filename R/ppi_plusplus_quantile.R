#===============================================================================
# PPI++ QUANTILE ESTIMATION
#===============================================================================

#=== EMPIRICAL CDF =============================================================

#' Empirical CDF of the Data
#'
#' Computes the empirical CDF of the data.
#'
#' @param Y (matrix): n x 1 matrix of observed data.
#'
#' @param grid (matrix): Grid of values to compute the CDF at.
#'
#' @param w (vector, optional): n-vector of sample weights.
#'
#' @returns (list): Empirical CDF and its standard deviation at the specified
#' grid points.
#'
#' @examples
#' # Need example code
#'
#' @export

compute_cdf <- function(Y, grid, w = NULL) {

  n <- length(Y)

  if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

  indicators <- matrix((Y <= rep(grid, each = n)) * w,

                       ncol = length(grid))

  cdf_mn <- apply(indicators, 2, mean)

  cdf_sd <- apply(indicators, 2, sd) * sqrt((n - 1) / n)

  return(list(cdf_mn, cdf_sd))
}

#=== DIFFERENCE IN EMPIRICAL CDFS ==============================================

#' Empirical CDF Difference
#'
#' Computes the difference between the empirical CDFs of the data and the
#' predictions.
#'
#' @param Y (matrix): n x 1 matrix of observed data.
#'
#' @param Yhat (matrix): n x 1 matrix of predictions.
#'
#' @param grid (matrix): Grid of values to compute the CDF at.
#'
#' @param w (vector, optional): n-vector of sample weights.
#'
#' @returns (list): Difference between the empirical CDFs of the data and the
#' predictions and its standard deviation at the specified grid points.
#'
#' @examples
#' # Need example code
#'
#' @export

compute_cdf_diff <- function(Y, Yhat, grid, w = NULL) {

  n <- length(Y)

  if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

  indicators_Y <- matrix((Y <= rep(grid, each = n)) * w,

                         ncol = length(grid))

  indicators_Yhat <- matrix((Yhat <= rep(grid, each = length(Yhat))) * w,

                            ncol = length(grid))

  diff_mn <- apply(indicators_Y - indicators_Yhat, 2, mean)

  diff_sd <- apply(indicators_Y - indicators_Yhat, 2, sd) * sqrt((n - 1) / n)

  return(list(diff_mn, diff_sd))
}

#=== RECTIFIED CDF =============================================================

#' Rectified CDF
#'
#' Computes the rectified CDF of the data.
#'
#' @param Y_l (vector): Gold-standard labels.
#'
#' @param f_l (vector): Predictions corresponding to the gold-standard labels.
#'
#' @param f_u (vector): Predictions corresponding to the unlabeled data.
#'
#' @param grid (vector): Grid of values to compute the CDF at.
#'
#' @param w_l (vector, optional): Sample weights for the labeled data set.
#'
#' @param w_u (vector, optional): Sample weights for the unlabeled data set.
#'
#' @returns (vector): Rectified CDF of the data at the specified grid points.
#'
#' @examples
#' # Need example code
#'
#' @export

rectified_cdf <- function(Y_l, f_l, f_u, grid, w_l = NULL, w_u = NULL) {

  n <- length(Y_l)

  N <- length(f_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  cdf_f_u <- compute_cdf(f_u, grid, w = w_u)[[1]]

  cdf_rectifier <- compute_cdf_diff(Y_l, f_l, grid, w = w_l)[[1]]

  return(cdf_f_u + cdf_rectifier)
}

#=== RECTIFIED P-VALUE =========================================================

#' Rectified P-Value
#'
#' Computes a rectified p-value.
#'
#' @param rectifier (float or vector): Rectifier value.
#'
#' @param rectifier_std (float or vector): Rectifier standard deviation.
#'
#' @param imputed_mean (float or vector): Imputed mean.
#'
#' @param imputed_std (float or vector): Imputed standard deviation.
#'
#' @param null (float, optional): Value of the null hypothesis to be tested.
#' Defaults to `0`.
#'
#' @param alternative (str, optional): Alternative hypothesis, either
#' 'two-sided', 'larger' or 'smaller'.
#'
#' @examples
#' # Need example code
#'
#'
#' @returns (float or vector): P-value.

rectified_p_value <- function(rectifier, rectifier_std,

                              imputed_mean, imputed_std, null = 0, alternative = "two-sided") {

  rectified_point_estimate <- imputed_mean + rectifier

  rectified_std <- pmax(sqrt(imputed_std^2 + rectifier_std^2), 1e-16)

  p_value <- zstat_generic(rectified_point_estimate, 0, rectified_std,

                           alternative, null)[[2]]

  return(p_value)
}

#=== PPI++ QUANTILE ESTIMATION POINT ESTIMATE ==================================

#' PPI++ Quantile Estimation Point Estimate
#'
#' Computes the prediction-powered point estimate of the quantile.
#'
#' @param Y_l (vector): Gold-standard labels.
#'
#' @param f_l (vector): Predictions corresponding to the gold-standard labels.
#'
#' @param f_u (vector): Predictions corresponding to the unlabeled data.
#'
#' @param q (float): Quantile to estimate.
#'
#' @param exact_grid (bool, optional): Whether to compute the exact solution
#' (TRUE) or an approximate solution based on a linearly spaced grid of 5000
#' values (FALSE).
#'
#' @param w_l (vector, optional): Sample weights for the labeled data set.
#'
#' @param w_u (vector, optional): Sample weights for the unlabeled data set.
#'
#' @returns (float): Prediction-powered point estimate of the quantile.
#'
#' @examples
#'
#' dat <- simdat()
#'
#' form <- Y - Yhat ~ X1
#'
#' Y_l <- dat[dat$set == "tst", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "tst", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set == "val", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_quantile_est(Y_l, f_l, f_u, q = 0.5)
#'
#' @export

ppi_plusplus_quantile_est <- function(Y_l, f_l, f_u, q, exact_grid = FALSE,

                                  w_l = NULL, w_u = NULL) {

  Y_l <- c(Y_l)

  f_l <- c(f_l)

  f_u <- c(f_u)

  n <- length(Y_l)

  N <- length(f_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  grid <- c(Y_l, f_l, f_u)

  if (exact_grid) {

    grid <- sort(unique(grid))

  } else {

    grid <- seq(min(grid), max(grid), length.out = 5000)
  }

  rect_cdf <- rectified_cdf(Y_l, f_l, f_u, grid, w_l = w_l, w_u = w_u)

  minimizer <- which.min(abs(rect_cdf - q))

  return(grid[minimizer])
}

#=== PPI++ QUANTILE ESTIMATION =================================================

#' PPI++ Quantile Estimation
#'
#' Computes the prediction-powered confidence interval for the quantile.
#'
#' @param Y_l (vector): Gold-standard labels.
#'
#' @param f_l (vector): Predictions corresponding to the gold-standard labels.
#'
#' @param f_u (vector): Predictions corresponding to the unlabeled data.
#'
#' @param q (float): Quantile to estimate. Must be in the range (0, 1).
#'
#' @param alpha (float, optional): Error level; the confidence interval will
#' target a coverage of 1 - alpha. Must be in the range (0, 1).
#'
#' @param exact_grid (bool, optional): Whether to use the exact grid of values
#' (TRUE) or a linearly spaced grid of 5000 values (FALSE).
#'
#' @param w_l (vector, optional): Sample weights for the labeled data set.
#'
#' @param w_u (vector, optional): Sample weights for the unlabeled data set.
#'
#' @returns (list): Lower and upper bounds of the prediction-powered confidence  # Fix this to be estimate and sd
#' interval for the quantile.
#'
#' @examples
#'
#' dat <- simdat()
#'
#' form <- Y - Yhat ~ X1
#'
#' Y_l <- dat[dat$set == "tst", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "tst", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set == "val", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_quantile(Y_l, f_l, f_u, q = 0.5)
#'
#' @export

ppi_plusplus_quantile <- function(Y_l, f_l, f_u, q, alpha = 0.05,

                              exact_grid = FALSE, w_l = NULL, w_u = NULL) {

  Y_l <- c(Y_l)

  f_l <- c(f_l)

  f_u <- c(f_u)

  n <- length(Y_l)

  N <- length(f_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  grid <- c(Y_l, f_l, f_u)

  if (exact_grid) {

    grid <- sort(unique(grid))

  } else {

    grid <- seq(min(grid), max(grid), length.out = 5000)
  }

  cdf_f_u <- compute_cdf(f_u, grid, w_u)

  cdf_rectifier <- compute_cdf_diff(Y_l, f_l, grid, w_l)

  rec_p_value <- rectified_p_value(

    cdf_rectifier[[1]], cdf_rectifier[[2]] / sqrt(n),

    cdf_f_u[[1]], cdf_f_u[[2]] / sqrt(N),

    null = q, alternative = "two-sided")

  result_grid <- grid[rec_p_value > alpha]

  return(c(result_grid[1], result_grid[length(result_grid)]))
}

#=== END =======================================================================
