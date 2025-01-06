#===============================================================================
# PPI++ QUANTILE ESTIMATION
#===============================================================================

#--- PPI++ QUANTILE ESTIMATION - POINT ESTIMATE --------------------------------

#' PPI++ Quantile Estimation (Point Estimate)
#'
#' @description
#' Helper function for PPI++ quantile estimation (point estimate)
#'
#' @details
#' PPI++: Efficient Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://arxiv.org/abs/2311.01453>`
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param q (float): Quantile to estimate. Must be in the range (0, 1).
#'
#' @param exact_grid (bool, optional): Whether to compute the exact solution
#' (TRUE) or an approximate solution based on a linearly spaced grid of 5000
#' values (FALSE).
#'
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to a vector of ones.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to a vector of ones.
#'
#' @return (float): Prediction-powered point estimate of the quantile.
#'
#' @examples
#'
#' dat <- simdat(model = "quantile")
#'
#' form <- Y - f ~ 1
#'
#' Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_quantile_est(Y_l, f_l, f_u, q = 0.5)
#'
#' @export

ppi_plusplus_quantile_est <- function(Y_l, f_l, f_u, q,

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

  rect_cdf <- rectified_cdf(Y_l, f_l, f_u, grid, w_l = w_l, w_u = w_u)

  minimizer <- which.min(abs(rect_cdf - q))

  return(grid[minimizer])
}

#--- PPI++ QUANTILE ESTIMATION - INFERENCE -------------------------------------

#' PPI++ Quantile Estimation
#'
#' @description
#' Helper function for PPI++ quantile estimation
#'
#' @details
#' PPI++: Efficient Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://arxiv.org/abs/2311.01453>
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param q (float): Quantile to estimate. Must be in the range (0, 1).
#'
#' @param alpha (scalar): type I error rate for hypothesis testing - values in
#' (0, 1); defaults to 0.05.
#'
#' @param exact_grid (bool, optional): Whether to compute the exact solution
#' (TRUE) or an approximate solution based on a linearly spaced grid of 5000
#' values (FALSE).
#'
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to a vector of ones.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to a vector of ones.
#'
#' @return tuple: Lower and upper bounds of the prediction-powered confidence
#' interval for the quantile.
#'
#' @examples
#'
#' dat <- simdat(model = "quantile")
#'
#' form <- Y - f ~ X1
#'
#' Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_quantile(Y_l, f_l, f_u, q = 0.5)
#'
#' @export

ppi_plusplus_quantile <- function(Y_l, f_l, f_u, q,

  alpha = 0.05, exact_grid = FALSE, w_l = NULL, w_u = NULL) {

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
