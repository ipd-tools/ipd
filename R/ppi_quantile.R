#===============================================================================
# PPI QUANTILE ESTIMATION
#===============================================================================

#--- PPI QUANTILE ESTIMATION ---------------------------------------------------

#' PPI Quantile Estimation
#'
#' @description
#' Helper function for PPI quantile estimation
#'
#' @details
#' Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://www.science.org/doi/10.1126/science.adi6000>
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
#' ppi_quantile(Y_l, f_l, f_u, q = 0.5)
#'
#' @export

ppi_quantile <- function(Y_l, f_l, f_u, q,

  alpha = 0.05, exact_grid = FALSE) {

  Y_l <- c(Y_l)

  f_l <- c(f_l)

  f_u <- c(f_u)

  n <- length(Y_l)

  N <- length(f_u)

  grid <- c(Y_l, f_l, f_u)

  if (exact_grid) {

    grid <- sort(unique(grid))

  } else {

    grid <- seq(min(grid), max(grid), length.out = 5000)
  }

  cdf_f_u <- compute_cdf(f_u, grid, w = NULL)

  cdf_rectifier <- compute_cdf_diff(Y_l, f_l, grid, w = NULL)

  rec_p_value <- rectified_p_value(

    cdf_rectifier[[1]], cdf_rectifier[[2]] / sqrt(n),

    cdf_f_u[[1]], cdf_f_u[[2]] / sqrt(N),

    null = q, alternative = "two-sided")

  result_grid <- grid[rec_p_value > alpha]

  return(c(result_grid[1], result_grid[length(result_grid)]))
}

#=== END =======================================================================
