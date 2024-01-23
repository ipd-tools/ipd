#===============================================================================
# PPI++ ORDINARY LEAST SQUARES
#===============================================================================

#=== PPI++ OLS POINT ESTIMATE ==================================================

#' PPI++ OLS Point Estimate (helper function)
#'
#' @description
#' Computes the prediction-powered point estimate of the OLS coefficients.
#'
#' @param X_l (matrix): n x p matrix of covariates in the labeled data.
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param X_u (matrix): N x p matrix of covariates in the unlabeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param lhat (float, optional): Power-tuning parameter.
#' The default value `NULL` will estimate the optimal value from data.
#' Setting `lhat = 1` recovers PPI with no power tuning.
#' Setting `lhat = 0` recovers the classical point estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`.
#' If `None`, it optimizes the total variance over all coordinates.
#' Must be in {1, ..., p} where p is the shape of the estimand.
#'
#' @param w_l (vector, optional): n-vector of sample weights for the
#' labeled data.
#'
#' @param w_u (vector, optional): N-vector of sample weights for the
#' unlabeled data.
#'
#' @returns (vector): p-vector of prediction-powered point estimates of the
#' OLS coefficients.
#'
#' @examples
#'
#' dat <- simdat()
#'
#' form <- Y - Yhat ~ X1
#'
#' X_l <- model.matrix(form, data = dat[dat$set == "labeled",])
#'
#' Y_l <- dat[dat$set == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set == "unlabeled",])
#'
#' f_u <- dat[dat$set == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_ols_est(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_plusplus_ols_est <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, w_l = NULL, w_u = NULL) {

  n <- nrow(f_l)

  N <- nrow(f_u)

  p <- ncol(X_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  use_u <- is.null(lhat) || lhat != 0

  if (is.null(lhat)) {

    theta_hat <- wls(X_u, f_u, w = w_u)

    delta_hat <- wls(X_l, Y_l - f_l, w = w_l)

    est <- theta_hat + delta_hat

    stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u, w_l, w_u, use_u)

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, coord, clip = T)

    return(ppi_plusplus_ols_est(X_l, Y_l, f_l, X_u, f_u,

      lhat = lhat, coord = coord, w_l = w_l, w_u = w_u))

  } else {

    theta_hat <- wls(X_u, lhat * f_u, w = w_u)

    delta_hat <- wls(X_l, Y_l - lhat * f_l, w = w_l)

    est <- theta_hat + delta_hat

    return(est)
  }
}

#=== PPI++ OLS =================================================================

#' PPI++ OLS Estimator and Inference
#'
#' @description
#' Computes the prediction-powered estimator and confidence interval for the
#' OLS coefficients using the PPI++ algorithm.
#'
#' @param X_l (matrix): n x p matrix of covariates in the labeled data.
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param X_u (matrix): N x p matrix of covariates in the unlabeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param lhat (float, optional): Power-tuning parameter.
#' The default value `NULL` will estimate the optimal value from data.
#' Setting `lhat = 1` recovers PPI with no power tuning.
#' Setting `lhat = 0` recovers the classical point estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`.
#' If `None`, it optimizes the total variance over all coordinates.
#' Must be in {1, ..., p} where p is the shape of the estimand.
#'
#' @param w_l (vector, optional): n-vector of sample weights for the
#' labeled data.
#'
#' @param w_u (vector, optional): N-vector of sample weights for the
#' unlabeled data.
#'
#' @returns (list): A list containing the following:
#'
#' \describe{
#'    \item{est}{(vector): p-vector of PPI++ OLS coefficient estimates.}
#'    \item{se}{(vector): p-vector of standard errors of the coefficients.}
#' }
#'
#' @examples
#'
#' dat <- simdat()
#'
#' form <- Y - Yhat ~ X1
#'
#' X_l <- model.matrix(form, data = dat[dat$set == "labeled",])
#'
#' Y_l <- dat[dat$set == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set == "unlabeled",])
#'
#' f_u <- dat[dat$set == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_plusplus_ols <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, w_l = NULL, w_u = NULL) {

  n <- nrow(f_l)

  N <- nrow(f_u)

  p <- ncol(X_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  use_u <- is.null(lhat) || lhat != 0

  est <- ppi_plusplus_ols_est(X_l, Y_l, f_l, X_u, f_u, lhat, coord, w_l, w_u)

  stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u, w_l, w_u, use_u)

  if (is.null(lhat)) {

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, coord, clip = T)

    return(ppi_plusplus_ols(X_l, Y_l, f_l, X_u, f_u,

      lhat = lhat, coord = coord, w_l = w_l, w_u = w_u))
  }

  var_u <- cov(lhat * stats$grads_hat_unlabeled)

  var_l <- cov(stats$grads - lhat * stats$grads_hat)

  Sigma_hat <- stats$inv_hessian %*% (n/N * var_u + var_l) %*% stats$inv_hessian

  return(list(est = est, se = sqrt(diag(Sigma_hat) / n)))
}

#=== END =======================================================================
