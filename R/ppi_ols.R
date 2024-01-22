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
#' @param n (int): Number of labeled observations.
#'
#' @param p (int): Number of covariates (features) of interest.
#'
#' @param N (int): Number of unlabeled observations.
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
#' X_l <- model.matrix(form, data = dat[dat$set == "tst",])
#'
#' Y_l <- dat[dat$set == "tst", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "tst", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set == "val",])
#'
#' f_u <- dat[dat$set == "val", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' n <- nrow(X_l)
#'
#' p <- ncol(X_l)
#'
#' N <- nrow(X_u)
#'
#' ppi_ols_est(X_l, Y_l, f_l, X_u, f_u, n, p, N)
#'
#' @import stats
#'
#' @export

ppi_ols_est <- function(X_l, Y_l, f_l, X_u, f_u, n, p, N,

                             lhat = NULL, coord = NULL, w_l = NULL, w_u = NULL) {

  X_u <- as.matrix(X_u)

  n <- nrow(X_l)
  N <- nrow(X_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  use_u <- is.null(lhat) || lhat != 0

  if (is.null(lhat)) {

    imputed_theta <- wls(X_u, f_u, w = w_u)

    rectifier <- wls(X_l, Y_l - f_l, w = w_l)

    est <- imputed_theta + rectifier

    stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u, w_l, w_u, use_u)

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

                          stats$grads_hat_unlabeled, stats$inv_hessian, coord, clip = T)

    return(ppi_ols_est(X_l, Y_l, f_l, X_u, f_u,

                            lhat = lhat, coord = coord, w_l = w_l, w_u = w_u))

  } else {

    imputed_theta <- wls(X_u, lhat * f_u, w = w_u)

    rectifier <- wls(X_l, Y_l - lhat * f_l, w = w_l)

    est <- imputed_theta + rectifier

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
#' @param n (int): Number of labeled observations.
#'
#' @param p (int): Number of covariates (features) of interest.
#'
#' @param N (int): Number of unlabeled observations.
#'
#' @param alpha (float): Significance level in \[0,1\]
#'
#' @param alternative (string): Alternative hypothesis, either 'two-sided',
#' 'larger' or 'smaller'.
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
#' X_l <- model.matrix(form, data = dat[dat$set == "tst",])
#'
#' Y_l <- dat[dat$set == "tst", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "tst", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set == "val",])
#'
#' f_u <- dat[dat$set == "val", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' n <- nrow(X_l)
#'
#' p <- ncol(X_l)
#'
#' N <- nrow(X_u)
#'
#' ppi_ols(X_l, Y_l, f_l, X_u, f_u, n, p, N)
#'
#' @import stats
#'
#' @export

ppi_ols <- function(X_l, Y_l, f_l, X_u, f_u, n, p, N,

                         alpha = 0.05, alternative = "two-sided", lhat = NULL,

                         coord = NULL, w_l = NULL, w_u = NULL) {

  X_l <- as.matrix(X_l)                                                          # Need in wrapper?
  X_u <- as.matrix(X_u)                                                          # Need in wrapper?

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  use_u <- is.null(lhat) || lhat != 0

  est <- ppi_ols_est(

    X_l, Y_l, f_l, X_u, f_u, n, p, N, lhat, coord, w_l, w_u)

  stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u, w_l, w_u, use_u)

  inv_hessian <- stats$inv_hessian

  if (is.null(lhat)) {

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

                          stats$grads_hat_unlabeled, inv_hessian, coord, clip = T)

    return(ppi_ols(X_l, Y_l, f_l, X_u, f_u, n, p, N,

                        alpha = alpha, alternative = alternative, lhat = lhat, coord = coord,

                        w_l = w_l, w_u = w_u))
  }

  var_u <- cov(lhat * stats$grads_hat_unlabeled)

  var_l <- cov(stats$grads - lhat * stats$grads_hat)

  Sigma_hat <- inv_hessian %*% (n/N * var_u + var_l) %*% inv_hessian

  return(list(est = est, se = sqrt(diag(Sigma_hat) / n)))
}

#=== END =======================================================================
