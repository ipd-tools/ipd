#===============================================================================
# PPI++ ORDINARY LEAST SQUARES
#===============================================================================

#--- PPI++ OLS - POINT ESTIMATE ------------------------------------------------

#' PPI++ OLS (Point Estimate)
#'
#' @description
#' Helper function for PPI++ OLS estimation (point estimate)
#'
#' @details
#' PPI++: Efficient Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://arxiv.org/abs/2311.01453>
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
#' @param lhat (float, optional): Power-tuning parameter (see
#' <https://arxiv.org/abs/2311.01453>). The default value, \code{NULL},
#' will estimate the optimal value from the data. Setting \code{lhat = 1}
#' recovers PPI with no power tuning, and setting \code{lhat = 0} recovers
#' the classical point estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize
#' \code{lhat = 1}. If \code{NULL}, it optimizes the total variance over all
#' coordinates. Must be in (1, ..., d) where d is the dimension of the estimand.
#'
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to a vector of ones.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to a vector of ones.
#'
#' @return (vector): vector of prediction-powered point estimates of the OLS
#' coefficients.
#'
#' @examples
#'
#' dat <- simdat(model = "ols")
#'
#' form <- Y - f ~ X1
#'
#' X_l <- model.matrix(form, data = dat[dat$set_label == "labeled",])
#'
#' Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled",])
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
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

#--- PPI++ OLS - INFERENCE -----------------------------------------------------

#' PPI++ OLS
#'
#' @description
#' Helper function for PPI++ OLS estimation
#'
#' @details
#' PPI++: Efficient Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://arxiv.org/abs/2311.01453>`
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
#' @param lhat (float, optional): Power-tuning parameter (see
#' <https://arxiv.org/abs/2311.01453>). The default value, \code{NULL},
#' will estimate the optimal value from the data. Setting \code{lhat = 1}
#' recovers PPI with no power tuning, and setting \code{lhat = 0} recovers
#' the classical point estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize
#' \code{lhat = 1}. If \code{NULL}, it optimizes the total variance over all
#' coordinates. Must be in (1, ..., d) where d is the dimension of the estimand.
#'
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to a vector of ones.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to a vector of ones.
#'
#' @return (list): A list containing the following:
#'
#' \describe{
#'    \item{est}{(vector): vector of PPI++ OLS regression coefficient
#'    estimates.}
#'    \item{se}{(vector): vector of standard errors of the coefficients.}
#'    \item{lambda}{(float): estimated power-tuning parameter.}
#'    \item{rectifier_est}{(vector): vector of the rectifier OLS
#'    regression coefficient estimates.}
#'    \item{var_u}{(matrix): covariance matrix for the gradients in the
#'    unlabeled data.}
#'    \item{var_l}{(matrix): covariance matrix for the gradients in the
#'    labeled data.}
#'    \item{grads}{(matrix): matrix of gradients for the
#'    labeled data.}
#'    \item{grads_hat_unlabeled}{(matrix): matrix of predicted gradients for
#'    the unlabeled data.}
#'    \item{grads_hat}{(matrix): matrix of predicted gradients for the
#'    labeled data.}
#'    \item{inv_hessian}{(matrix): inverse Hessian matrix.}
#' }
#'
#' @examples
#'
#' dat <- simdat(model = "ols")
#'
#' form <- Y - f ~ X1
#'
#' X_l <- model.matrix(form, data = dat[dat$set_label == "labeled",])
#'
#' Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled",])
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_plusplus_ols <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, w_l = NULL, w_u = NULL) {

  n <- NROW(f_l)

  N <- NROW(f_u)

  p <- NCOL(X_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  use_u <- is.null(lhat) || lhat != 0

  delta_hat <- wls(X_l, Y_l - f_l, w = w_l)

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

  return(list(est = est, se = sqrt(diag(Sigma_hat) / n), lambda = lhat,

    rectifier_est = delta_hat, var_u = var_u, var_l = var_l,

    grads = stats$grads, grads_hat_unlabeled = stats$grads_hat_unlabeled,

    grads_hat = stats$grads_hat, inv_hessian = stats$inv_hessian))
}

#=== END =======================================================================
