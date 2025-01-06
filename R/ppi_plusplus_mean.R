#===============================================================================
# PPI++ MEAN ESTIMATION
#===============================================================================

#--- PPI++ MEAN ESTIMATION - POINT ESTIMATE ------------------------------------

#' PPI++ Mean Estimation (Point Estimate)
#'
#' @description
#' Helper function for PPI++ mean estimation (point estimate)
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
#' @return float or ndarray: Prediction-powered point estimate of the mean.
#'
#' @examples
#'
#' dat <- simdat(model = "mean")
#'
#' form <- Y - f ~ 1
#'
#' Y_l <- dat[dat$set_label == "labeled",   all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled",   all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_mean_est(Y_l, f_l, f_u)
#'
#' @import stats
#'
#' @export

ppi_plusplus_mean_est <- function(Y_l, f_l, f_u,

  lhat = NULL, coord = NULL, w_l = NULL, w_u = NULL) {

  n <- nrow(Y_l)

  N <- nrow(f_u)

  p <- if (length(dim(f_l)) > 1) dim(f_l)[2] else 1

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  if (is.null(lhat)) {

    est <- mean(w_u * f_u) + mean(w_l * (Y_l - f_l))

    grads <- w_l * (Y_l - est)

    grads_hat <- w_l * (f_l - est)

    grads_hat_unlabeled <- w_u * (f_u - est)

    inv_hessian <- diag(p)

    lhat <- calc_lhat_glm(

      grads, grads_hat, grads_hat_unlabeled, inv_hessian, coord, clip = T)

    return(ppi_plusplus_mean_est(Y_l, f_l, f_u, lhat, coord, w_l, w_u))

  } else {

    return(mean(w_u * lhat * f_u) + mean(w_l * (Y_l - lhat * f_l)))
  }
}

#--- PPI++ MEAN ESTIMATION - INFERENCE -----------------------------------------

#' PPI++ Mean Estimation
#'
#' @description
#' Helper function for PPI++ mean estimation
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
#' @param alpha (scalar): type I error rate for hypothesis testing - values in
#' (0, 1); defaults to 0.05.
#'
#' @param alternative (string): Alternative hypothesis. Must be one of
#' \code{"two-sided"}, \code{"less"}, or \code{"greater"}.
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
#' @return tuple: Lower and upper bounds of the prediction-powered confidence
#' interval for the mean.
#'
#' @examples
#'
#' dat <- simdat(model = "mean")
#'
#' form <- Y - f ~ 1
#'
#' Y_l <- dat[dat$set_label == "labeled",   all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled",   all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_plusplus_mean(Y_l, f_l, f_u)
#'
#' @import stats
#'
#' @export

ppi_plusplus_mean <- function(Y_l, f_l, f_u,

  alpha = 0.05, alternative = "two-sided",

  lhat = NULL, coord = NULL, w_l = NULL, w_u = NULL) {

  #- Compute Dimensions of Inputs

  n <- ifelse(is.null(dim(Y_l)), length(Y_l), nrow(Y_l))

  N <- ifelse(is.null(dim(f_u)), length(f_u), nrow(f_u))

  p <- if (length(dim(f_l)) > 1) dim(f_l)[2] else 1

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  if (is.null(lhat)) {

    est <- ppi_plusplus_mean_est(Y_l, f_l, f_u, 1, w_l, w_u)

    grads <- w_l * (Y_l - est)

    grads_hat <- w_l * (f_l - est)

    grads_hat_unlabeled <- w_u * (f_u - est)

    inv_hessian <- diag(p)

    lhat <- calc_lhat_glm(

      grads, grads_hat, grads_hat_unlabeled, inv_hessian, coord, clip = T)

    return(

      ppi_plusplus_mean(Y_l, f_l, f_u,

        alpha, alternative,lhat, coord, w_l, w_u))
  }

  est <- ppi_plusplus_mean_est(Y_l, f_l, f_u, lhat, coord, w_l, w_u)

  imputed_std <- sd(w_u * (lhat * f_u)) * sqrt((N - 1) / N) / sqrt(N)

  rectifier_std <- sd(w_l * (Y_l - lhat * f_l)) * sqrt((n - 1) / n) / sqrt(n)

  return(zconfint_generic(

    est, sqrt(imputed_std^2 + rectifier_std^2), alpha, alternative))
}

#=== END =======================================================================
