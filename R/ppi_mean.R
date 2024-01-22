#===============================================================================
# PPI++ MEAN ESTIMATION
#===============================================================================

#=== MEAN ESTIMATION ===========================================================

#' ...Need title...
#'
#' @description
#' Computes the prediction-powered point estimate of the p-dimensional mean.
#'
#' @details
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__ A. N. Angelopoulos, J. C.
#' Duchi, and T. Zrnic. PPI++: Efficient Prediction Powered Inference.
#' arxiv:2311.01453, 2023.
#'
#' @param Y_l (ndarray): Gold-standard labels.
#'
#' @param f_l (ndarray): Predictions corresponding to the gold-standard
#' labels.
#'
#' @param f_u (ndarray): Predictions corresponding to the unlabeled
#' data.
#'
#' @param lhat (float, optional): Power-tuning parameter (see
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None`
#' will estimate the optimal value from data. Setting `lhat=1` recovers PPI
#' with no power tuning, and setting `lhat=0` recovers the classical point
#' estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`. If
#' `None`, it optimizes the total variance over all coordinates. Must be in
#' {1, ..., d} where d is the dimension of the estimand.
#'
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to all ones vector.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to all ones vector.
#'
#' @returns float or ndarray: Prediction-powered point estimate of the mean.
#'
#' @examples
#'
#' #need examples
#'
#' @import stats
#'
#' @export

ppi_mean_est <- function(Y_l, f_l, f_u, lhat = NULL, coord = NULL,

                              w_l = NULL, w_u = NULL) {

  n <- ifelse(is.null(dim(Y_l)), length(Y_l), nrow(Y_l))
  N <- ifelse(is.null(dim(f_u)), length(f_u), nrow(f_u))
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

    return(ppi_mean_est(Y_l, f_l, f_u, lhat, coord, w_l, w_u))

  } else {

    return(mean(w_u * lhat * f_u) + mean(w_l * (Y_l - lhat * f_l)))
  }
}

#' ...Need title...
#'
#' @description
#' Computes the prediction-powered confidence interval for a d-dimensional mean.
#'
#' @details
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__ A. N. Angelopoulos, J. C.
#' Duchi, and T. Zrnic. PPI++: Efficient Prediction Powered Inference.
#' arxiv:2311.01453, 2023.
#'
#' @param Y_l (ndarray): Gold-standard labels.
#'
#' @param f_l (ndarray): Predictions corresponding to the gold-standard
#' labels.
#'
#' @param f_u (ndarray): Predictions corresponding to the unlabeled
#' data.
#'
#' @param alpha (float, optional): Error level; the confidence interval will target a coverage of 1 - alpha. Must be in (0, 1).
#'
#' @param alternative (str, optional): Alternative hypothesis, either 'two-sided', 'larger' or 'smaller'.
#'
#' @param lhat (float, optional): Power-tuning parameter (see
#' `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None`
#' will estimate the optimal value from data. Setting `lhat=1` recovers PPI
#' with no power tuning, and setting `lhat=0` recovers the classical point
#' estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`. If
#' `None`, it optimizes the total variance over all coordinates. Must be in
#' {1, ..., d} where d is the dimension of the estimand.
#'
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to all ones vector.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to all ones vector.
#'
#' @returns tuple: Lower and upper bounds of the prediction-powered confidence
#' interval for the mean.
#'
#' @examples
#'
#' #need examples
#'
#' @import stats
#'
#' @export

ppi_mean <- function(Y_l, f_l, f_u, alpha = 0.05, alternative = "two-sided",

                          lhat = NULL, coord = NULL, w_l = NULL, w_u = NULL) {

  #- Compute Dimensions of Inputs

  n <- ifelse(is.null(dim(Y_l)), length(Y_l), nrow(Y_l))
  N <- ifelse(is.null(dim(f_u)), length(f_u), nrow(f_u))
  p <- if (length(dim(f_l)) > 1) dim(f_l)[2] else 1

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  if (is.null(lhat)) {

    est <- ppi_mean_est(Y_l, f_l, f_u, 1, w_l, w_u)

    grads <- w_l * (Y_l - est)

    grads_hat <- w_l * (f_l - est)

    grads_hat_unlabeled <- w_u * (f_u - est)

    inv_hessian <- diag(p)

    lhat <- calc_lhat_glm(

      grads, grads_hat, grads_hat_unlabeled, inv_hessian, coord, clip = T)

    return(ppi_mean(Y_l, f_l, f_u, lhat, coord, w_l, w_u))
  }

  est <- ppi_mean_est(Y_l, f_l, f_u, lhat, coord, w_l, w_u)

  imputed_std <- sd(w_u * (lhat * f_u)) * sqrt((N - 1) / N) / sqrt(N)

  rectifier_std <- sd(w_l * (Y_l - lhat * f_l)) * sqrt((n - 1) / n) / sqrt(n)

  return(zconfint_generic(

    est, sqrt(imputed_std^2 + rectifier_std^2), alpha, alternative))
}

#=== END =======================================================================
