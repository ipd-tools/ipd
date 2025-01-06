#===============================================================================
# PPI++ LOGISTIC REGRESSION
#===============================================================================

#--- PPI++ LOGISTIC REGRESSION - POINT ESTIMATE --------------------------------

#' PPI++ Logistic Regression (Point Estimate)
#'
#' @description
#' Helper function for PPI++ logistic regression (point estimate)
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
#' @param opts (list, optional): Options to pass to the optimizer.
#' See ?optim for details.
#'
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to a vector of ones.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to a vector of ones.
#'
#' @return (vector): vector of prediction-powered point estimates of the
#' logistic regression coefficients.
#'
#' @examples
#'
#' dat <- simdat(model = "logistic")
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
#' ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_plusplus_logistic_est <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, opts = NULL, w_l = NULL, w_u = NULL) {

  n <- nrow(f_l)

  N <- nrow(f_u)

  p <- ncol(X_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  if (is.null(opts) || !("factr" %in% names(opts))) {

    opts <- list(factr = 1e-15)
  }

  theta <- coef(glm(Y_l ~ . - 1,

    data = data.frame(Y_l, X_l), family = binomial))

  theta <- matrix(theta, ncol = 1)

  lhat_curr <- ifelse(is.null(lhat), 1, lhat)

  #-- Rectified Logistic Regression Loss Function

  rectified_logistic_loss <- function(theta) {

    sum(w_u * (-f_u * (X_u %*% theta) + log1pexp(X_u %*% theta))) *

      lhat_curr / N - sum(

        w_l * (-f_l * (X_l %*% theta) + log1pexp(X_l %*% theta))) *

      lhat_curr / n + sum(

        w_l * (-Y_l * (X_l %*% theta) + log1pexp(X_l %*% theta))) / n
  }

  #-- Rectified Logistic Regression Gradient

  rectified_logistic_grad <- function(theta) {

    lhat_curr / N * t(X_u) %*% (w_u * (plogis(X_u %*% theta) - f_u)) -

      lhat_curr / n * t(X_l) %*% (w_l * (plogis(X_l %*% theta) - f_l)) +

      1 / n * t(X_l) %*% (w_l * (plogis(X_l %*% theta) - Y_l))
  }

  est <- optim(par = theta, fn = rectified_logistic_loss,

    gr = rectified_logistic_grad, method = "L-BFGS-B",

    control = list(factr = opts$factr))$par

  if (is.null(lhat)) {

    stats <- logistic_get_stats(est, X_l, Y_l, f_l, X_u, f_u, w_l, w_u)

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, clip = TRUE)

    return(ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u,

      opts = opts, lhat = lhat, coord = coord, w_l = w_l, w_u = w_u))

  } else {

    return(est)
  }
}

#--- PPI++ LOGISTIC REGRESSION - INFERENCE -------------------------------------

#' PPI++ Logistic Regression
#'
#' @description
#' Helper function for PPI++ logistic regression
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
#' @param opts (list, optional): Options to pass to the optimizer.
#' See ?optim for details.
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
#'    \item{est}{(vector): vector of PPI++ logistic regression coefficient
#'    estimates.}
#'    \item{se}{(vector): vector of standard errors of the coefficients.}
#'    \item{lambda}{(float): estimated power-tuning parameter.}
#'    \item{rectifier_est}{(vector): vector of the rectifier logistic
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
#' dat <- simdat(model = "logistic")
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
#' ppi_plusplus_logistic(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_plusplus_logistic <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, opts = NULL, w_l = NULL, w_u = NULL) {

  n <- nrow(f_l)

  N <- nrow(f_u)

  p <- ncol(X_u)

  w_l <- if (is.null(w_l)) rep(1, n) else w_l / sum(w_l) * n

  w_u <- if (is.null(w_u)) rep(1, N) else w_u / sum(w_u) * N

  use_u <- is.null(lhat) || lhat != 0

  theta0 <- coef(glm(Y_l ~ . - 1,

    data = data.frame(Y_l, X_l), family = binomial))

  est <- ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u,

    opts = opts, lhat = lhat, coord = coord, w_l = w_l, w_u = w_u)

  stats <- logistic_get_stats(est, X_l, Y_l, f_l, X_u, f_u, w_l, w_u,

    use_u = use_u)

  if (is.null(lhat)) {

    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat,

      stats$grads_hat_unlabeled, stats$inv_hessian, clip = TRUE)

    return(ppi_plusplus_logistic(X_l, Y_l, f_l, X_u, f_u,

      lhat = lhat, coord = coord, opts = opts, w_l = w_l, w_u = w_u))
  }

  var_u <- cov(lhat * stats$grads_hat_unlabeled)

  var_l <- cov(stats$grads - lhat * stats$grads_hat)

  Sigma_hat <- stats$inv_hessian %*% (n/N * var_u + var_l) %*% stats$inv_hessian

  return(list(est = est, se = sqrt(diag(Sigma_hat) / n), lambda = lhat,

    rectifier_est = theta0 - est, var_u = var_u, var_l = var_l,

    grads = stats$grads, grads_hat_unlabeled = stats$grads_hat_unlabeled,

    grads_hat = stats$grads_hat, inv_hessian = stats$inv_hessian))
}

#=== END =======================================================================
