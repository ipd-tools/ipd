#===============================================================================
# PPI++ LOGISTIC REGRESSION
#===============================================================================

#=== PPI++ LOGISTIC REGRESSION POINT ESTIMATE ==================================

#' PPI++ Logistic Regression Point Estimate
#'
#' @description
#' Computes the prediction-powered point estimate of the logistic regression
#' coefficients.
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
#' @param opts (list, optional): Options to pass to the optimizer.
#' See ?optim for details.
#'
#' @param w_l (vector, optional): n-vector of sample weights for the
#' labeled data.
#'
#' @param w_u (vector, optional): N-vector of sample weights for the
#' unlabeled data.
#'
#' @returns (vector): p-vector of prediction-powered point estimates of the
#' logistic regression coefficients.
#'
#' @examples
#'
#' dat <- simdat_logistic()
#'
#' form <- Y - Yhat ~ Xc
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

  rectified_logistic_loss <- function(theta) {

    sum(w_u * (-f_u * (X_u %*% theta) + log1pexp(X_u %*% theta))) *

      lhat_curr / N - sum(

        w_l * (-f_l * (X_l %*% theta) + log1pexp(X_l %*% theta))) *

      lhat_curr / n + sum(

        w_l * (-Y_l * (X_l %*% theta) + log1pexp(X_l %*% theta))) / n
  }

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

#=== PPI++ LOGISTIC GRADIENT AND HESSIAN =======================================

#' Logistic Regression Gradient and Hessian
#'
#' @description
#' Computes the statistics needed for the logstic regression-based
#' prediction-powered inference.
#'
#' @param est (vector): Point estimates of the coefficients.
#'
#' @param X_l (matrix): Covariates for the labeled data set.
#'
#' @param Y_l (vector): Labels for the labeled data set.
#'
#' @param f_l (vector): Predictions for the labeled data set.
#'
#' @param X_u (matrix): Covariates for the unlabeled data set.
#'
#' @param f_u (vector): Predictions for the unlabeled data set.
#'
#' @param w_l (vector, optional): Sample weights for the labeled data set.
#'
#' @param w_u (vector, optional): Sample weights for the unlabeled data set.
#'
#' @param use_u (bool, optional): Whether to use the unlabeled data set.
#'
#' @returns (list): A list containing the following:
#'
#' \describe{
#'    \item{grads}{(matrix): n x p matrix gradient of the loss function with
#'    respect to the coefficients.}
#'    \item{grads_hat}{(matrix): n x p matrix gradient of the loss function
#'    with respect to the coefficients, evaluated using the labeled
#'    predictions.}
#'    \item{grads_hat_unlabeled}{(matrix): N x p matrix gradient of the loss
#'    function with respect to the coefficients, evaluated using the unlabeled
#'    predictions.}
#'    \item{inv_hessian}{(matrix): p x p matrix inverse Hessian of the loss
#'    function with respect to the coefficients.}
#' }
#'
#' @examples
#'
#' dat <- simdat_logistic()
#'
#' form <- Y - Yhat ~ Xc
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
#' est <- ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u)
#'
#' stats <- logistic_get_stats(est, X_l, Y_l, f_l, X_u, f_u)
#'
#' @export

logistic_get_stats <- function(est, X_l, Y_l, f_l, X_u, f_u,

  w_l = NULL, w_u = NULL, use_u = TRUE) {

  n <- nrow(f_l)

  N <- nrow(f_u)

  p <- ncol(X_u)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  mu_l <- plogis(X_l %*% est)

  mu_u <- plogis(X_u %*% est)

  hessian <- matrix(0, nrow = p, ncol = p)

  grads_hat_unlabeled <- matrix(0, nrow = N, ncol = p)

  if (use_u) {

    for (i in 1:N) {

      hessian <- hessian + w_u[i] / (N + n) * mu_u[i] * (1 - mu_u[i]) *

        tcrossprod(X_u[i, ])

      grads_hat_unlabeled[i, ] <- w_u[i] * X_u[i, ] * (mu_u[i] - f_u[i])
    }
  }

  grads <- matrix(0, nrow = n, ncol = p)

  grads_hat <- matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {

    if (use_u) {

      hessian <- hessian + w_l[i] / (N + n) * mu_l[i] * (1 - mu_l[i]) *

        tcrossprod(X_l[i, ])

    } else {

      hessian <- hessian + w_l[i] / n * mu_l[i] * (1 - mu_l[i]) *

        tcrossprod(X_l[i, ])
    }

    grads[i, ] <- w_l[i] * X_l[i, ] * (mu_l[i] - Y_l[i])

    grads_hat[i, ] <- w_l[i] * X_l[i, ] * (mu_l[i] - f_l[i])
  }

  inv_hessian <- solve(hessian)

  return(

    list(grads = grads, grads_hat = grads_hat,

      grads_hat_unlabeled = grads_hat_unlabeled, inv_hessian = inv_hessian))
}

#=== PPI++ LOGISTIC REGRESSION =================================================

#' PPI++ Logistic Regression Estimator and Inference
#'
#' @description
#' Computes the prediction-powered confidence interval for the logistic
#' regression coefficients using the PPI++ algorithm.
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
#' @param opts (list, optional): Options to pass to the optimizer.
#' See ?optim for details.
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
#'    \item{est}{(vector): p-vector of PPI++ logistic regression coefficient
#'    estimates.}
#'    \item{se}{(vector): p-vector of standard errors of the coefficients.}
#' }
#'
#' @examples
#'
#' dat <- simdat_logistic()
#'
#' form <- Y - Yhat ~ Xc
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

  return(list(est = est, se = sqrt(diag(Sigma_hat) / n)))
}

#=== END =======================================================================
