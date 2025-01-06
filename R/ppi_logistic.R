#===============================================================================
# PPI LOGISTIC REGRESSION
#===============================================================================

#--- PPI LOGISTIC REGRESSION ---------------------------------------------------

#' PPI Logistic Regression
#'
#' @description
#' Helper function for PPI logistic regression
#'
#' @details
#' Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://www.science.org/doi/10.1126/science.adi6000>
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
#' @param opts (list, optional): Options to pass to the optimizer.
#' See ?optim for details.
#'
#' @return (list): A list containing the following:
#'
#' \describe{
#'    \item{est}{(vector): vector of PPI logistic regression coefficient
#'    estimates.}
#'    \item{se}{(vector): vector of standard errors of the coefficients.}
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
#' ppi_logistic(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_logistic <- function(X_l, Y_l, f_l, X_u, f_u, opts = NULL) {

  n <- nrow(f_l)

  N <- nrow(f_u)

  p <- ncol(X_u)

  theta0 <- coef(glm(Y_l ~ . - 1,

    data = data.frame(Y_l, X_l), family = binomial))

  est <- ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u,

    opts = opts, lhat = 1)

  stats <- logistic_get_stats(est, X_l, Y_l, f_l, X_u, f_u, use_u = T)

  var_u <- cov(stats$grads_hat_unlabeled)

  var_l <- cov(stats$grads - stats$grads_hat)

  Sigma_hat <- stats$inv_hessian %*% (n/N * var_u + var_l) %*% stats$inv_hessian

  return(list(est = est, se = sqrt(diag(Sigma_hat) / n),

    rectifier_est = theta0 - est, var_u = var_u, var_l = var_l,

    grads = stats$grads, grads_hat_unlabeled = stats$grads_hat_unlabeled,

    grads_hat = stats$grads_hat, inv_hessian = stats$inv_hessian))
}

#=== END =======================================================================
