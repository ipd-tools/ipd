#===============================================================================
# PSPA POISSON REGRESSION
#===============================================================================

#' PSPA Poisson Regression
#'
#' @description
#' Helper function for PSPA Poisson regression
#'
#' @details
#' Post-prediction adaptive inference
#' (Miao et al. 2023) <https://arxiv.org/abs/2311.14220>
#'
#' @param X_l (matrix): n x p matrix of covariates in the labeled data.
#'
#' @param Y_l (vector): n-vector of count labeled outcomes.
#'
#' @param f_l (vector): n-vector of binary predictions in the labeled data.
#'
#' @param X_u (matrix): N x p matrix of covariates in the unlabeled data.
#'
#' @param f_u (vector): N-vector of binary predictions in the unlabeled data.
#'
#' @param weights (array): p-dimensional array of weights vector for variance
#' reduction. PSPA will estimate the weights if not specified.
#'
#' @param alpha (scalar): type I error rate for hypothesis testing - values in
#' (0, 1); defaults to 0.05
#'
#' @return A list of outputs: estimate of inference model parameters and
#' corresponding standard error.
#'
#' @examples
#'
#' dat <- simdat(model = "poisson")
#'
#' form <- Y - f ~ X1
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
#' pspa_poisson(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

pspa_poisson <- function(X_l, Y_l, f_l, X_u, f_u,

  weights = NA, alpha = 0.05) {

  fit <- pspa_y(X_lab = X_l, X_unlab = X_u,

    Y_lab = Y_l, Yhat_lab = f_l, Yhat_unlab = f_u,

    intercept = T,

    weights = weights, alpha = alpha, method = "poisson")

  fit <- as.data.frame(fit)

  est <- fit$Estimate

  se <- fit$Std.Error

  return(list(est = est, se = se))
}

#=== END =======================================================================
