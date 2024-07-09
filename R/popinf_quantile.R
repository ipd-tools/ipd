#===============================================================================
# POP-INF QUANTILE ESTIMATION
#===============================================================================

#--- POP-INF QUANTILE ESTIMATION -----------------------------------------------

#' POP-Inf Quantile Estimation
#'
#' @description
#' Helper function for POP-Inf quantile estimation
#'
#' @details
#' Assumption-lean and data-adaptive post-prediction inference
#' (Miao et al. 2023) <https://arxiv.org/abs/2311.14220>
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param q (float): Quantile to estimate. Must be in the range (0, 1).
#'
#' @param weights (array): 1-dimensional array of weights vector for variance
#' reduction. POP-Inf will estimate the weights if not specified.
#'
#' @param alpha (scalar): type I error rate for hypothesis testing - values in
#' (0, 1); defaults to 0.05.
#'
#' @param delta (scalar):tolerance for assessing convergence; defaults to 0.05.
#'
#' @param K (integer): maximum number of iterations; defaults to 100.
#'
#' @returns A list of outputs: estimate of inference model parameters and
#' corresponding standard error.
#'
#' @examples
#'
#' dat <- simdat(model = "quantile")
#'
#' form <- Y - f ~ 1
#'
#' Y_l <- dat[dat$set == "labeled",   all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "labeled",   all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' popinf_quantile(Y_l, f_l, f_u, q = 0.5)
#'
#' @import stats POPInf
#'
#' @export

popinf_quantile <- function(Y_l, f_l, f_u, q,

  weights = NA, alpha = 0.05, delta = 0.05, K = 100) {

  fit <- POPInf::pop_M(Y_lab = Y_l, Yhat_lab = f_l, Yhat_unlab = f_u,

    quant = q, intercept = T, max_iterations = K,

    convergence_threshold = delta, weights = weights,

    alpha = alpha, method = "quantile")

  fit <- as.data.frame(fit)

  est <- fit$Estimate

  se <- fit$Std.Error

  return(list(est = est, se = se))
}

#=== END =======================================================================
