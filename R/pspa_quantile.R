#===============================================================================
# PSPA QUANTILE ESTIMATION
#===============================================================================

#--- PSPA QUANTILE ESTIMATION -----------------------------------------------

#' PSPA Quantile Estimation
#'
#' @description
#' Helper function for PSPA quantile estimation
#'
#' @details
#' Post-prediction adaptive inference
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
#' reduction. PSPA will estimate the weights if not specified.
#'
#' @param alpha (scalar): type I error rate for hypothesis testing - values in
#' (0, 1); defaults to 0.05.
#'
#' @return A list of outputs: estimate of inference model parameters and
#' corresponding standard error.
#'
#' @examples
#'
#' dat <- simdat(model = "quantile")
#'
#' form <- Y - f ~ 1
#'
#' Y_l <- dat[dat$set_label == "labeled",   all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled",   all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' pspa_quantile(Y_l, f_l, f_u, q = 0.5)
#'
#' @import stats
#'
#' @export

pspa_quantile <- function(Y_l, f_l, f_u, q,

  weights = NA, alpha = 0.05) {

  fit <- pspa_y(Y_lab = Y_l, Yhat_lab = f_l, Yhat_unlab = f_u,

    quant = q, intercept = T,

    weights = weights,

    alpha = alpha, method = "quantile")

  fit <- as.data.frame(fit)

  est <- fit$Estimate

  se <- fit$Std.Error

  return(list(est = est, se = se))
}

#=== END =======================================================================
