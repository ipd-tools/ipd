#===============================================================================
# POSTPI ANALYTIC ORDINARY LEAST SQUARES
#===============================================================================

#--- POSTPI ANALYTIC OLS -------------------------------------------------------

#' PostPI OLS (Analytic Correction)
#'
#' @description
#' Helper function for PostPI OLS estimation (analytic correction)
#'
#' @details
#' Methods for correcting inference based on outcomes predicted by machine
#' learning (Wang et al., 2020)
#' <https://www.pnas.org/doi/abs/10.1073/pnas.2001238117>
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
#' @param scale_se (boolean): Logical argument to scale relationship model
#' error variance. Defaults to TRUE; FALSE option is retained for posterity.
#'
#' @param n_t (integer, optional) Size of the dataset used to train the
#' prediction function (necessary if \code{n_t} < \code{nrow(X_l)}.
#' Defaults to \code{Inf}.
#'
#' @return A list of outputs: estimate of the inference model parameters and
#' corresponding standard error estimate.
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
#' postpi_analytic_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

postpi_analytic_ols <- function(X_l, Y_l, f_l, X_u, f_u,

  scale_se = TRUE, n_t = Inf) {

  #- 1. Estimate Relationship Model

  fit_rel <- lm(Y_l ~ f_l)

  #- 2. Estimate Inference Model

  fit_inf <- lm(f_u ~ X_u - 1)

  #- 3. Coefficient Estimator

  est <- solve(crossprod(X_u)) %*% t(X_u) %*%

    (coef(fit_rel)[1] + coef(fit_rel)[2] * X_u %*% coef(fit_inf))

  #- 4. SE of Coefficient Estimator

  if (scale_se) {

    se <- sqrt(diag(solve(crossprod(X_u)) *

      (sigma(fit_rel)^2 * nrow(X_u) / min(nrow(X_l), n_t) +

      (coef(fit_rel)[2]^2) * sigma(fit_inf)^2)))

  } else {

    se <- sqrt(diag(solve(crossprod(X_u)) *

      (sigma(fit_rel)^2 + (coef(fit_rel)[2]^2) * sigma(fit_inf)^2)))
  }

  #- Output

  return(list(est = as.vector(est), se = as.vector(se)))
}

#=== END =======================================================================
