#===============================================================================
#
#  PROGRAM: postpi_analytic_ols.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#           Kentaro Hoffman (khoffm3@uw.edu)
#           Awan Afiaz (aafiaz@uw.edu)
#
#
#  PURPOSE: Implementation of analytic correction algorithm from Wang et al.(2020)
#
#           Methods for correcting inference based on outcomes predicted by
#           machine learning
#
#  INPUTS:  rel_form, inf_form, dat
#
#  OUTPUTS:  Post-prediction inference functions for various target estimands:
#
#  Updated: 2024-01-20
#
#===============================================================================

#=== ANALYTIC CORRECTION =======================================================

#' PostPI Linear Regression using Wang et al. (2020) Analytic Correction
#'
#' @description
#' A short description.
#'
#' @details
#' Additional details...
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
#' @param X_in_rel (boolean, optional): Logical argument for whether or not to
#' include the features of interest in the relationship model.
#' Defaults to FALSE.
#'
#' @returns A list of outputs: estimate of the inference model parameters and
#' corresponding standard error estimate.
#'
#'
#' @examples
#'
#' rel_form <- Y ~ Yhat
#'
#' inf_form <- Yhat ~ X1
#'
#' dat <- simdat()
#'
#' postpi_analytic_ols(rel_form, inf_form, dat = dat)
#'
#' @export
#'
#' @import stats
#'

#-- PostPI - ANALYTIC for OLS

postpi_analytic_ols <- function(X_l, Y_l, f_l, X_u, f_u, X_in_rel = FALSE) {

  #- 0. Setup

  inf_form <- reformulate(

    all.vars(formula)[-(1:2)], response = all.vars(formula)[2])

  if(X_in_rel) {

    rel_form <- reformulate(

      all.vars(formula)[-1], response = all.vars(formula)[1])

  } else {

    rel_form <- reformulate(

      all.vars(formula)[2], response = all.vars(formula)[1])
  }

  #- 1. Estimate Relationship Model

  fit_rel <- lm(rel_form, data = as.data.frame(cbind(Y_l, f_l, X_l)))

  #- 2. Estimate Inference Model

  fit_inf <- lm(inf_form, data = as.data.frame(cbind(f_u, X_u)))

  #- 3. Coefficient Estimator

  X_val <- model.matrix(inf_form, data = as.data.frame(cbind(f_u, X_u)))

  est <- solve(crossprod(X_val)) %*% t(X_val) %*%

    (coef(fit_rel)[1] + coef(fit_rel)[2]*X_val %*% coef(fit_inf))

  #- 4. SE of Coefficient Estimator

  se <- sqrt(diag(solve(crossprod(X_val))*(sigma(fit_rel)^2 +

    (coef(fit_rel)[2]^2)*sigma(fit_inf)^2)))

  #- Output

  return(list(est = as.vector(est), se = as.vector(se)))
}

