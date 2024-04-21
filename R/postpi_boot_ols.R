#===============================================================================
#
#  PROGRAM: postpi_boot_ols.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#           Kentaro Hoffman (khoffm3@uw.edu)
#           Awan Afiaz (aafiaz@uw.edu)
#
#
#  PURPOSE: Implementation of bootstrap correction algorithm from Wang et al.(2020)
#
#           Methods for correcting inference based on outcomes predicted by
#           machine learning
#
#  INPUTS:  rel_form, inf_form, dat, nboot
#
#  OUTPUTS:  Post-prediction inference functions for various target estimands:
#
#  Updated: 2024-01-20
#
#===============================================================================

#=== BOOTSTRAP CORRECTION =======================================================

#' IPD Linear Regression using Wang et al. (2020) Bootstrap Correction
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
#' @param nboot number of bootstraps
#'
#' @param rel_func relationship function model form (lm, gam, glm, etc)
#'
#' @param scale_se (boolean): Logical argument to scale relationship model error variance (defaults to TRUE; retained for posterity).
#'
#' @returns A list of outputs: estimate of inference model parameters and corresponding standard error based on both parametric and non-parametric bootstrap methods.
#'
#' @examples
#'
#' dat <- simdat()
#'
#' form <- Y - Yhat ~ X1
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
#' nboot <- 100
#'
#' postpi_boot_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @export
#'
#' @import stats
#' @importFrom ranger ranger
#' @importFrom gam gam


#-- PostPI - BOOTSTRAP

postpi_boot_ols <- function(X_l, Y_l, f_l, X_u, f_u,

  nboot = 100, rel_func = "lm", scale_se = T) {

  #-- 1. Estimate Prediction Model (Done in Data Step)

  #-- 2. Estimate Relationship Model

  if (rel_func == "lm") {

    fit_rel <- lm(Y_l ~ f_l)

  } else if (rel_func == "rf") {

    fit_rel <- ranger(Y_l ~ f_l, keep.inbag = T)

  } else if (rel_func == "gam") {

    fit_rel <- gam(Y_l ~ f_l)

  } else {

    stop("Currently only 'lm', 'rf', and 'gam' are supported")
  }

  #-- 3. Bootstrap

  set.seed(12345)

  n <- nrow(X_l)

  N <- nrow(X_u)

  ests_b <- sapply(1:nboot, function(b) {

    #-   i. Sample Predicted Values and Covariates with Replacement

    idx_b <- sample(1:N, N, replace = T)

    X_u_b <- X_u[idx_b, ]

    #-  ii. Simulate Values from Relationship Model

    if (rel_func == "lm") {

      if (scale_se) {

        Y_u_b <- rnorm(N, predict(fit_rel, as.data.frame(X_u_b)),

          sigma(fit_rel) * sqrt(N / n))


      } else {

        Y_u_b <- rnorm(N, predict(fit_rel, as.data.frame(X_u_b)),

          sigma(fit_rel))
      }

    } else if (rel_func == "rf") {

      rel_preds <- predict(fit_rel, data = as.data.frame(X_u_b), type = "se")

      if (scale_se) {

        Y_u_b <- rnorm(N, rel_preds$predictions, rel_preds$se * sqrt(N / n))

      } else {

        Y_u_b <- rnorm(N, rel_preds$predictions, rel_preds$se)
      }

    } else if (rel_func == "gam") {

      if (scale_se) {

        Y_u_b <- rnorm(N, predict(fit_rel, as.data.frame(X_u_b)),

          sigma(fit_rel) * sqrt(N / n))

      } else {

        Y_u_b <- rnorm(N, predict(fit_rel, as.data.frame(X_u_b)),

          sigma(fit_rel))
      }

    } else {

      stop("Currently only 'lm', 'rf', and 'gam' are supported")
    }

    #- iii. Fit Inference Model on Simulated Outcomes

    fit_inf_b <- lm(Y_u_b ~ X_u_b - 1)

    #-  iv. Extract Coefficient Estimator

    #-   v. Extract SE of Estimator

    return(summary(fit_inf_b)$coefficients[, 1:2])
  })

  #-- 4. Estimate Inference Model Coefficient

  est <- apply(ests_b[1:(nrow(ests_b)/2),], 1, median)

  #-- 5. Estimate Inference Model SE

  #- a. Parametric Bootstrap

  se_p <- apply(ests_b[(nrow(ests_b)/2 + 1):nrow(ests_b),], 1, median)

  #- b. Nonparametric Bootstrap

  se_n <- apply(ests_b[1:(nrow(ests_b)/2),], 1, sd)

  #-- Output

  return(list(est = est, se_par = se_p, se_npar = se_n))
}
