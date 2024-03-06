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
#' @description
#' A short description...
#'
#' @details
#' Additional details...
#'
#' @param rel_form A formula defining the relationship model. This should be
#' of the form Y ~ Yhat, where Y is the name of the column corresponding to
#' the observed outcome in the labeled data and Yhat is the name of the column
#' corresponding to the predicted outcome in the labeled data.
#'
#' @param inf_form A formula defining the inference model. This should be of
#' the form Yhat ~ X, where Yhat is the name of the column corresponding to the
#' predicted outcome in the unlabeled data, and X generally corresponds to the
#' features of interest (e.g., X1 + X2).
#'
#' @param dat data in the form of the simdat function
#'
#' @param nboot number of bootstraps
#'
#' @param rel_func relationship function model form (lm, gam, glm, etc)
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
#' postpi_boot_ols(X_l, Y_l, f_l, X_u, f_u, nboot)
#'
#' @export
#'
#' @import stats
#' @importFrom ranger ranger
#' @importFrom gam gam


#-- PostPI - BOOTSTRAP

postpi_boot_ols <- function(X_l, Y_l, f_l, X_u, f_u,

  nboot = 100, rel_func = "lm") {

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

  n_u <- nrow(X_u)

  ests_b <- sapply(1:nboot, function(b) {

    #-   i. Sample Predicted Values and Covariates with Replacement

    idx_b <- sample(1:n_u, n_u, replace = T)

    X_u_b <- X_u_b[idx_b, ]

    #-  ii. Simulate Values from Relationship Model

    if (rel_func == "lm") {

      Y_u_b <- rnorm(n_u, predict(fit_rel, X_u_b), sigma(fit_rel))

    } else if (rel_func == "rf") {

      rel_preds <- predict(fit_rel, data = X_u_b, type = "se")

      Y_u_b <- rnorm(n_u, rel_preds$predictions, rel_preds$se)

    } else if (rel_func == "gam") {

      Y_u_b <- rnorm(n_u, predict(fit_rel, X_u_b), sigma(fit_rel))

    } else {

      stop("Currently only 'lm', 'rf', and 'gam' are supported")
    }

    #- iii. Fit Inference Model on Simulated Outcomes

    fit_inf_b <- lm(Y_u_b ~ X_u_b)

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
