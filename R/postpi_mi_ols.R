#===============================================================================
#
#  PROGRAM: postpi_mi_ols.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#           Kentaro Hoffman (khoffm3@uw.edu)
#           Awan Afiaz (aafiaz@uw.edu)
#
#
#  PURPOSE: Implementation of the Multiple Imputation algorithm of PostPI
#
#           (Insert paper name here)
#
#  INPUTS:  rel_form, inf_form, dat, m
#
#  OUTPUTS:  Post-prediction inference functions for various target estimands.
#
#  Updated: 2024-01-20
#
#===============================================================================

#=== MULTIPLE IMPUTATION =======================================================

#' PostPI Linear Regression using Salerno et al. (2024) Multiple Imputation
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
#' @param m number of imputations to be performed
#'
#' @param rel_func relationship function model form (lm, gam, glm, etc)
#'
#' @param scale_sigmar logical: scale sigma_r or not
#'
#' @returns @returns A list of outputs: estimate of inference model parameters and corresponding standard error.
#'
#' @examples
#'
#' rel_form <- Y ~ Yhat
#'
#' inf_form <- Yhat ~ X1
#'
#' dat <- simdat()
#'
#' m <- 100
#'
#' postpi_mi_ols(rel_form, inf_form, dat = dat, m)
#'
#' @export
#'
#' @import stats
#' @importFrom ranger ranger
#' @importFrom gam gam

#-- PostPI - MULTIPLE IMPUTATION

postpi_mi_ols <- function(rel_form, inf_form, dat, m = 100, rel_func = "lm",

                    scale_sigmar = F) {

  #-- 1. Estimate Prediction Model (Done in Data Step)

  #-- 2. Estimate Relationship Model

  if (rel_func == "lm") {

    fit_rel <- lm(rel_form, data = dat[dat$set == "labeled",])

  } else if (rel_func == "rf") {

    fit_rel <- ranger::ranger(rel_form, data = dat[dat$set == "labeled",], keep.inbag = T)

  } else if (rel_func == "gam") {

    fit_rel <- gam::gam(rel_form, data = dat[dat$set == "labeled",])

  } else {

    stop("Currently only 'lm', 'rf', and 'gam' are supported")
  }

  #-- 3. Multiple Imputation

  set.seed(12345)

  dat_val <- dat[dat$set == "unlabeled",]

  n_val <- nrow(dat_val)

  inf_form_m <- reformulate(all.vars(inf_form)[-1], response = "Y_tilde_m")

  if (scale_sigmar) {

    m_rel <- model.matrix(rel_form, data = dat_val)
  }

  ests_m <- sapply(1:m, function(b) {

    #-   i. Impute Outcome Values from Relationship Model

    if (rel_func == "lm") {

      if (scale_sigmar) {                                                        # Testing this out

        sigma_star <- sqrt(sum((fit_rel$residuals)^2) /

                             rchisq(1, fit_rel$df.residual))

        v_star <- solve(as.matrix(crossprod(qr.R(fit_rel$qr))))

        beta_star <- fit_rel$coefficients +

          (t(chol(v_star)) %*% rnorm(length(fit_rel$coefficients)))*sigma_star

        dat_val$Y_tilde_m <- m_rel %*% beta_star + rnorm(n_val)*sigma_star       # Testing this out

      } else {

        dat_val$Y_tilde_m <- rnorm(

          n_val, predict(fit_rel, dat_val), sigma(fit_rel))

      }

    } else if (rel_func == "rf") {

      rel_preds <- predict(fit_rel, data = dat_val, type = "se")

      dat_val$Y_tilde_m <- rnorm(n_val, rel_preds$predictions, rel_preds$se)

    } else if (rel_func == "gam") {

      dat_val$Y_tilde_m <- rnorm(

        n_val, predict(fit_rel, dat_val), sigma(fit_rel))

    } else {

      stop("Currently only 'lm', 'rf', and 'gam' are supported")
    }

    #- iii. Fit Inference Model on Simulated Outcomes

    fit_inf_m <- lm(inf_form_m, data = dat_val)

    #-  iv. Extract Coefficient Estimator

    #-   v. Extract SE of Estimator

    return(summary(fit_inf_m)$coefficients[, 1:2])
  })

  coefs_m <- ests_m[1:(nrow(ests_m)/2),] |> t()

  ses_m   <- ests_m[(nrow(ests_m)/2 + 1):nrow(ests_m),] |> t()

  #-- 4. Estimate Inference Model Coefficient

  est <- colMeans(coefs_m)

  #-- 5. Estimate Inference Model SE with Rubin's Rules

  #- a. Within-Imputation Variance

  var_w <- colMeans(ses_m^2)

  #- b. Between-Imputation Variance

  var_b <- apply(coefs_m, 1, function(Qi) { (Qi - est)^2 }) |>

    t() |> apply(2, function(param) { sum(param) / (m - 1) })

  #- c. Total Variance

  var_t <- var_w + (1 + 1/m)*var_b

  #-- Output

  return(list(est = est, se = sqrt(var_t)))
}
