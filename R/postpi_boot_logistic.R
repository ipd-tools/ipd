#===============================================================================
#  POSTPI BOOTSTRAP LOGISTIC REGRESSION
#===============================================================================

#--- POSTPI BOOTSTRAP LOGISTIC REGRESSION --------------------------------------

#' PostPI Logistic Regression (Bootstrap Correction)
#'
#' @description
#' Helper function for PostPI logistic regression (bootstrap correction)
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
#' @param nboot (integer): Number of bootstrap samples. Defaults to 100.
#'
#' @param se_type (string): Which method to calculate the standard errors.
#' Options include "par" (parametric) or "npar" (nonparametric).
#' Defaults to "par".
#'
#' @param seed (optional) An \code{integer} seed for random number generation.
#'
#' @return A list of outputs: estimate of inference model parameters and
#' corresponding standard error based on both parametric and non-parametric
#' bootstrap methods.
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
#' postpi_boot_logistic(X_l, Y_l, f_l, X_u, f_u, nboot = 200)
#'
#' @import stats
#'
#' @import caret
#'
#' @export

postpi_boot_logistic <- function(X_l, Y_l, f_l, X_u, f_u,

  nboot = 100, se_type = "par", seed = NULL) {

  #-- 1. Estimate Prediction Model (Done in Data Step)

  #-- 2. Estimate Relationship Model

  fit_rel <- caret::train(

    factor(Y) ~ factor(f), data = data.frame(Y = Y_l, f = f_l),

    method = "knn")

  #-- 3. Bootstrap

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(X_l)

  N <- nrow(X_u)

  ests_b <- sapply(1:nboot, function(b) {

    #-   i. Sample Predicted Values and Covariates with Replacement

    idx_b <- sample(1:N, N, replace = T)

    f_u_b <- f_u[idx_b, ]

    X_u_b <- X_u[idx_b, ]

    #-  ii. Simulate Values from Relationship Model

    prob_b <- predict(fit_rel, data.frame(f = f_u_b), type = "prob")[,2]

    Y_u_b  <- rbinom(N, 1, prob_b)

    #- iii. Fit Inference Model on Simulated Outcomes

    fit_inf_b <- glm(Y_u_b ~ X_u_b - 1, family = binomial)

    #-  iv. Extract Coefficient Estimator

    #-   v. Extract SE of Estimator

    return(summary(fit_inf_b)$coefficients[, 1:2])
  })

  #-- 4. Estimate Inference Model Coefficient

  est <- apply(ests_b[1:(nrow(ests_b)/2),], 1, median)

  #-- 5. Estimate Inference Model SE

  if (se_type == "par") {

    #- a. Parametric Bootstrap

    se <- apply(ests_b[(nrow(ests_b)/2 + 1):nrow(ests_b),], 1, median)

  } else if (se_type == "npar") {

    #- b. Nonparametric Bootstrap

    se <- apply(ests_b[1:(nrow(ests_b)/2),], 1, sd)

  } else {

    stop("se_type must be one of 'par' or 'npar'")
  }

  #-- Output

  return(list(est = est, se = se))
}
