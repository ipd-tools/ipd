#===============================================================================
#  POSTPI BOOTSTRAP ORDINARY LEAST SQUARES
#===============================================================================

#--- POSTPI BOOTSTRAP OLS ------------------------------------------------------

#' PostPI OLS (Bootstrap Correction)
#'
#' @description
#' Helper function for PostPI OLS estimation (bootstrap correction)
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
#' @param rel_func (string): Method for fitting the relationship model.
#' Options include "lm" (linear model), "rf" (random forest), and "gam"
#' (generalized additive model). Defaults to "lm".
#'
#' @param scale_se (boolean): Logical argument to scale relationship model
#' error variance. Defaults to TRUE; FALSE option is retained for posterity.
#'
#' @param n_t (integer, optional) Size of the dataset used to train the
#' prediction function (necessary if \code{n_t} < \code{nrow(X_l)}.
#' Defaults to \code{Inf}.
#'
#' @param seed (optional) An \code{integer} seed for random number generation.
#'
#' @return A list of outputs: estimate of inference model parameters and
#' corresponding standard error based on both parametric and non-parametric
#' bootstrap methods.
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
#' postpi_boot_ols(X_l, Y_l, f_l, X_u, f_u, nboot = 200)
#'
#' @import stats
#'
#' @importFrom ranger ranger
#'
#' @importFrom gam gam
#'
#' @importFrom splines ns
#'
#' @export

postpi_boot_ols <- function(X_l, Y_l, f_l, X_u, f_u,

  nboot = 100, se_type = "par", rel_func = "lm", scale_se = TRUE, n_t = Inf,

  seed = NULL) {

  #-- 1. Estimate Prediction Model (Done in Data Step)

  #-- 2. Estimate Relationship Model

  if (rel_func == "lm") {

    fit_rel <- lm(Y ~ f, data = data.frame(Y = Y_l, f = f_l))

  } else if (rel_func == "rf") {

    fit_rel <- ranger(Y ~ f, data = data.frame(Y = Y_l, f = f_l),

      keep.inbag = T)

  } else if (rel_func == "gam") {

    fit_rel <- lm(Y ~ splines::ns(f, df = 3),

      data = data.frame(Y = Y_l, f = f_l))

  } else {

    stop("Currently only 'lm', 'rf', and 'gam' are supported")
  }

  #-- 3. Bootstrap

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(X_l)

  N <- nrow(X_u)

  est_b <- se_b <- matrix(nrow = nboot, ncol = ncol(X_u))

  for(b in 1:nboot) {

    #-   i. Sample Predicted Values and Covariates with Replacement

    idx_b <- sample(1:N, N, replace = T)

    f_u_b <- f_u[idx_b, ]

    X_u_b <- X_u[idx_b, ]

    #-  ii. Simulate Values from Relationship Model

    if (rel_func == "lm") {

      if (scale_se) {

        Y_u_b <- rnorm(N, predict(fit_rel, data.frame(f = f_u_b)),

          sigma(fit_rel) * sqrt(N / min(n, n_t)))

      } else {

        Y_u_b <- rnorm(N, predict(fit_rel, data.frame(f = f_u_b)),

          sigma(fit_rel))
      }

    } else if (rel_func == "rf") {

      rel_preds <- predict(fit_rel, data = data.frame(f = f_u_b), type = "se")

      if (scale_se) {

        Y_u_b <- rnorm(N, rel_preds$predictions,

          rel_preds$se * sqrt(N / min(n, n_t)))

      } else {

        Y_u_b <- rnorm(N, rel_preds$predictions, rel_preds$se)
      }

    } else if (rel_func == "gam") {

      if (scale_se) {

        Y_u_b <- rnorm(N, predict(fit_rel, data.frame(f = f_u_b)),

          sigma(fit_rel) * sqrt(N / min(n, n_t)))

      } else {

        Y_u_b <- rnorm(N, predict(fit_rel, data.frame(f = f_u_b)),

          sigma(fit_rel))
      }

    } else {

      stop("Currently only 'lm', 'rf', and 'gam' are supported")
    }

    #- iii. Fit Inference Model on Simulated Outcomes

    fit_inf_b <- lm(Y_u_b ~ X_u_b - 1)

    #-  iv. Extract Coefficient Estimator

    est_b[b,] <- summary(fit_inf_b)$coefficients[, 1]

    #-   v. Extract SE of Estimator

    se_b[b,]  <- summary(fit_inf_b)$coefficients[, 2]
  }

  #-- 4. Estimate Inference Model Coefficient

  est <- apply(est_b, 2, mean)

  #-- 5. Estimate Inference Model SE

  if (se_type == "par") {

    #- a. Parametric Bootstrap

    se <- apply(se_b, 2, mean)

  } else if (se_type == "npar") {

    #- b. Nonparametric Bootstrap

    se <- apply(est_b, 2, sd)

  } else {

    stop("se_type must be one of 'par' or 'npar'")
  }

  #-- Output

  return(list(est = est, se = se))
}
