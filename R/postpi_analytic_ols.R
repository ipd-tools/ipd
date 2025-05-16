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
#' @param original (boolean): Logical argument to use original method from
#' Wang et al. (2020). Defaults to FALSE; TRUE retained for posterity.
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
#' X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])
#'
#' Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>
#'   matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>
#'   matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>
#'   matrix(ncol = 1)
#'
#' postpi_analytic_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

postpi_analytic_ols <- function(
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    original = FALSE) {

    n_l <- nrow(X_l)

    n_u <- nrow(X_u)

    #-- 1. Relationship Model

    fit_rel <- lm(Y_l ~ f_l)

    eta_hat <- residuals(fit_rel)

    gamma0_hat <- coef(fit_rel)[1]

    gamma1_hat <- coef(fit_rel)[2]

    #-- 2. Point Estimates

    if (original) {

        #- a. Naive Inference Model

        fit_inf <- lm(f_u ~ X_u - 1)

        #- b. Estimator

        est <- solve(crossprod(X_u)) %*% t(X_u) %*%

            (gamma0_hat + gamma1_hat * X_u %*% coef(fit_inf))

    } else {

        #- a. Means-Based Covariances

        C_Xf_u <- crossprod(X_u, f_u) / n_u

        C_Xe_l <- crossprod(X_l, eta_hat) / n_l

        M_XX_u <- crossprod(X_u) / n_u

        #- b. Estimator

        est <- solve(M_XX_u) %*% (gamma1_hat * C_Xf_u + C_Xe_l)
    }

    #-- 3. Standard Errors

    if (original) {

        var_hat <- solve(crossprod(X_u)) *

            (sigma(fit_rel)^2 + (gamma1_hat^2) * sigma(fit_inf)^2)

    } else {

        S1 <- crossprod(X_u * c(f_u))  / n_u  -  C_Xf_u %*% t(C_Xf_u)

        S2 <- crossprod(X_l * eta_hat) / n_l  -  C_Xe_l %*% t(C_Xe_l)

        Gamma <- gamma1_hat^2 * S1 + (n_u / n_l) * S2

        var_hat <- solve(M_XX_u) %*% Gamma %*% solve(M_XX_u) / n_u
    }

    se <- sqrt(diag(var_hat))

    #-- 4. Output

    return(list(est = as.vector(est), se = as.vector(se)))
}
