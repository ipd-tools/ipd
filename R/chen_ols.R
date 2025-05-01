#===============================================================================
# CHEN & CHEN ORDINARY LEAST SQUARES
#===============================================================================

#--- CHEN & CHEN OLS -----------------------------------------------------------

#' Chen & Chen OLS
#'
#' @description
#' Helper function for Chen & Chen OLS estimation
#'
#' @details
#' Another look at inference after prediction (Gronsbell et al., 2025)
#' <https://arxiv.org/pdf/2411.19908>
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
#' @return (list): A list containing the following:
#'
#' \describe{
#'    \item{est}{(vector): vector of Chen & Chen OLS regression coefficient
#'    estimates.}
#'    \item{se}{(vector): vector of standard errors of the coefficients.}
#' }
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
#' #chen_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

chen_ols <- function(
    Y_l,
    X_l,
    f_l,
    f_u,
    X_u) {

    f_a <- c(f_l, f_u)

    X_a <- rbind(X_l, X_u)

    #-- Sample Sizes

    n_l <- length(Y_l)

    n_total <- length(f_a)

    #-- Fit Models

    #- Model for Observed Outcomes

    beta_model <- lm(Y_l ~ ., data = cbind(Y_l, X_l))

    #- Model for Predictions in Labeled Set

    gamma_labeled_model <- lm(f_l ~ ., data = cbind(f_l, X_l))

    #- Model for Predictions in All Data

    gamma_all_model <- lm(f_a ~ ., data = cbind(f_a, X_a))

    #-- Extract Coefficients

    #- Coefficients for Observed Outcomes

    beta_coeff <- coef(beta_model)

    #- Coefficients for Predictions in Labeled Data

    gamma_labeled_coeff <- coef(gamma_labeled_model)

    #- Coefficients for Predictions in All Data

    gamma_all_coeff <- coef(gamma_all_model)


    ########### HOW TO HANDLE THIS WITH CURRENT WRAPPER FUNCTION? ##############

    # Add intercept for labeled covariates
    X_l_int <- cbind(1, as.matrix(X_l))

    #-- Compute Derivative Matrices

    D1 <- -crossprod(X_l_int, X_l_int) / n_l

    D2 <- -crossprod(X_l_int, X_l_int) / n_l

    #-- Compute Score Matrices

    score_beta <- X_l_int * residuals(beta_model)

    score_gamma <- X_l_int * residuals(gamma_labeled_model)

    #-- Covariance Matrices

    C11 <- crossprod(score_beta, score_beta) / n_l

    C12 <- crossprod(score_beta, score_gamma) / n_l

    C22 <- crossprod(score_gamma, score_gamma) / n_l

    #-- Matrix Inverses

    D1_inv <- solve(D1)

    C22_inv <- solve(C22)

    #-- Compute Correction Terms

    D1_C12 <- D1_inv %*% C12

    D1_C12_C22 <- D1_C12 %*% C22_inv

    #-- Calculate Adjusted Coefficient Estimates

    theta_hat <- beta_coeff -

        D1_C12_C22 %*% D2 %*% (gamma_labeled_coeff - gamma_all_coeff)

    #-- Compute Variance-Covariance Matrix for the Adjusted Estimates

    Omega <- (D1_inv %*% C11 %*% D1_inv -

        (1 - n_l / n_total) * D1_C12_C22 %*% t(D1_C12)) / n_l

    #-- Standard Errors for the Adjusted Estimates

    se_beta <- sqrt(diag(Omega))

    #-- Return Results

    return(c(est = theta_hat, se = se_beta))
}

#=== END =======================================================================
