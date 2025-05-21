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
#' @param intercept (Logical): Do the design matrices include intercept
#' columns? Default is \code{TRUE}.
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
#' chen_ols(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE)
#'
#' @import stats
#'
#' @export

chen_ols <- function(
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    intercept = TRUE) {

    f_a <- rbind(f_l, f_u)
    X_a <- rbind(X_l, X_u)

    n_l <- nrow(Y_l)
    n_a <- nrow(f_a)

    if (intercept){

        fit_beta_l  <- lm(Y_l ~ X_l - 1)
        fit_gamma_l <- lm(f_l ~ X_l - 1)
        fit_gamma_a <- lm(f_a ~ X_a - 1)

    } else {

        fit_beta_l  <- lm(Y_l ~ X_l)
        fit_gamma_l <- lm(f_l ~ X_l)
        fit_gamma_a <- lm(f_a ~ X_a)
    }

    coef_beta_l  <- coef(fit_beta_l)
    coef_gamma_l <- coef(fit_gamma_l)
    coef_gamma_a <- coef(fit_gamma_a)

    D1 <- -crossprod(X_l, X_l) / n_l
    D2 <- -crossprod(X_l, X_l) / n_l

    score_beta  <- X_l * residuals(fit_beta_l)
    score_gamma <- X_l * residuals(fit_gamma_l)

    C11 <- crossprod(score_beta,  score_beta)  / n_l
    C12 <- crossprod(score_beta,  score_gamma) / n_l
    C22 <- crossprod(score_gamma, score_gamma) / n_l

    D1_inv  <- solve(D1)
    C22_inv <- solve(C22)

    D1_C12     <- D1_inv %*% C12
    D1_C12_C22 <- D1_C12 %*% C22_inv

    theta_hat <- coef_beta_l -
        D1_C12_C22 %*% D2 %*% (coef_gamma_l - coef_gamma_a)

    Omega <- (D1_inv %*% C11 %*% D1_inv -
        (1 - n_l / n_a) * D1_C12_C22 %*% t(D1_C12)) / n_l

    se_beta <- sqrt(diag(Omega))

    #-- Return Results

    return(list(est = theta_hat, se = se_beta))
}
