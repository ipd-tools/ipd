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
#' postpi_boot_logistic(X_l, Y_l, f_l, X_u, f_u, nboot = 200)
#'
#' @import stats
#'
#' @import caret
#'
#' @export

postpi_boot_logistic <- function(
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    nboot = 100,
    se_type = "par") {

    fit_rel <- caret::train(factor(Y) ~ factor(f),

        data = data.frame(Y = Y_l, f = f_l), method = "knn")

    N <- nrow(X_u)
    p <- ncol(X_u)

    boot_mat <- matrix(NA_real_, nrow = 2*p, ncol = nboot)

    for (b in seq_len(nboot)) {

        idx   <- sample.int(N, N, replace = TRUE)
        f_u_b <- f_u[idx]
        X_u_b <- X_u[idx, , drop = FALSE]

        prob_b <- predict(fit_rel,

            newdata = data.frame(f = as.numeric(f_u_b)), type = "prob")[,2]

        Y_u_b <- rbinom(N, 1, prob_b)

        fit_inf <- glm(Y_u_b ~ X_u_b - 1, family = binomial)

        cm <- summary(fit_inf)$coefficients[, seq_len(2)]

        boot_mat[, b] <- c(cm[,1], cm[,2])
    }

    est_vec <- apply(boot_mat[seq_len(p), ], 1, median)

    if (se_type == "par") {

        se_vec <- apply(boot_mat[(p+1):(2*p), ], 1, median)

    } else if (se_type == "npar") {

        se_vec <- apply(boot_mat[      seq_len(p), ], 1, sd)

    } else {

        stop("`se_type` must be either 'par' or 'npar'")
    }

    list(est = est_vec, se = se_vec)
}
