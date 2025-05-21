#--- PSPA M-ESTIMATION FOR PREDICTED -------------------------------------------

#' PSPA M-Estimation for ML-predicted labels
#'
#' \code{pspa_y} function conducts post-prediction M-Estimation.
#'
#' @param X_l Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_u Array or data.frame containing observed or predicted
#' covariates in unlabeled data.
#'
#' @param Y_l Array or data.frame of observed outcomes in
#' labeled data.
#'
#' @param f_l Array or data.frame of predicted outcomes in
#' labeled data.
#'
#' @param f_u Array or data.frame of predicted outcomes in
#' unlabeled data.
#'
#' @param alpha Specifies the confidence level as 1 - alpha for
#' confidence intervals.
#'
#' @param weights weights vector PSPA linear regression (d-dimensional, where
#' d equals the number of covariates).
#'
#' @param quant quantile for quantile estimation
#'
#' @param intercept Boolean indicating if the input covariates' data contains
#' the intercept (TRUE if the input data contains)
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return  A summary table presenting point estimates, standard error,
#' confidence intervals (1 - alpha), P-values, and weights.

pspa_y <- function(
    X_l = NA,
    X_u = NA,
    Y_l,
    f_l,
    f_u,
    alpha = 0.05,
    weights = NA,
    quant = NA,
    intercept = FALSE,
    method) {

    if (method %in% c("ols", "logistic", "poisson")) {

        if (intercept) {

            X_l <- as.matrix(X_l)
            X_u <- as.matrix(X_u)

        } else {

            X_l <- cbind(1, as.matrix(X_l))
            X_u <- cbind(1, as.matrix(X_u))
        }
    }

    Y_l <- as.matrix(as.numeric(unlist(Y_l)))

    f_l <- as.matrix(as.numeric(unlist(f_l)))

    f_u <- as.matrix(as.numeric(unlist(f_u)))

    n <- nrow(Y_l)

    N <- nrow(f_u)

    q <- if (method %in% c("mean", "quantile")) 1 else ncol(X_l)

    est <- est_ini(X_l, Y_l, quant, method)

    if (is.na(sum(weights))) {

        current_w <- rep(0, q)

        optimized_weight <- optim_weights(X_l, X_u, Y_l, f_l, f_u,

            current_w, est, quant, method)

        current_w <- optimized_weight

    } else {

        current_w <- weights
    }

    est <- optim_est(X_l, X_u, Y_l, f_l, f_u, current_w, est, quant, method)

    final_sigma <- Sigma_cal(X_l, X_u, Y_l, f_l, f_u,

        current_w, est, quant, method)

    standard_errors <- sqrt(diag(final_sigma) / n)

    p_values <- 2 * pnorm(abs(est / standard_errors), lower.tail = FALSE)

    lower_ci <- est - qnorm(1 - alpha / 2) * standard_errors

    upper_ci <- est + qnorm(1 - alpha / 2) * standard_errors

    output_table <- data.frame(

        Estimate  = est,
        Std.Error = standard_errors,
        Lower.CI  = lower_ci,
        Upper.CI  = upper_ci,
        P.value   = p_values,
        Weight    = current_w
    )

    colnames(output_table) <- c("Estimate", "Std.Error", "Lower.CI",

        "Upper.CI", "P.value", "Weight")

    return(output_table)
}

#--- BREAD FOR PSPA ------------------------------------------------------------

#' Calculation of the matrix A based on single dataset
#'
#' \code{A} function for the calculation of the matrix A based on single dataset
#'
#' @param X Array or data.frame containing covariates
#'
#' @param Y Array or data.frame of outcomes
#'
#' @param quant quantile for quantile estimation
#'
#' @param theta parameter theta
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return  matrix A based on single dataset

A <- function(
    X,
    Y,
    quant = NA,
    theta,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    p <- ncol(X)
    n <- nrow(X)

    if (method == "ols") {

        A <- crossprod(X) / n

    } else if (method == "quantile") {

        Y_vec <- if (is.list(Y)) unlist(Y) else Y

        A <- vapply(theta,

            function(a) density(Y_vec, from = a, to = a, n = 1)$y,

            FUN.VALUE = numeric(1))

    } else if (method == "mean") {

        A <- 1

    } else if (method %in% c("logistic", "poisson")) {

        mid <- sqrt(diag(as.vector(link_Hessian(X %*% theta, method)))) %*% X

        A <- crossprod(mid) / n
    }

    A
}

#--- INITIAL ESTIMATES ---------------------------------------------------------

#' Initial estimation
#'
#' \code{est_ini} function for initial estimation
#'
#' @param X Array or data.frame containing covariates
#'
#' @param Y Array or data.frame of outcomes
#'
#' @param quant quantile for quantile estimation
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return initial estimator

est_ini <- function(
    X,
    Y,
    quant = NA,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    if (method == "ols") {

        est <- lm(Y ~ X - 1)$coef

    } else if (method == "quantile") {

        est <- quantile(Y, quant)

    } else if (method == "mean") {

        est <- mean(Y)

    } else if (method == "logistic") {

        est <- glm(Y ~ X - 1, family = binomial)$coef

    } else if (method == "poisson") {

        est <- glm(Y ~ X - 1, family = poisson)$coef
    }

    return(est)
}

#--- ESTIMATING EQUATION -------------------------------------------------------

#' Estimating equation
#'
#' \code{psi} function for estimating equation
#'
#' @param X Array or data.frame containing covariates
#'
#' @param Y Array or data.frame of outcomes
#'
#' @param theta parameter theta
#'
#' @param quant quantile for quantile estimation
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return estimating equation

psi <- function(
    X,
    Y,
    theta,
    quant = NA,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    if (method == "quantile") {

        psi <- t(as.matrix(-quant + 1 * (as.numeric(Y) <= as.vector(theta))))

    } else if (method == "mean") {

        psi <- t(as.matrix(as.vector(theta) - as.numeric(Y)))

    } else if (method %in% c("ols", "logistic", "poisson")) {

        t <- X %*% theta

        res <- Y - link_grad(t, method)

        psi <- -t(as.vector(res) * X)
    }

    return(psi)
}

#--- SAMPLE EXPECTATION OF PSI -------------------------------------------------

#' Sample expectation of psi
#'
#' \code{mean_psi} function for sample expectation of psi
#'
#' @param X Array or data.frame containing covariates
#'
#' @param Y Array or data.frame of outcomes
#'
#' @param theta parameter theta
#'
#' @param quant quantile for quantile estimation
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return sample expectation of psi

mean_psi <- function(
    X,
    Y,
    theta,
    quant = NA,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    psi <- as.matrix(rowMeans(psi(X, Y, theta, quant, method)))

    return(psi)
}

#--- SAMPLE EXPECTATION OF PSPA PSI --------------------------------------------

#' Sample expectation of PSPA psi
#'
#' \code{mean_psi_pop} function for sample expectation of PSPA psi
#'
#' @param X_l Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_u Array or data.frame containing observed or predicted
#' covariates in unlabeled data.
#'
#' @param Y_l Array or data.frame of observed outcomes in labeled data.
#'
#' @param f_l Array or data.frame of predicted outcomes in labeled data.
#'
#' @param f_u Array or data.frame of predicted outcomes in
#' unlabeled data.
#'
#' @param w weights vector PSPA linear regression (d-dimensional, where
#' d equals the number of covariates).
#'
#' @param theta parameter theta
#'
#' @param quant quantile for quantile estimation
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return sample expectation of PSPA psi

mean_psi_pop <- function(
    X_l,
    X_u,
    Y_l,
    f_l,
    f_u,
    w,
    theta,
    quant = NA,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    if (method %in% c("ols", "logistic", "poisson")) {

        psi_pop <- mean_psi(X_l, Y_l, theta, quant, method) + diag(w) %*%

            (mean_psi(X_u, f_u, theta, quant, method) -

                mean_psi(X_l, f_l, theta, quant, method))

    } else if (method %in% c("mean", "quantile")) {

        psi_pop <- mean_psi(X_l, Y_l, theta, quant, method) + w *

            (mean_psi(X_u, f_u, theta, quant, method) -

                mean_psi(X_l, f_l, theta, quant, method))
    }

    return(psi_pop)
}

#=== STATISTICS FOR GLMS =======================================================

#--- GRADIENT OF THE LINK FUNCTION ---------------------------------------------

#' Gradient of the link function
#'
#' \code{link_grad} function for gradient of the link function
#'
#' @param t t
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return gradient of the link function

link_grad <- function(
    t,
    method = c("ols", "logistic", "poisson")) {

    method <- match.arg(method)

    if (method == "ols") {

        grad <- t

    } else if (method == "logistic") {

        grad <- 1 / (1 + exp(-t))

    } else if (method == "poisson") {

        grad <- exp(t)
    }

    return(grad)
}

#--- HESSIAN OF THE LINK FUNCTION ----------------------------------------------

#' Hessians of the link function
#'
#' \code{link_Hessian} function for Hessians of the link function
#'
#' @param t t
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return Hessians of the link function

link_Hessian <- function(
    t,
    method = c("logistic", "poisson")) {

    method <- match.arg(method)

    if (method == "logistic") {

        hes <- exp(t) / (1 + exp(t))^2

    } else if (method == "poisson") {

        hes <- exp(t)
    }

    return(hes)
}

#--- COVARIANCE MATRIX OF ESTIMATING EQUATION ----------------------------------

#' Variance-covariance matrix of the estimation equation
#'
#' \code{Sigma_cal} function for variance-covariance matrix of the
#' estimation equation
#'
#' @param X_l Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_u Array or data.frame containing observed or predicted
#' covariates in unlabeled data.
#'
#' @param Y_l Array or data.frame of observed outcomes in labeled data.
#'
#' @param f_l Array or data.frame of predicted outcomes in labeled data.
#'
#' @param f_u Array or data.frame of predicted outcomes in
#' unlabeled data.
#'
#' @param w weights vector PSPA linear regression (d-dimensional, where
#' d equals the number of covariates).
#'
#' @param theta parameter theta
#'
#' @param quant quantile for quantile estimation
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return variance-covariance matrix of the estimation equation

Sigma_cal <- function(
    X_l,
    X_u,
    Y_l,
    f_l,
    f_u,
    w,
    theta,
    quant = NA,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    psi_Y_l <- psi(X_l, Y_l, theta, quant, method)
    psi_f_l <- psi(X_l, f_l, theta, quant, method)
    psi_f_u <- psi(X_u, f_u, theta, quant, method)

    n <- nrow(Y_l)
    N <- nrow(f_u)

    q <- if (method %in% c("mean", "quantile")) 1 else q <- ncol(X_l)

    M1 <- cov(t(psi_Y_l))

    M2 <- cov(t(psi_f_l))

    M3 <- cov(t(psi_f_u))

    M4 <- cov(t(psi_Y_l), t(psi_f_l))

    if (method == "ols") {

        A <- A(rbind(X_l, X_u), c(Y_l, f_u), quant, theta, method)

        A_inv <- solve(A)

    } else if (method %in% c("mean", "quantile")) {

        A <- A(X_l, Y_l, quant, theta, method)

        A_inv <- 1 / A

    } else {

        A_lab      <- A(X_l, Y_l, quant, theta, method)

        Ahat_lab   <- A(X_l, f_l, quant, theta, method)

        Ahat_unlab <- A(X_u, f_u, quant, theta, method)

        A <- A_lab + diag(w) %*% (Ahat_unlab - Ahat_lab)

        A_inv <- solve(A)
    }

    rho <- n / N

    if (method %in% c("mean", "quantile")) {

        Sigma <- A_inv * (M1 + w * (M2 + rho * M3) * w - 2 * w * M4) * A_inv

    } else {

        Sigma <- A_inv %*% (M1 + diag(w) %*% (M2 + rho * M3) %*% diag(w) -

            2 * diag(w) %*% M4) %*% A_inv
    }

    return(Sigma)
}

#--- ONE-STEP UPDATE -----------------------------------------------------------

#' One-step update for obtaining estimator
#'
#' \code{optim_est} function for One-step update for obtaining estimator
#'
#' @param X_l Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_u Array or data.frame containing observed or
#' predicted covariates in unlabeled data.
#'
#' @param Y_l Array or data.frame of observed outcomes in labeled data.
#'
#' @param f_l Array or data.frame of predicted outcomes in labeled data.
#'
#' @param f_u Array or data.frame of predicted outcomes in
#' unlabeled data.
#'
#' @param w weights vector PSPA linear regression (d-dimensional, where
#' d equals the number of covariates).
#'
#' @param theta parameter theta
#'
#' @param quant quantile for quantile estimation
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return estimator

optim_est <- function(
    X_l,
    X_u,
    Y_l,
    f_l,
    f_u,
    w,
    theta,
    quant = NA,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    if (method == "mean") {

        theta <- mean(Y_l) + as.vector(w) * (mean(f_u) - mean(f_l))

    } else if (method == "quantile") {

        A_lab <- A(X_l, Y_l, quant, theta, method)

        Ahat_lab <- A(X_l, f_l, quant, theta, method)

        Ahat_unlab <- A(X_u, f_u, quant, theta, method)

        A <- A_lab + w * (Ahat_unlab - Ahat_lab)

        A_inv <- 1 / A

        theta <- theta - A_inv * mean_psi_pop(X_l, X_u, Y_l, f_l, f_u, w,

            theta, quant = quant, method = method)

    } else if (method == "ols") {

        A <- A(rbind(X_l, X_u), c(Y_l, f_u), quant, theta, method)

        A_inv <- solve(A)

        theta <- theta - A_inv %*% mean_psi_pop(X_l, X_u, Y_l, f_l, f_u, w,

            theta, quant = quant, method = method)

    } else {

        A_lab <- A(X_l, Y_l, quant, theta, method)

        Ahat_lab <- A(X_l, f_l, quant, theta, method)

        Ahat_unlab <- A(X_u, f_u, quant, theta, method)

        A <- A_lab + diag(w) %*% (Ahat_unlab - Ahat_lab)

        A_inv <- solve(A)

        theta <- theta - A_inv %*% mean_psi_pop(X_l, X_u, Y_l, f_l, f_u, w,

            theta, quant = quant, method = method)
    }

    return(theta)
}

#--- ONE-STEP UPDATE - WEIGHT VECTOR -------------------------------------------

#' One-step update for obtaining the weight vector
#'
#' \code{optim_weights} function for One-step update for obtaining estimator
#'
#' @param X_l Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_u Array or data.frame containing observed or
#' predicted covariates in unlabeled data.
#'
#' @param Y_l Array or data.frame of observed outcomes in labeled data.
#'
#' @param f_l Array or data.frame of predicted outcomes in labeled data.
#'
#' @param f_u Array or data.frame of predicted outcomes in
#' unlabeled data.
#'
#' @param w weights vector PSPA linear regression (d-dimensional, where
#' d equals the number of covariates).
#'
#' @param theta parameter theta
#'
#' @param quant quantile for quantile estimation
#'
#' @param method indicates the method to be used for M-estimation.
#' Options include "mean", "quantile", "ols", "logistic", and "poisson".
#'
#' @return weights

optim_weights <- function(
    X_l,
    X_u,
    Y_l,
    f_l,
    f_u,
    w,
    theta,
    quant = NA,
    method = c("ols", "quantile", "mean", "logistic", "poisson")) {

    method <- match.arg(method)

    psi_Y_l <- psi(X_l, Y_l, theta, quant = quant, method = method)
    psi_f_l <- psi(X_l, f_l, theta, quant = quant, method = method)
    psi_f_u <- psi(X_u, f_u, theta, quant = quant, method = method)

    n <- nrow(Y_l)
    N <- nrow(f_u)

    q <- if (method %in% c("mean", "quantile")) 1 else q <- ncol(X_l)

    A_lab_inv <- solve(A(X_l, Y_l, quant, theta, method))

    M2 <- cov(t(psi_f_l))

    M3 <- cov(t(psi_f_u))

    M4 <- cov(t(psi_Y_l), t(psi_f_l))

    rho <- n / N

    w <- diag(A_lab_inv %*% M4 %*% A_lab_inv) /

        diag(A_lab_inv %*% (M2 + rho * M3) %*% A_lab_inv)

    w[w < 0] <- 0

    w[w > 1] <- 1

    return(w)
}

#--- SIMULATE PSPA DATA --------------------------------------------------------

#' Simulate the data for testing the functions
#'
#' \code{sim_data_y} for simulation with ML-predicted Y
#'
#' @param r imputation correlation
#'
#' @param binary simulate binary outcome or not
#'
#' @return simulated data
#'
#' @import MASS
#'
#' @importFrom randomForest randomForest

sim_data_y <- function(
    r = 0.9,
    binary = FALSE) {

    n_t <- 500
    n_l <- 500
    n_u <- 5000

    sigma_Y <- sqrt(5)

    mu <- c(0, 0)

    Sigma <- matrix(c(1, 0, 0, 1), 2, 2)

    n_data <- n_u + n_l + n_t

    data <- as.data.frame(MASS::mvrnorm(n_data, mu, Sigma))

    colnames(data) <- c("X1", "X2")

    beta_1 <- beta_2 <- r * sigma_Y / sqrt(2 * 3)

    data$epsilon <- rnorm(n_data, 0, sqrt(1 - r^2)) * sigma_Y

    data$Y <- data$X1 * beta_1 + data$X2 * beta_2 + data$X1^2 *

        beta_1 + data$X2^2 * beta_1 + data$epsilon

    if (binary) {

        data$Y <- ifelse(data$Y > median(unlist(data$Y)), 1, 0)

        dat_t <- data[seq_len(n_t), ]
        dat_l <- data[(n_t + 1):(n_l + n_t), ]
        dat_u <- data[(n_l + n_t + 1):n_data, ]

        dat_t$Y <- as.factor(dat_t$Y)

        train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = dat_t)

        dat_l$Y_hat <- predict(train_fit, newdata = dat_l)
        dat_u$Y_hat <- predict(train_fit, newdata = dat_u)

        X_l <- as.data.frame(dat_l$X1)
        X_u <- as.data.frame(dat_u$X1)
        Y_l <- as.data.frame(dat_l$Y)
        f_l <- as.data.frame(as.numeric(dat_l$Y_hat) - 1)
        f_u <- as.data.frame(as.numeric(dat_u$Y_hat) - 1)

    } else {

        dat_t <- data[seq_len(n_t), ]
        dat_l <- data[(n_t + 1):(n_l + n_t), ]
        dat_u <- data[(n_l + n_t + 1):n_data, ]

        train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = dat_t)

        dat_l$Y_hat <- predict(train_fit, newdata = dat_l)
        dat_u$Y_hat <- predict(train_fit, newdata = dat_u)

        X_l <- as.data.frame(dat_l$X1)
        X_u <- as.data.frame(dat_u$X1)
        Y_l <- as.data.frame(dat_l$Y)
        f_l <- as.data.frame(dat_l$Y_hat)
        f_u <- as.data.frame(dat_u$Y_hat)
    }

    colnames(X_l) <- "X1"
    colnames(X_u) <- "X1"
    colnames(Y_l) <- "Y"
    colnames(f_l) <- "Y_hat"
    colnames(f_u) <- "Y_hat"

    out <- list(X_l = X_l, X_u = X_u,Y_l = Y_l, f_l = f_l, f_u = f_u)

    return(out)
}
