#--- EMPIRICAL CDF -------------------------------------------------------------

#' Empirical CDF of the Data
#'
#' Computes the empirical CDF of the data.
#'
#' @param Y (matrix): n x 1 matrix of observed data.
#'
#' @param grid (matrix): Grid of values to compute the CDF at.
#'
#' @param w (vector, optional): n-vector of sample weights.
#'
#' @return (list): Empirical CDF and its standard deviation at the specified
#' grid points.

compute_cdf <- function(
    Y,
    grid,
    w = NULL) {

    n <- length(Y)

    if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

    indicators <- matrix((Y <= rep(grid, each = n)) * w, ncol = length(grid))

    cdf_mn <- apply(indicators, 2, mean)

    cdf_sd <- apply(indicators, 2, sd) * sqrt((n - 1) / n)

    return(list(cdf_mn, cdf_sd))
}

#--- DIFFERENCE IN EMPIRICAL CDFS ----------------------------------------------

#' Empirical CDF Difference
#'
#' Computes the difference between the empirical CDFs of the data and the
#' predictions.
#'
#' @param Y (matrix): n x 1 matrix of observed data.
#'
#' @param f (matrix): n x 1 matrix of predictions.
#'
#' @param grid (matrix): Grid of values to compute the CDF at.
#'
#' @param w (vector, optional): n-vector of sample weights.
#'
#' @return (list): Difference between the empirical CDFs of the data and the
#' predictions and its standard deviation at the specified grid points.

compute_cdf_diff <- function(
    Y,
    f,
    grid,
    w = NULL) {

    n <- length(Y)

    if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

    ind_Y <- matrix((Y <= rep(grid, each = n)) * w, ncol = length(grid))

    ind_f <- matrix((f <= rep(grid, each = length(f))) * w, ncol = length(grid))

    diff_mn <- apply(ind_Y - ind_f, 2, mean)

    diff_sd <- apply(ind_Y - ind_f, 2, sd) * sqrt((n - 1) / n)

    return(list(diff_mn, diff_sd))
}

#--- RECTIFIED CDF -------------------------------------------------------------

#' Rectified CDF
#'
#' Computes the rectified CDF of the data.
#'
#' @param Y_l (vector): Gold-standard labels.
#'
#' @param f_l (vector): Predictions corresponding to the gold-standard labels.
#'
#' @param f_u (vector): Predictions corresponding to the unlabeled data.
#'
#' @param grid (vector): Grid of values to compute the CDF at.
#'
#' @param w_l (vector, optional): Sample weights for the labeled data set.
#'
#' @param w_u (vector, optional): Sample weights for the unlabeled data set.
#'
#' @return (vector): Rectified CDF of the data at the specified grid points.

rectified_cdf <- function(
    Y_l,
    f_l,
    f_u,
    grid,
    w_l = NULL,
    w_u = NULL) {

    n <- length(Y_l)
    N <- length(f_u)

    if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

    if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

    cdf_f_u <- compute_cdf(f_u, grid, w = w_u)[[1]]

    cdf_rectifier <- compute_cdf_diff(Y_l, f_l, grid, w = w_l)[[1]]

    return(cdf_f_u + cdf_rectifier)
}

#--- RECTIFIED P-VALUE ---------------------------------------------------------

#' Rectified P-Value
#'
#' Computes a rectified p-value.
#'
#' @param rectifier (float or vector): Rectifier value.
#'
#' @param rectifier_std (float or vector): Rectifier standard deviation.
#'
#' @param imputed_mean (float or vector): Imputed mean.
#'
#' @param imputed_std (float or vector): Imputed standard deviation.
#'
#' @param null (float, optional): Value of the null hypothesis to be tested.
#' Defaults to `0`.
#'
#' @param alternative (str, optional): Alternative hypothesis, either
#' 'two-sided', 'larger' or 'smaller'.
#'
#' @return (float or vector): The rectified p-value.

rectified_p_value <- function(
    rectifier,
    rectifier_std,
    imputed_mean,
    imputed_std,
    null = 0,
    alternative = "two-sided") {

    rectified_point_estimate <- imputed_mean + rectifier

    rectified_std <- pmax(sqrt(imputed_std^2 + rectifier_std^2), 1e-16)

    p_value <- zstat_generic(rectified_point_estimate, 0, rectified_std,

        alternative, null)[[2]]

    return(p_value)
}

#--- ORDINARY LEAST SQUARES ----------------------------------------------------

#' Ordinary Least Squares
#'
#' @description
#' Computes the ordinary least squares coefficients.
#'
#' @param X (matrix): n x p matrix of covariates.
#'
#' @param Y (vector): p-vector of outcome values.
#'
#' @param return_se (bool, optional): Whether to return the standard errors of
#' the coefficients.
#'
#' @return (list): A list containing the following:
#'
#' \describe{
#'    \item{theta}{(vector): p-vector of ordinary least squares estimates of
#'    the coefficients.}
#'    \item{se}{(vector): If return_se == TRUE, return the p-vector of
#'    standard errors of the coefficients.}
#' }

ols <- function(
    X,
    Y,
    return_se = FALSE) {

    fit <- lm(Y ~ X - 1)

    theta <- coef(fit)

    if (return_se) {

        se <- sqrt(diag(vcov(fit)))

        return(list(theta = theta, se = se))

    } else {

        return(theta)
    }
}

#--- WEIGHTED LEAST SQUARES ----------------------------------------------------

#' Weighted Least Squares
#'
#' @description
#' Computes the weighted least squares estimate of the coefficients.
#'
#' @param X (matrix): n x p matrix of covariates.
#'
#' @param Y (vector): p-vector of outcome values.
#'
#' @param w (vector, optional): n-vector of sample weights.
#'
#' @param return_se (bool, optional): Whether to return the standard errors of
#' the coefficients.
#'
#' @return (list): A list containing the following:
#'
#' \describe{
#'    \item{theta}{(vector): p-vector of weighted least squares estimates of
#'    the coefficients.}
#'    \item{se}{(vector): If return_se == TRUE, return the p-vector of
#'    standard errors of the coefficients.}
#' }

wls <- function(
    X,
    Y,
    w = NULL,
    return_se = FALSE) {

    if (is.null(w) || all(w == 1)) {

        return(ols(X, Y, return_se = return_se))
    }

    fit <- lm(Y ~ X - 1, weights = w)

    theta <- coef(fit)

    if (return_se) {

        se <- sqrt(diag(vcov(fit)))

        return(list(theta = theta, se = se))

    } else {

        return(theta)
    }
}

#--- OLS GRADIENT AND HESSIAN --------------------------------------------------

#' OLS Gradient and Hessian
#'
#' @description
#' Computes the statistics needed for the OLS-based
#' prediction-powered inference.
#'
#' @param est (vector): Point estimates of the coefficients.
#'
#' @param X_l (matrix): Covariates for the labeled data set.
#'
#' @param Y_l (vector): Labels for the labeled data set.
#'
#' @param f_l (vector): Predictions for the labeled data set.
#'
#' @param X_u (matrix): Covariates for the unlabeled data set.
#'
#' @param f_u (vector): Predictions for the unlabeled data set.
#'
#' @param w_l (vector, optional): Sample weights for the labeled data set.
#'
#' @param w_u (vector, optional): Sample weights for the unlabeled data set.
#'
#' @param use_u (boolean, optional): Whether to use the unlabeled data set.
#'
#' @return (list): A list containing the following:
#'
#' \describe{
#'    \item{grads}{(matrix): n x p matrix gradient of the loss function with
#'    respect to the coefficients.}
#'    \item{grads_hat}{(matrix): n x p matrix gradient of the loss function
#'    with respect to the coefficients, evaluated using the labeled
#'    predictions.}
#'    \item{grads_hat_unlabeled}{(matrix): N x p matrix gradient of the loss
#'    function with respect to the coefficients, evaluated using the unlabeled
#'    predictions.}
#'    \item{inv_hessian}{(matrix): p x p matrix inverse Hessian of the loss
#'    function with respect to the coefficients.}
#' }

ols_get_stats <- function(
    est,
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    w_l = NULL,
    w_u = NULL,
    use_u = TRUE) {

    n <- nrow(f_l)
    N <- nrow(f_u)
    p <- ncol(X_l)

    if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

    if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

    hessian <- matrix(0, nrow = p, ncol = p)

    grads_hat_unlabeled <- matrix(0, nrow = N, ncol = p)

    if (use_u) {

        for (i in seq_len(N)) {

            hessian <- hessian + w_u[i] / (N + n) * tcrossprod(X_u[i, ])

            grads_hat_unlabeled[i, ] <- w_u[i] * X_u[i, ] *

                (sum(X_u[i, ] * est) - f_u[i])
        }
    }

    grads <- matrix(0, nrow = n, ncol = p)

    grads_hat <- matrix(0, nrow = n, ncol = p)

    for (i in seq_len(n)) {

        if (use_u) {

            hessian <- hessian + w_l[i] / (N + n) * tcrossprod(X_l[i, ])

        } else {

            hessian <- hessian + w_l[i] / n * tcrossprod(X_l[i, ])
        }

        grads[i, ] <- w_l[i] * X_l[i, ] * (sum(X_l[i, ] * est) - Y_l[i])

        grads_hat[i, ] <- w_l[i] * X_l[i, ] * (sum(X_l[i, ] * est) - f_l[i])
    }

    inv_hessian <- solve(hessian)

    return(list(grads = grads, grads_hat = grads_hat,

        grads_hat_unlabeled = grads_hat_unlabeled, inv_hessian = inv_hessian))
}

#--- ESTIMATE POWER TUNING PARAMETER -------------------------------------------

#' Estimate PPI++ Power Tuning Parameter
#'
#' @description
#' Calculates the optimal value of lhat for the prediction-powered confidence
#' interval for GLMs.
#'
#' @param grads (matrix): n x p matrix gradient of the loss function with
#' respect to the parameter evaluated at the labeled data.
#'
#' @param grads_hat (matrix): n x p matrix gradient of the loss function with
#' respect to the model parameter evaluated using predictions on the labeled
#' data.
#'
#' @param grads_hat_unlabeled (matrix): N x p matrix gradient of the loss
#' function with respect to the parameter evaluated using predictions on the
#' unlabeled data.
#'
#' @param inv_hessian (matrix): p x p matrix inverse of the Hessian of the
#' loss function with respect to the parameter.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`.
#' If `None`, it optimizes the total variance over all coordinates.
#' Must be in (1, ..., d) where d is the shape of the estimand.
#'
#' @param clip (boolean, optional): Whether to clip the value of lhat to be
#' non-negative. Defaults to `False`.
#'
#' @return (float): Optimal value of `lhat` in \[0,1\].

calc_lhat_glm <- function(
    grads,
    grads_hat,
    grads_hat_unlabeled,
    inv_hessian,
    coord = NULL,
    clip = FALSE) {

    if (is.null(dim(grads))) {

        grads <- matrix(grads, ncol = 1)
    }

    if (is.null(dim(grads_hat))) {

        grads_hat <- matrix(grads_hat, ncol = 1)
    }

    if (is.null(dim(grads_hat_unlabeled))) {

        grads_hat_unlabeled <- matrix(grads_hat_unlabeled, ncol = 1)
    }

    n <- nrow(grads)
    N <- nrow(grads_hat_unlabeled)
    p <- ncol(inv_hessian)

    cov_grads <- matrix(0, nrow = p, ncol = p)

    for (i in seq_len(n)) {

        cov_grads <- cov_grads + (1 / n) * (

            outer(grads[i, ] - colMeans(grads),

                grads_hat[i, ] - colMeans(grads_hat)) +

            outer(grads_hat[i, ] - colMeans(grads_hat),

                grads[i, ] - colMeans(grads)))
    }

    var_grads_hat <- cov(rbind(grads_hat, grads_hat_unlabeled))

    vhat <- if (is.null(coord)) {

        inv_hessian

    } else {

        inv_hessian %*% diag(p)[, coord]
    }

    if (p > 1) {

        num <- ifelse(is.null(coord),

            sum(diag(vhat %*% cov_grads %*% vhat)),

            vhat %*% cov_grads %*% vhat)

        denom <- ifelse(is.null(coord),

            2 * (1 + n / N) * sum(diag(vhat %*% var_grads_hat %*% vhat)),

            2 * (1 + n / N) * vhat %*% var_grads_hat %*% vhat)
    } else {

        num <- vhat * cov_grads * vhat

        denom <- 2 * (1 + n / N) * vhat * var_grads_hat * vhat
    }

    lhat <- num / denom

    if (clip) {

        lhat <- pmax(0, pmin(lhat, 1))
    }

    return(as.numeric(lhat))
}

#--- NORMAL Z-STATISTIC --------------------------------------------------------

#' Compute Z-Statistic and P-Value
#'
#' @description
#' Computes the z-statistic and the corresponding p-value for a given test.
#'
#' @param value1 (numeric): The first value or sample mean.
#'
#' @param value2 (numeric): The second value or sample mean.
#'
#' @param std_diff (numeric): The standard error of the difference between the
#' two values.
#'
#' @param alternative (character): The alternative hypothesis. Can be one of
#' "two-sided" (or "2-sided", "2s"), "larger" (or "l"), or "smaller" (or "s").
#'
#' @param diff (numeric, optional): The hypothesized difference between the
#' two values. Default is 0.
#'
#' @return (list): A list containing the following:
#' \describe{
#'    \item{zstat}{(numeric): The computed z-statistic.}
#'    \item{pvalue}{(numeric): The corresponding p-value for the test.}
#' }

zstat_generic <- function(
    value1,
    value2,
    std_diff,
    alternative,
    diff = 0) {

    zstat <- (value1 - value2 - diff) / std_diff

    if (alternative %in% c("two-sided", "2-sided", "2s")) {

        pvalue <- 2 * (1 - pnorm(abs(zstat)))

    } else if (alternative %in% c("larger", "l")) {

        pvalue <- 1 - pnorm(zstat)

    } else if (alternative %in% c("smaller", "s")) {

        pvalue <- pnorm(zstat)

    } else {

        stop("Invalid alternative")
    }

    return(list(zstat = zstat, pvalue = pvalue))
}

#--- NORMAL CONFIDENCE INTERVALS -----------------------------------------------

#' Normal Confidence Intervals
#'
#' @description
#' Calculates normal confidence intervals for a given alternative at a given
#' significance level.
#'
#' @param mean (float): Estimated normal mean.
#'
#' @param std_mean (float): Estimated standard error of the mean.
#'
#' @param alpha (float): Significance level in \[0,1\]
#'
#' @param alternative (string): Alternative hypothesis, either 'two-sided',
#' 'larger' or 'smaller'.
#'
#' @return (vector): Lower and upper (1 - alpha) * 100% confidence limits.

zconfint_generic <- function(
    mean,
    std_mean,
    alpha,
    alternative) {

    if (alternative %in% c("two-sided", "2-sided", "2s")) {

        zcrit <- qnorm(1 - alpha / 2)
        lower <- mean - zcrit * std_mean
        upper <- mean + zcrit * std_mean

    } else if (alternative %in% c("larger", "l")) {

        zcrit <- qnorm(alpha)
        lower <- mean + zcrit * std_mean
        upper <- Inf

    } else if (alternative %in% c("smaller", "s")) {

        zcrit <- qnorm(1 - alpha)
        lower <- -Inf
        upper <- mean + zcrit * std_mean

    } else {

        stop("Invalid alternative")
    }

    return(cbind(lower, upper))
}

#--- LOG1P EXPONENTIAL ---------------------------------------------------------

#' Log1p Exponential
#'
#' @description
#' Computes the natural logarithm of 1 plus the exponential of the input,
#' to handle large inputs.
#'
#' @param x (vector): A numeric vector of inputs.
#'
#' @return (vector): A numeric vector where each element is the result of
#' log(1 + exp(x)).

log1pexp <- function(x) {

    idxs <- x > 10

    out <- numeric(length(x))

    out[idxs] <- x[idxs]

    out[!idxs] <- log1p(exp(x[!idxs]))

    return(out)
}

#=== PPI++ LOGISTIC REGRESSION =================================================

#=== PPI++ LOGISTIC GRADIENT AND HESSIAN =======================================

#' Logistic Regression Gradient and Hessian
#'
#' @description
#' Computes the statistics needed for the logstic regression-based
#' prediction-powered inference.
#'
#' @param est (vector): Point estimates of the coefficients.
#'
#' @param X_l (matrix): Covariates for the labeled data set.
#'
#' @param Y_l (vector): Labels for the labeled data set.
#'
#' @param f_l (vector): Predictions for the labeled data set.
#'
#' @param X_u (matrix): Covariates for the unlabeled data set.
#'
#' @param f_u (vector): Predictions for the unlabeled data set.
#'
#' @param w_l (vector, optional): Sample weights for the labeled data set.
#'
#' @param w_u (vector, optional): Sample weights for the unlabeled data set.
#'
#' @param use_u (bool, optional): Whether to use the unlabeled data set.
#'
#' @return (list): A list containing the following:
#'
#' \describe{
#'    \item{grads}{(matrix): n x p matrix gradient of the loss function with
#'    respect to the coefficients.}
#'    \item{grads_hat}{(matrix): n x p matrix gradient of the loss function
#'    with respect to the coefficients, evaluated using the labeled
#'    predictions.}
#'    \item{grads_hat_unlabeled}{(matrix): N x p matrix gradient of the loss
#'    function with respect to the coefficients, evaluated using the unlabeled
#'    predictions.}
#'    \item{inv_hessian}{(matrix): p x p matrix inverse Hessian of the loss
#'    function with respect to the coefficients.}
#' }

logistic_get_stats <- function(
    est,
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    w_l = NULL,
    w_u = NULL,
    use_u = TRUE) {

    n <- nrow(f_l)
    N <- nrow(f_u)
    p <- ncol(X_u)

    if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

    if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

    mu_l <- plogis(X_l %*% est)

    mu_u <- plogis(X_u %*% est)

    hessian <- matrix(0, nrow = p, ncol = p)

    grads_hat_unlabeled <- matrix(0, nrow = N, ncol = p)

    if (use_u) {

        for (i in seq_len(N)) {

            hessian <- hessian + w_u[i] / (N + n) * mu_u[i] * (1 - mu_u[i]) *

                tcrossprod(X_u[i, ])

            grads_hat_unlabeled[i, ] <- w_u[i] * X_u[i, ] * (mu_u[i] - f_u[i])
        }
    }

    grads <- matrix(0, nrow = n, ncol = p)

    grads_hat <- matrix(0, nrow = n, ncol = p)

    for (i in seq_len(n)) {

        if (use_u) {

            hessian <- hessian + w_l[i] / (N + n) * mu_l[i] * (1 - mu_l[i]) *

                tcrossprod(X_l[i, ])

        } else {

            hessian <- hessian + w_l[i] / n * mu_l[i] * (1 - mu_l[i]) *

                tcrossprod(X_l[i, ])
        }

        grads[i, ] <- w_l[i] * X_l[i, ] * (mu_l[i] - Y_l[i])

        grads_hat[i, ] <- w_l[i] * X_l[i, ] * (mu_l[i] - f_l[i])
    }

    inv_hessian <- solve(hessian)

    return(list(grads = grads, grads_hat = grads_hat,

        grads_hat_unlabeled = grads_hat_unlabeled, inv_hessian = inv_hessian))
}
