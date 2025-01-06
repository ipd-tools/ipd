#===============================================================================
# HELPERS
#===============================================================================

#=== PPI++ QUANTILE ESTIMATION =================================================

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
#'
#' @examples
#'
#' Y <- c(1, 2, 3, 4, 5)
#'
#' grid <- seq(0, 6, by = 0.5)
#'
#' compute_cdf(Y, grid)
#'
#' @export

compute_cdf <- function(Y, grid, w = NULL) {

  n <- length(Y)

  if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

  indicators <- matrix((Y <= rep(grid, each = n)) * w,

                       ncol = length(grid))

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
#'
#' @examples
#'
#' Y <- c(1, 2, 3, 4, 5)
#'
#' f <- c(1.1, 2.2, 3.3, 4.4, 5.5)
#'
#' grid <- seq(0, 6, by = 0.5)
#'
#' compute_cdf_diff(Y, f, grid)
#'
#' @export

compute_cdf_diff <- function(Y, f, grid, w = NULL) {

  n <- length(Y)

  if (is.null(w)) w <- rep(1, n) else w <- w / sum(w) * n

  indicators_Y <- matrix((Y <= rep(grid, each = n)) * w,

    ncol = length(grid))

  indicators_f <- matrix((f <= rep(grid, each = length(f))) * w,

    ncol = length(grid))

  diff_mn <- apply(indicators_Y - indicators_f, 2, mean)

  diff_sd <- apply(indicators_Y - indicators_f, 2, sd) * sqrt((n - 1) / n)

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
#'
#' @examples
#'
#' Y_l <- c(1, 2, 3, 4, 5)
#'
#' f_l <- c(1.1, 2.2, 3.3, 4.4, 5.5)
#'
#' f_u <- c(1.2, 2.3, 3.4)
#'
#' grid <- seq(0, 6, by = 0.5)
#'
#' rectified_cdf(Y_l, f_l, f_u, grid)
#'
#' @export

rectified_cdf <- function(Y_l, f_l, f_u, grid, w_l = NULL, w_u = NULL) {

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
#'
#' @examples
#'
#' rectifier <- 0.7
#'
#' rectifier_std <- 0.5
#'
#' imputed_mean <- 1.5
#'
#' imputed_std <- 0.3
#'
#' rectified_p_value(rectifier, rectifier_std, imputed_mean, imputed_std)
#'
#' @export

rectified_p_value <- function(rectifier, rectifier_std,

  imputed_mean, imputed_std, null = 0, alternative = "two-sided") {

  rectified_point_estimate <- imputed_mean + rectifier

  rectified_std <- pmax(sqrt(imputed_std^2 + rectifier_std^2), 1e-16)

  p_value <- zstat_generic(rectified_point_estimate, 0, rectified_std,

    alternative, null)[[2]]

  return(p_value)
}

#=== PPI++ OLS =================================================================

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
#'
#' @examples
#'
#' n <- 1000
#'
#' X <- rnorm(n, 1, 1)
#'
#' Y <- X + rnorm(n, 0, 1)
#'
#' ols(X, Y, return_se = TRUE)
#'
#' @export

ols <- function(X, Y, return_se = FALSE) {

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
#'
#' @examples
#'
#' n <- 1000
#'
#' X <- rnorm(n, 1, 1)
#'
#' w <- rep(1, n)
#'
#' Y <- X + rnorm(n, 0, 1)
#'
#' wls(X, Y, w = w, return_se = TRUE)
#'
#' @export

wls <- function(X, Y, w = NULL, return_se = FALSE) {

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
#' est <- ppi_plusplus_ols_est(X_l, Y_l, f_l, X_u, f_u)
#'
#' stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u, use_u = TRUE)
#'
#' @export

ols_get_stats <- function(est, X_l, Y_l, f_l, X_u, f_u,

  w_l = NULL, w_u = NULL, use_u = TRUE) {

  n <- nrow(f_l)

  N <- nrow(f_u)

  p <- ncol(X_l)

  if (is.null(w_l)) w_l <- rep(1, n) else w_l <- w_l / sum(w_l) * n

  if (is.null(w_u)) w_u <- rep(1, N) else w_u <- w_u / sum(w_u) * N

  hessian <- matrix(0, nrow = p, ncol = p)

  grads_hat_unlabeled <- matrix(0, nrow = N, ncol = p)

  if (use_u) {

    for (i in 1:N) {

      hessian <- hessian + w_u[i] / (N + n) * tcrossprod(X_u[i, ])

      grads_hat_unlabeled[i, ] <- w_u[i] * X_u[i, ] *

        (sum(X_u[i, ] * est) - f_u[i])
    }
  }

  grads <- matrix(0, nrow = n, ncol = p)

  grads_hat <- matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {

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
#' est <- ppi_plusplus_ols_est(X_l, Y_l, f_l, X_u, f_u)
#'
#' stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u)
#'
#' calc_lhat_glm(stats$grads, stats$grads_hat, stats$grads_hat_unlabeled,
#'
#'   stats$inv_hessian, coord = NULL, clip = FALSE)
#'
#' @return (float): Optimal value of `lhat` in \[0,1\].
#'
#' @export

calc_lhat_glm <- function(grads, grads_hat, grads_hat_unlabeled, inv_hessian,

  coord = NULL, clip = FALSE) {

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

  for (i in 1:n) {

    cov_grads <- cov_grads + (1/n) * (outer(grads[i,] - colMeans(grads),

      grads_hat[i,] - colMeans(grads_hat)) +

        outer(grads_hat[i,] - colMeans(grads_hat),

          grads[i,] - colMeans(grads)))
  }

  var_grads_hat <- cov(rbind(grads_hat, grads_hat_unlabeled))

  if (is.null(coord)) {

    vhat <- inv_hessian

  } else {

    vhat <- inv_hessian %*% diag(p)[, coord]
  }

  if (p > 1) {

    num <- ifelse(is.null(coord),

      sum(diag(vhat %*% cov_grads %*% vhat)), vhat %*% cov_grads %*% vhat)

    denom <- ifelse(is.null(coord),

      2*(1 + n/N) * sum(diag(vhat %*% var_grads_hat %*% vhat)),

      2*(1 + n/N) * vhat %*% var_grads_hat %*% vhat)

  } else {

    num <- vhat * cov_grads * vhat

    denom <- 2*(1 + n/N) * vhat * var_grads_hat * vhat
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
#'
#' @examples
#'
#' value1 <- 1.5
#'
#' value2 <- 1.0
#'
#' std_diff <- 0.2
#'
#' alternative <- "two-sided"
#'
#' result <- zstat_generic(value1, value2, std_diff, alternative)
#'
#' @export

zstat_generic <- function(value1, value2, std_diff, alternative, diff = 0) {

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
#'
#' @examples
#'
#' n <- 1000
#'
#' Y <- rnorm(n, 1, 1)
#'
#' se_Y <-  sd(Y) / sqrt(n)
#'
#' zconfint_generic(Y, se_Y, alpha = 0.05, alternative = "two-sided")
#'
#' @export

zconfint_generic <- function(mean, std_mean, alpha, alternative) {

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
#'
#' @examples
#'
#' x <- c(-1, 0, 1, 10, 100)
#'
#' log1pexp(x)
#'
#' @export

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
#' est <- ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u)
#'
#' stats <- logistic_get_stats(est, X_l, Y_l, f_l, X_u, f_u)
#'
#' @export

logistic_get_stats <- function(est, X_l, Y_l, f_l, X_u, f_u,

  w_l = NULL, w_u = NULL, use_u = TRUE) {

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

    for (i in 1:N) {

      hessian <- hessian + w_u[i] / (N + n) * mu_u[i] * (1 - mu_u[i]) *

        tcrossprod(X_u[i, ])

      grads_hat_unlabeled[i, ] <- w_u[i] * X_u[i, ] * (mu_u[i] - f_u[i])
    }
  }

  grads <- matrix(0, nrow = n, ncol = p)

  grads_hat <- matrix(0, nrow = n, ncol = p)

  for (i in 1:n) {

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

  return(

    list(grads = grads, grads_hat = grads_hat,

         grads_hat_unlabeled = grads_hat_unlabeled, inv_hessian = inv_hessian))
}

#=== PSPA ======================================================================

#-------------------------------------------------------------------------------
# PSPA: Post-Prediction Adaptive Inference
# Jiacheng Miao, Xinran Miao, Yixuan Wu, Jiwei Zhao, and Qiongshi Lu
# Available from https://arxiv.org/abs/2311.14220
#-------------------------------------------------------------------------------

#--- PSPA M-ESTIMATION FOR PREDICTED -------------------------------------------

#' PSPA M-Estimation for ML-predicted labels
#'
#' \code{pspa_y} function conducts post-prediction M-Estimation.
#'
#' @param X_lab Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_unlab Array or data.frame containing observed or predicted
#' covariates in unlabeled data.
#'
#' @param Y_lab Array or data.frame of observed outcomes in
#' labeled data.
#'
#' @param Yhat_lab Array or data.frame of predicted outcomes in
#' labeled data.
#'
#' @param Yhat_unlab Array or data.frame of predicted outcomes in
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
#'
#' @examples
#'
#' data <- sim_data_y()
#'
#' X_lab <- data$X_lab
#'
#' X_unlab <- data$X_unlab
#'
#' Y_lab <- data$Y_lab
#'
#' Yhat_lab <- data$Yhat_lab
#'
#' Yhat_unlab <- data$Yhat_unlab
#'
#' pspa_y(X_lab = X_lab, X_unlab = X_unlab,
#'
#'  Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
#'
#'  alpha = 0.05, method = "ols")
#'
#' @export

pspa_y <- function(X_lab = NA, X_unlab = NA, Y_lab, Yhat_lab, Yhat_unlab,

  alpha = 0.05, weights = NA, quant = NA, intercept = FALSE, method) {

  ## Common Values

  if (method %in% c("ols", "logistic", "poisson")) {

    if (intercept) {

      X_lab <- as.matrix(X_lab)
      X_unlab <- as.matrix(X_unlab)

    } else {

      X_lab <- cbind(1, as.matrix(X_lab))
      X_unlab <- cbind(1, as.matrix(X_unlab))
    }
  }

  Y_lab <- as.matrix(as.numeric(unlist(Y_lab)))

  Yhat_lab <- as.matrix(as.numeric(unlist(Yhat_lab)))

  Yhat_unlab <- as.matrix(as.numeric(unlist(Yhat_unlab)))

  n <- nrow(Y_lab)

  N <- nrow(Yhat_unlab)

  if (method %in% c("mean", "quantile")) {

    q <- 1

  } else {

    q <- ncol(X_lab)
  }

  ## Initial Values
  est <- est_ini(X_lab, Y_lab, quant, method)

  if (is.na(sum(weights))) {

    current_w <- rep(0, q)

    optimized_weight <- optim_weights(X_lab, X_unlab,

      Y_lab, Yhat_lab, Yhat_unlab, current_w, est, quant, method)

    current_w <- optimized_weight

  } else {

    current_w <- weights
  }

  ## Calculate Final Coefficients and Standard Errors

  est <- optim_est(X_lab, X_unlab,

    Y_lab, Yhat_lab, Yhat_unlab, current_w, est, quant, method)

  final_sigma <- Sigma_cal(X_lab, X_unlab,

    Y_lab, Yhat_lab, Yhat_unlab, current_w, est, quant, method)

  standard_errors <- sqrt(diag(final_sigma) / n)

  p_values <- 2 * pnorm(abs(est / standard_errors), lower.tail = FALSE)

  ## Create Output Table

  lower_ci <- est - qnorm(1 - alpha / 2) * standard_errors

  upper_ci <- est + qnorm(1 - alpha / 2) * standard_errors

  output_table <- data.frame(

    Estimate = est, Std.Error = standard_errors,
    Lower.CI = lower_ci, Upper.CI = upper_ci,
    P.value = p_values, Weight = current_w)

  colnames(output_table) <- c(

    "Estimate", "Std.Error", "Lower.CI", "Upper.CI", "P.value", "Weight")

  return(output_table)
}

#=== PSPA UTILS ================================================================

#-------------------------------------------------------------------------------
# General functions for M-estimation and Z-estimation
#-------------------------------------------------------------------------------

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
#'
#' @export

A <- function(X, Y, quant = NA, theta, method) {

  if (method == "ols") {

    n <- nrow(X)

    A <- (1 / n) * t(X) %*% X

  } else if (method == "quantile") {

    ## Kernel Density Estimation

    A <- sapply(theta,

      function(a, b) density(b, from = a, to = a, n = 1)$y, unlist(Y))

  } else if (method == "mean") {

    A <- 1

  } else if (method %in% c("logistic", "poisson")) {

    n <- nrow(X)

    mid <- sqrt(diag(as.vector(link_Hessian(X %*% theta, method)))) %*% X

    A <- 1 / n * t(mid) %*% mid
  }

  return(A)
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
#' @return initial estimatior
#'
#' @export

est_ini <- function(X, Y, quant = NA, method) {

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
#' @return esimating equation
#'
#' @export

psi <- function(X, Y, theta, quant = NA, method) {

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
#'
#' @export

mean_psi <- function(X, Y, theta, quant = NA, method) {

  psi <- as.matrix(rowMeans(psi(X, Y, theta, quant, method)))

  return(psi)
}

#--- SAMPLE EXPECTATION OF PSPA PSI --------------------------------------------

#' Sample expectation of PSPA psi
#'
#' \code{mean_psi_pop} function for sample expectation of PSPA psi
#'
#' @param X_lab Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_unlab Array or data.frame containing observed or predicted
#' covariates in unlabeled data.
#'
#' @param Y_lab Array or data.frame of observed outcomes in labeled data.
#'
#' @param Yhat_lab Array or data.frame of predicted outcomes in labeled data.
#'
#' @param Yhat_unlab Array or data.frame of predicted outcomes in
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
#'
#' @export

mean_psi_pop <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab,

  w, theta, quant = NA, method) {

  if (method %in% c("ols", "logistic", "poisson")) {

    psi_pop <- mean_psi(X_lab, Y_lab, theta, quant, method) + diag(w) %*% (

      mean_psi(X_unlab, Yhat_unlab, theta, quant, method) -

      mean_psi(X_lab,   Yhat_lab,   theta, quant, method))

  } else if (method %in% c("mean", "quantile")) {

    psi_pop <- mean_psi(X_lab, Y_lab, theta, quant, method) + w * (

      mean_psi(X_unlab, Yhat_unlab, theta, quant, method) -

      mean_psi(X_lab,   Yhat_lab,   theta, quant, method))
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
#'
#' @export

link_grad <- function(t, method) {

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
#'
#' @export

link_Hessian <- function(t, method) {

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
#' @param X_lab Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_unlab Array or data.frame containing observed or predicted
#' covariates in unlabeled data.
#'
#' @param Y_lab Array or data.frame of observed outcomes in labeled data.
#'
#' @param Yhat_lab Array or data.frame of predicted outcomes in labeled data.
#'
#' @param Yhat_unlab Array or data.frame of predicted outcomes in
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
#'
#' @export

Sigma_cal <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab,

  w, theta, quant = NA, method) {

  psi_y_lab <- psi(

    X_lab, Y_lab, theta, quant = quant, method = method)

  psi_yhat_lab <- psi(

    X_lab, Yhat_lab, theta, quant = quant, method = method)

  psi_yhat_unlab <- psi(

    X_unlab, Yhat_unlab, theta, quant = quant, method = method)

  n <- nrow(Y_lab)

  N <- nrow(Yhat_unlab)

  if (method %in% c("mean", "quantile")) {

    q <- 1

  } else {

    q <- ncol(X_lab)
  }

  M1 <- cov(t(psi_y_lab))

  M2 <- cov(t(psi_yhat_lab))

  M3 <- cov(t(psi_yhat_unlab))

  M4 <- cov(t(psi_y_lab), t(psi_yhat_lab))

  if (method == "ols") {

    A <- A(rbind(X_lab, X_unlab), c(Y_lab, Yhat_unlab), quant, theta, method)

    A_inv <- solve(A)

  } else if (method %in% c("mean", "quantile")) {

    A <- A(X_lab, Y_lab, quant, theta, method)

    A_inv <- 1/A

  } else {

    A_lab <- A(X_lab, Y_lab, quant, theta, method)

    Ahat_lab <- A(X_lab, Yhat_lab, quant, theta, method)

    Ahat_unlab <- A(X_unlab, Yhat_unlab, quant, theta, method)

    A <- A_lab + diag(w) %*% (Ahat_unlab - Ahat_lab)

    A_inv <- solve(A)
  }

  rho <- n / N

  if (method %in% c("mean", "quantile")) {

    Sigma <- A_inv * (M1 + w * (M2 + rho * M3) * w - 2 * w * M4) * A_inv

  } else {

    Sigma <- A_inv %*% (M1 + diag(w) %*%

      (M2 + rho * M3) %*% diag(w) - 2 * diag(w) %*% M4) %*% A_inv
  }

  return(Sigma)
}

#--- ONE-STEP UPDATE -----------------------------------------------------------

#' One-step update for obtaining estimator
#'
#' \code{optim_est} function for One-step update for obtaining estimator
#'
#' @param X_lab Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_unlab Array or data.frame containing observed or
#' predicted covariates in unlabeled data.
#'
#' @param Y_lab Array or data.frame of observed outcomes in labeled data.
#'
#' @param Yhat_lab Array or data.frame of predicted outcomes in labeled data.
#'
#' @param Yhat_unlab Array or data.frame of predicted outcomes in
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
#'
#' @export

optim_est <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab,

  w, theta, quant = NA, method) {

  if (method == "mean") {

    theta <- mean(Y_lab) + as.vector(w) * (mean(Yhat_unlab) - mean(Yhat_lab))

  } else if (method == "quantile") {

    A_lab <- A(X_lab, Y_lab, quant, theta, method)

    Ahat_lab <- A(X_lab, Yhat_lab, quant, theta, method)

    Ahat_unlab <- A(X_unlab, Yhat_unlab, quant, theta, method)

    A <- A_lab + w * (Ahat_unlab - Ahat_lab)

    A_inv <- 1 / A

    theta <- theta - A_inv * mean_psi_pop(

      X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab,

      w, theta, quant = quant, method = method)

  } else if (method == "ols") {

    A <- A(rbind(X_lab, X_unlab), c(Y_lab, Yhat_unlab), quant, theta, method)

    A_inv <- solve(A)

    theta <- theta - A_inv %*% mean_psi_pop(

      X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab,

      w, theta, quant = quant, method = method)

  } else {

    A_lab <- A(X_lab, Y_lab, quant, theta, method)

    Ahat_lab <- A(X_lab, Yhat_lab, quant, theta, method)

    Ahat_unlab <- A(X_unlab, Yhat_unlab, quant, theta, method)

    A <- A_lab + diag(w) %*% (Ahat_unlab - Ahat_lab)

    A_inv <- solve(A)

    theta <- theta - A_inv %*% mean_psi_pop(

      X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab,

      w, theta, quant = quant, method = method)
  }

  return(theta)
}

#--- ONE-STEP UPDATE - WEIGHT VECTOR -------------------------------------------

#' One-step update for obtaining the weight vector
#'
#' \code{optim_weights} function for One-step update for obtaining estimator
#'
#' @param X_lab Array or data.frame containing observed covariates in
#' labeled data.
#'
#' @param X_unlab Array or data.frame containing observed or
#' predicted covariates in unlabeled data.
#'
#' @param Y_lab Array or data.frame of observed outcomes in labeled data.
#'
#' @param Yhat_lab Array or data.frame of predicted outcomes in labeled data.
#'
#' @param Yhat_unlab Array or data.frame of predicted outcomes in unlabeled data.
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
#'
#' @export

optim_weights <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab,

  w, theta, quant = NA, method) {

  ## Objective Function

  psi_y_lab <- psi(

    X_lab, Y_lab, theta, quant = quant, method = method)

  psi_yhat_lab <- psi(

    X_lab, Yhat_lab, theta, quant = quant, method = method)

  psi_yhat_unlab <- psi(

    X_unlab, Yhat_unlab, theta, quant = quant, method = method)

  n <- nrow(Y_lab)

  N <- nrow(Yhat_unlab)

  if (method %in% c("mean", "quantile")) {

    q <- 1

  } else {

    q <- ncol(X_lab)
  }

  A_lab_inv <- solve(A(X_lab, Y_lab, quant, theta, method))

  M2 <- cov(t(psi_yhat_lab))

  M3 <- cov(t(psi_yhat_unlab))

  M4 <- cov(t(psi_y_lab), t(psi_yhat_lab))

  rho <- n / N

  w <- diag(A_lab_inv %*% M4 %*% A_lab_inv) /

    diag(A_lab_inv %*% (M2 + rho * M3) %*% A_lab_inv)

  ## Lower Bound: 0; Upper Bound: 1

  w[w < 0] <- 0

  w[w > 1] <- 1

  return(w)
}

#=== SIMULATE PSPA DATA ========================================================

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
#' @import randomForest MASS
#'
#' @export

sim_data_y <- function(r = 0.9, binary = FALSE) {

  ## Input Parameters

  n_train <- 500

  n_lab <- 500

  n_unlab <- 5000

  sigma_Y <- sqrt(5)

  ## Simulate Data

  mu <- c(0, 0) # Mean Vector

  Sigma <- matrix(c(1, 0, 0, 1), 2, 2) # Covariance Matrix

  n_data <- n_unlab + n_lab + n_train

  data <- as.data.frame(MASS::mvrnorm(n_data, mu, Sigma))

  colnames(data) <- c("X1", "X2")

  beta_1 <- beta_2 <- r * sigma_Y / sqrt(2 * 3)

  data$epsilon <- rnorm(n_data, 0, sqrt(1 - r^2)) * sigma_Y

  data$Y <- data$X1 * beta_1 + data$X2 * beta_2 + data$X1^2 *

    beta_1 + data$X2^2 * beta_1 + data$epsilon

  if (binary) {

    data$Y <- ifelse(data$Y > median(unlist(data$Y)), 1, 0)

    ## Split Data

    train_data <- data[1:n_train, ]

    lab_data <- data[(n_train + 1):(n_lab + n_train), ]

    unlab_data <- data[(n_lab + n_train + 1):n_data, ]

    ## Fit Machine Learning Model

    train_data$Y <- as.factor(train_data$Y)

    train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)

    lab_data$Y_hat <- predict(train_fit, newdata = lab_data)

    unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

    X_lab <- as.data.frame(lab_data$X1)

    X_unlab <- as.data.frame(unlab_data$X1)

    Y_lab <- as.data.frame(lab_data$Y)

    Yhat_lab <- as.data.frame(as.numeric(lab_data$Y_hat) - 1)

    Yhat_unlab <- as.data.frame(as.numeric(unlab_data$Y_hat) - 1)

    colnames(X_lab) <- "X1"

    colnames(X_unlab) <- "X1"

    colnames(Y_lab) <- "Y"

    colnames(Yhat_lab) <- "Y_hat"

    colnames(Yhat_unlab) <- "Y_hat"

    out <- list(

      X_lab = X_lab, X_unlab = X_unlab,

      Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab)

  } else {

    ## Split Data

    train_data <- data[1:n_train, ]

    lab_data <- data[(n_train + 1):(n_lab + n_train), ]

    unlab_data <- data[(n_lab + n_train + 1):n_data, ]

    # Fit Machine Learning Model

    train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)

    lab_data$Y_hat <- predict(train_fit, newdata = lab_data)

    unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

    X_lab <- as.data.frame(lab_data$X1)

    X_unlab <- as.data.frame(unlab_data$X1)

    Y_lab <- as.data.frame(lab_data$Y)

    Yhat_lab <- as.data.frame(lab_data$Y_hat)

    Yhat_unlab <- as.data.frame(unlab_data$Y_hat)

    colnames(X_lab) <- "X1"

    colnames(X_unlab) <- "X1"

    colnames(Y_lab) <- "Y"

    colnames(Yhat_lab) <- "Y_hat"

    colnames(Yhat_unlab) <- "Y_hat"

    out <- list(

      X_lab = X_lab, X_unlab = X_unlab,

      Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab)
  }

  return(out)
}

#=== END =======================================================================
