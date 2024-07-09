#===============================================================================
# HELPERS
#===============================================================================

#=== PPI++ =====================================================================

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
#' @returns (list): A list containing the following:
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

ols <- function(X, Y, return_se = F) {

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
#' @returns (list): A list containing the following:
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

wls <- function(X, Y, w = NULL, return_se = F) {

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
#' @returns (list): A list containing the following:
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
#' dat <- simdat()
#'
#' form <- Y - f ~ X1
#'
#' X_l <- model.matrix(form, data = dat[dat$set == "labeled",])
#'
#' Y_l <- dat[dat$set == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set == "unlabeled",])
#'
#' f_u <- dat[dat$set == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
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
#' Must be in {1, ..., d} where d is the shape of the estimand.
#'
#' @param clip (boolean, optional): Whether to clip the value of lhat to be
#' non-negative. Defaults to `False`.
#'
#' @examples
#'
#' dat <- simdat()
#'
#' form <- Y - f ~ X1
#'
#' X_l <- model.matrix(form, data = dat[dat$set == "labeled",])
#'
#' Y_l <- dat[dat$set == "labeled", all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set == "labeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set == "unlabeled",])
#'
#' f_u <- dat[dat$set == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' est <- ppi_plusplus_ols_est(X_l, Y_l, f_l, X_u, f_u)
#'
#' stats <- ols_get_stats(est, X_l, Y_l, f_l, X_u, f_u)
#'
#' calc_lhat_glm(stats$grads, stats$grads_hat, stats$grads_hat_unlabeled,
#'
#'   stats$inv_hessian, coord = NULL, clip = FALSE)
#'
#' @returns (float): Optimal value of `lhat` in \[0,1\].
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

# """generic (normal) z-test based on summary statistic
#
#     The test statistic is :
#         tstat = (value1 - value2 - diff) / std_diff
#
#     and is assumed to be normally distributed.
#
#     Parameters
#     ----------
#     value1 : float or ndarray
#         Value, for example mean, of the first sample.
#     value2 : float or ndarray
#         Value, for example mean, of the second sample.
#     std_diff : float or ndarray
#         Standard error of the difference value1 - value2
#     alternative : str
#         The alternative hypothesis, H1, has to be one of the following
#
#            * 'two-sided' : H1: ``value1 - value2 - diff`` not equal to 0.
#            * 'larger' :   H1: ``value1 - value2 - diff > 0``
#            * 'smaller' :  H1: ``value1 - value2 - diff < 0``
#
#     diff : float
#         value of difference ``value1 - value2`` under the null hypothesis
#
#     Returns
#     -------
#     tstat : float or ndarray
#         Test statistic.
#     pvalue : float or ndarray
#         P-value of the hypothesis test assuming that the test statistic is
#         t-distributed with ``df`` degrees of freedom.
#     """

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
#' @returns (vector): Lower and upper (1 - alpha) * 100% confidence limits.
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


### NEED TO DOCUMENT ###########################################################

log1pexp <- function(x) {

  idxs <- x > 10

  out <- numeric(length(x))

  out[idxs] <- x[idxs]

  out[!idxs] <- log1p(exp(x[!idxs]))

  return(out)
}

#=== END =======================================================================
