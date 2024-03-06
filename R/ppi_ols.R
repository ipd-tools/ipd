#=== ORDINARY LEAST SQUARES ====================================================

#' PPI Linear Regression using Angelopoulos et al. (2023) Analytic Form
#'
#' @description
#' Computes the prediction-powered point estimate of the OLS coefficients.
#'
#' @details
#' Additional details...
#'
#' @param X_l (ndarray): Covariates corresponding to the gold-standard labels.
#'
#' @param Y_l (ndarray): Gold-standard labels.
#'
#' @param f_l (ndarray): Predictions corresponding to the gold-standard labels.
#'
#' @param X_u (ndarray): Covariates corresponding to the unlabeled data.
#'
#' @param f_u (ndarray): Predictions corresponding to the unlabeled data.
#'
#' @param lhat (float, optional): Power-tuning parameter (see `[ADZ23] <https://arxiv.org/abs/2311.01453>`__). The default value `None` will estimate the optimal value from data. Setting `lhat=1` recovers PPI with no power tuning, and setting `lhat=0` recovers the classical point estimate.
#'
#' @param coord (int, optional): Coordinate for which to optimize `lhat`. If `None`, it optimizes the total variance over all coordinates. Must be in {1, ..., d} where d is the shape of the estimand.
#'
#' @param w (ndarray, optional): Sample weights for the labeled data set.
#'
#' @param w_unlabeled (ndarray, optional): Sample weights for the unlabeled data set.
#'
#' @returns A list of outputs: estimate of inference model parameters and corresponding standard errors
#'
#' @examples
#'
#' dat <- simdat()
#'
#' form <- Y - Yhat ~ X1
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
#' ppi_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_ols <- function(X_l, Y_l, f_l, X_u, f_u,

  lhat = NULL, coord = NULL, w = NULL, w_unlabeled = NULL) {

  n <- NROW(X_l)

  p <- NCOL(X_l)

  N <- NROW(X_u)

  #- 1. Prediction-Powered Estimator

  theta_tilde_f <- solve(crossprod(X_u)) %*% t(X_u) %*% f_u

  delta_hat_f   <- solve(crossprod(X_l)) %*% t(X_l) %*% (f_l - Y_l)

  theta_hat_pp  <- theta_tilde_f - delta_hat_f

  #- 2. Meat and Bread for Imputed Estimate

  Sigma_tilde <- crossprod(X_u) / N

  M_tilde <- sapply(1:N, function(i) {

    (c(f_u[i] - crossprod(X_u[i,], theta_tilde_f)))^2 *

      tcrossprod(X_u[i,])}) |>

    rowMeans() |> matrix(nrow = p)

  iSigma_tilde <- solve(Sigma_tilde)

  #- 3. Sandwich Variance Estimator for Imputed Estimate

  V_tilde <- iSigma_tilde %*% M_tilde %*% iSigma_tilde

  #- 4. Meat and Bread for Empirical Rectifier

  Sigma <- crossprod(X_l) / n

  M <- sapply(1:n, function(i) {

    (c(f_l[i] - Y_l[i] - crossprod(X_l[i,], delta_hat_f)))^2 *

      tcrossprod(X_l[i,])}) |>

    rowMeans() |> matrix(nrow = p)

  iSigma <- solve(Sigma)

  #- 5. Sandwich Variance Estimator for Empirical Rectifier

  V <- iSigma %*% M %*% iSigma

  #- 6. Standard Error Estimates

  se <- sqrt(diag(V) / n + diag(V_tilde) / N)

  #- Output

  return(list(est = as.vector(theta_hat_pp), se = as.vector(se)))
}

#=== END =======================================================================
