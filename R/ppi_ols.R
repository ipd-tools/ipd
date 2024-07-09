#===============================================================================
# PPI ORDINARY LEAST SQUARES
#===============================================================================

#--- PPI OLS -------------------------------------------------------------------

#' PPI OLS
#'
#' @description
#' Helper function for prediction-powered inference for OLS estimation
#'
#' @details
#' Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://www.science.org/doi/10.1126/science.adi6000>
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
#' @param w_l (ndarray, optional): Sample weights for the labeled data set.
#' Defaults to a vector of ones.
#'
#' @param w_u (ndarray, optional): Sample weights for the unlabeled
#' data set. Defaults to a vector of ones.
#'
#' @returns A list of outputs: estimate of inference model parameters and corresponding standard errors
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
#' ppi_ols(X_l, Y_l, f_l, X_u, f_u)
#'
#' @import stats
#'
#' @export

ppi_ols <- function(X_l, Y_l, f_l, X_u, f_u, w_l = NULL, w_u = NULL) {

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
