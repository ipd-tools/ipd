#===============================================================================
#
#  PROGRAM: popinf_ols.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Implementation of algorithm from Miao et al. (2023)
#
#           Assumption-lean and Data-adaptive Post-Prediction Inference
#
#  INPUTS:  None
#
#  OUTPUS:   Assumption-lean and Data-adaptive Post-prediction inference
#
#            function (popinf_ols)
#
#  Notes:   1. Typo in Sigma in paper? Should be A^{-1} V A^{-1}?
#
#  Updated: 2023-12-12
#
#===============================================================================

#=== MIAO ET AL. (2023) ========================================================

#' Assumption-lean and data-adaptive post-prediction inference (Miao et al. 2023)
#'
#' @description
#' A short description...
#'
#' @details
#' Additional details...
#'
#' @param rec_form A formula defining the rectifier model. This should be of
#' the form Y - Yhat ~ X, where Y is the name of the column corresponding to
#' the observed outcome in the labeled data, Yhat is the name of the column
#' corresponding to the predicted outcome in the labeled data, and X generally
#' corresponds to the features of interest (e.g., X1 + X2).
#'
#' @param inf_form A formula defining the inference model. This should be of
#' the form Yhat ~ X, where Yhat is the name of the column corresponding to the
#' predicted outcome in the unlabeled data, and X generally corresponds to the
#' features of interest (e.g., X1 + X2).
#'
#' @param dat data in the form of the simdat function
#'
#' @param alpha scalar type I error rate for hypothesis testing - values in (0, 1); defaults to 0.05
#'
#' @param eta scalar step size for gradient descent; defaults to 0.001
#'
#' @param delta scalar tolerance for assessing convergence; defaults to 0.0001
#'
#' @param K integer maximum number of iterations; defaults to 100
#'
#' @returns A list of outputs: estimate of inference model parameters and corresponding standard error.
#'
#' @examples
#'
#' rec_form <- Y - Yhat ~ X1
#'
#' inf_form <- Yhat ~ X1
#'
#' dat <- simdat()
#'
#' popinf_ols(rec_form, inf_form, dat = dat)
#'
#' @export
#'
#' @import stats

popinf_ols <- function(rec_form, inf_form, dat,

                       alpha = 0.05, eta = 0.001, delta = 0.0001, K = 100) {

  #- 0. Inputs

  X_l <- model.matrix(rec_form, data = dat[dat$set == "tst",])

  Y_l <- dat[dat$set == "tst", all.vars(rec_form)[1]]

  f_l <- dat[dat$set == "tst", all.vars(rec_form)[2]]

  X_u <- model.matrix(inf_form, data = dat[dat$set == "val",])

  f_u <- dat[dat$set == "val", all.vars(inf_form)[1]]

  #- 1. Pre-Compute Common Values

  p <- ncol(X_l)

  n <- nrow(X_l)

  N <- nrow(X_u)

  XtX_l <- crossprod(X_l)

  XtX_u <- crossprod(X_u)

  XtY_l <- crossprod(X_l, Y_l)

  Xtf_l <- crossprod(X_l, f_l)

  Xtf_u <- crossprod(X_u, f_u)

  A <- XtX_l / n                                                                 # Or stack X_l and X_u?

  iA <- solve(A)

  #- 2. Initial Estimates

  rho <- n / N

  theta_0 <- theta <- solve(XtX_l) %*% XtY_l

  w_0 <- w <- rep(0, p)

  w_constr <- list(lb = 0, ub = 1)

  #- 3. Iteratively Update Parameter Estimates and Weights

  converge <- F

  k <- 1

  while ((!converge) & (k < K)) {

    #- a. Update Parameters

    theta_0 <- theta

    theta <- solve(XtX_l + w*(XtX_u - XtX_l)) %*% (XtY_l + w*(Xtf_u - Xtf_l))

    M1 <- sapply(1:n, function(i) {

      (c(Y_l[i] - crossprod(X_l[i,], theta)))^2 * tcrossprod(X_l[i,])}) |>

      rowMeans() |> matrix(nrow = ncol(X_l))

    M2 <- sapply(1:n, function(i) {

      (c(f_l[i] - crossprod(X_l[i,], theta)))^2 * tcrossprod(X_l[i,])}) |>

      rowMeans() |> matrix(nrow = ncol(X_l))

    M3 <- sapply(1:N, function(i) {

      (c(f_u[i] - crossprod(X_u[i,], theta)))^2 * tcrossprod(X_u[i,])}) |>

      rowMeans() |> matrix(nrow = ncol(X_u))

    M4 <- sapply(1:n, function(i) {                                              # Check this

      c(f_l[i] - crossprod(X_l[i,], theta)) *

        c(Y_l[i] - crossprod(X_l[i,], theta)) * tcrossprod(X_u[i,])}) |>

      rowMeans() |> matrix(nrow = ncol(X_l))

    V <- M1 + tcrossprod(w) * (M2 + rho * M3) - 2 * diag(w) * M4

    Sigma <- iA %*% V %*% iA

    #- b. Update Weights

    w_0 <- w

    w <- sapply(1:p, function(j) {

      w_j <- optim(w_0[j],

                   fn = function(x) {

                     return(M1[j,j] + (M2[j,j] + rho*M3[j,j])*x^2 - 2*M4[j,j]*x)
                   },

                   control = list(fnscale = -1), lower = 0, upper = 1,

                   method = "Brent")$par

      cat("j:", w_j)

      return(w_j)
    })

    #- c. Check Convergence

    if (all(abs(theta - theta_0) < delta)) {

      converge <- T

      k <- K

    } else {

      k <- k + 1
    }
  }

  return(list(est = theta, se = sqrt(diag(Sigma) / (n - p - 1)), w = w))
}

#=== QC ========================================================================

# simdat <- function(n = c(300, 300, 300), beta1 = 1) {
#
#   X1 <- rnorm(sum(n), 1)
#   X2 <- rnorm(sum(n), 1)
#   X3 <- rnorm(sum(n), 1)
#   X4 <- rnorm(sum(n), 2)
#
#   Y <- c(beta1*X1 + 0.5*X2 + 3*smooth(X3) + 4*smooth(X4) + rnorm(sum(n)))
#
#   set <- rep(c("trn", "tst", "val"), n)
#
#   dat <- data.frame(X1, X2, X3, X4, Y, Yhat = NA, set)
#
#   fit_gam <- gam(Y ~ s(X1) + s(X2) + s(X3) + s(X4), data = dat[set == "trn",])
#
#   dat[set == "tst", "Yhat"] <- predict(fit_gam, newdat = dat[set == "tst",])
#
#   dat[set == "val", "Yhat"] <- predict(fit_gam, newdat = dat[set == "val",])
#
#   return(dat)
# }
#
# rec_form <- Y - Yhat ~ X1
#
# inf_form <- Yhat ~ X1
#
# dat <- simdat()
#
# popinf_ols(rec_form, inf_form, dat = dat)

#=== END =======================================================================
