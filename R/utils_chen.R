# utils_chen.R
#===============================================================================
# Chen-Chen, PDC, and PPI-PP HELPER FUNCTIONS
#===============================================================================

#--- EXPIT / INVERSE-LOGIT -----------------------------------------------------

#' expit / inverse-logit
#'
#' Computes logistic inverse link \eqn{\mathrm{expit}(x) = 1 / (1 + e^{-x})}.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of same length as `x`.
#'
#' @keywords internal
#'
#' @noRd
expit <- function(x) {

  1 / (1 + exp(-x))
}

#--- DERIVATIVE OF EXPIT -------------------------------------------------------

#' Derivative of expit / inverse-logit
#'
#' Computes \eqn{\mathrm{expit}(x)(1-\mathrm{expit}(x))}.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of same length as `x`.
#'
#' @keywords internal
#'
#' @noRd
expit_derivative <- function(x) {

  p <- expit(x)

  p * (1 - p)
}

#--- RESIDUALS + HESSIAN FOR GLM SCORE FUNCTIONS -------------------------------

#' Compute GLM residuals and (average) Hessian for canonical link families
#'
#' Internal helper used by Chen-Chen / PDC / PPI-PP estimators.
#'
#' @param Y Numeric outcome vector of length n.
#'
#' @param X_int Numeric design matrix including intercept column, n x p.
#'
#' @param theta Numeric coefficient vector of length p.
#'
#' @param family Character string specifying the family. Supported:
#'   `"gaussian"`, `"binomial"`, `"poisson"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{residuals}{Numeric vector of residuals, length n.}
#'   \item{hessian}{Numeric matrix, the average Hessian, p x p.}
#' }
#'
#' @details
#'
#' This function assumes canonical links:
#'
#' \itemize{
#'  \item Gaussian: identity, \eqn{\mu = \eta}
#'  \item Binomial: logit, \eqn{\mu = \mathrm{expit}(\eta)}
#'  \item Poisson: log link, \eqn{\mu = \exp(\eta)}
#' }
#'
#' The returned Hessian matrix is:
#'
#' \deqn{H = -(1/n) X^\top W X}
#'
#' where `W` is `d mu / d eta` for the binomial/poisson cases and `I` for
#' the Gaussian case.
#'
#' @keywords internal
#'
#' @noRd
compute_residuals <- function(Y, X_int, theta, family) {

  if (!family %in% c("gaussian", "binomial", "poisson")) {

    stop("Only Gaussian, binomial, and Poisson families are supported.",
         call. = FALSE)
  }

  eta_hat <- X_int %*% theta

  samp_size <- nrow(X_int)

  if (family == "gaussian") {

    residuals <- c(Y - eta_hat)

    hessian <- -crossprod(X_int) / samp_size

  } else if (family == "binomial") {

    mu_hat <- expit(eta_hat)

    residuals <- c(Y - mu_hat)

    expit_deriv <- c(expit_derivative(eta_hat))

    hessian <- -t(X_int) %*% (X_int * expit_deriv) / samp_size

  } else { # poisson

    mu_hat <- exp(eta_hat)

    residuals <- c(Y - mu_hat)

    exp_deriv <- c(exp(eta_hat))

    hessian <- -t(X_int) %*% (X_int * exp_deriv) / samp_size
  }

  list(residuals = residuals, hessian = hessian)
}

#--- INTERNAL FAMILY RESOLVER ---------------------------------------------------

#' Resolve a GLM family argument to a stats family object
#'
#' @param family Character ("gaussian","binomial","poisson"), a `family` object,
#'   or a family function (e.g., `gaussian`, `binomial`, `poisson`).
#'
#' @return A stats `family` object.
#'
#' @keywords internal
#'
#' @noRd
.resolve_family <- function(family) {

  if (inherits(family, "family")) return(family)

  if (is.function(family)) return(family())

  if (is.character(family) && length(family) == 1L) {

    return(switch(
      family,
      gaussian = stats::gaussian(),
      binomial = stats::binomial(),
      poisson  = stats::poisson(),
      stop("Unsupported family: '", family, "'.", call. = FALSE)
    ))
  }

  stop("`family` must be a character string, a family object, or a family function.", call. = FALSE)
}

#--- INTERNAL GLM MATRIX FITTER ------------------------------------------------

#' Fit a GLM using a design matrix
#'
#' Convenience helper to fit GLMs directly on design matrices so that the
#' method helpers can remain matrix-based.
#'
#' @param y Numeric response vector.
#' @param X Numeric design matrix.
#' @param family Character family ("gaussian"/"binomial"/"poisson") or `family` object.
#'
#' @return Numeric coefficient vector.
#'
#' @keywords internal
#' @noRd
.fit_glm_from_matrix <- function(y, X, family) {
  fam <- .resolve_family(family)
  stats::glm.fit(x = X, y = as.numeric(y), family = fam)$coefficients
}

#--- CHEN-CHEN MIN-MSE AUGMENTED ESTIMATOR -------------------------------------

#' Chen-Chen (min-MSE) Estimator
#'
#' Implements Chen-Chen / "min-MSE" augmented estimator.
#'
#' @param Y_lab Labeled outcomes, length n.
#'
#' @param pred_lab Labeled predictions, length n.
#'
#' @param X_lab_int Labeled design matrix (with intercept), n x p.
#'
#' @param X_all_int All-data design matrix (with intercept), N x p.
#'
#' @param beta_lab Coefs of target model on labeled data, length p.
#'
#' @param gamma_lab Coefs of pred~X on labeled data, length p.
#'
#' @param gamma_all Coefs of pred~X on all data, length p.
#'
#' @param family Character family for both target and auxiliary model.
#'  Supported: `"gaussian"`, `"binomial"`, `"poisson"`.
#'
#' @return A list with:
#' \describe{
#'   \item{theta_hat}{Corrected coefficient vector, length p.}
#'   \item{theta_hat_cov}{Estimated covariance matrix, p x p.}
#' }
#'
#' @keywords internal
#'
#' @noRd
compute_min_mse_est <- function(
    Y_lab,
    pred_lab,
    X_lab_int,
    X_all_int,
    beta_lab,
    gamma_lab,
    gamma_all,
    family) {

  n <- nrow(X_lab_int)
  N <- nrow(X_all_int)

  rho <- n / N

  beta_res  <- compute_residuals(Y_lab,    X_lab_int, beta_lab,  family)
  gamma_res <- compute_residuals(pred_lab, X_lab_int, gamma_all, family)

  residuals_beta  <- beta_res$residuals
  D1 <- beta_res$hessian

  residuals_gamma <- gamma_res$residuals
  D2 <- gamma_res$hessian

  score_beta  <- X_lab_int * residuals_beta
  score_gamma <- X_lab_int * residuals_gamma

  C11 <- crossprod(score_beta) / n
  C12 <- crossprod(score_beta, score_gamma) / n
  C22 <- crossprod(score_gamma) / n

  D1_inv  <- solve(D1)
  C22_inv <- solve(C22)

  D1C12 <- D1_inv %*% C12
  D1C12C22 <- D1C12 %*% C22_inv

  theta_hat <- as.vector(beta_lab - D1C12C22 %*% D2 %*% (gamma_lab - gamma_all))

  theta_hat_cov <- (D1_inv %*% C11 %*% D1_inv -
                      (1 - rho) * D1C12C22 %*% t(D1C12)) / n

  list(theta_hat = theta_hat, theta_hat_cov = theta_hat_cov)
}

#--- PDC AUGMENTED ESTIMATOR ---------------------------------------------------

#' PDC Estimator
#'
#' Implements the Prediction De-Correlated Inference (PDC) estimator. The
#' augmentation weights are obtained by regressing each component of the
#' target score onto an augmented basis of prediction score features, and
#' then applying the mean shift between labeled and all data.
#'
#' @param Y_lab Labeled outcomes, length n.
#'
#' @param pred_lab Labeled predictions, length n.
#'
#' @param pred_all Predictions on all data, length N.
#'
#' @param X_lab_int Labeled design matrix (with intercept), n x p.
#'
#' @param X_all_int All-data design matrix (with intercept), N x p.
#'
#' @param beta_lab Coefs of target model on labeled data, length p.
#'
#' @param family Character family.
#'  Supported: `"gaussian"`, `"binomial"`, `"poisson"`.
#'
#' @return A list with:
#' \describe{
#'   \item{theta_hat}{Corrected coefficient vector, length p.}
#'   \item{theta_hat_cov}{Estimated covariance matrix, p x p.}
#' }
#'
#' @keywords internal
#'
#' @noRd
compute_pdc_est <- function(
    Y_lab,
    pred_lab,
    pred_all,
    X_lab_int,
    X_all_int,
    beta_lab,
    family) {

  n <- nrow(X_lab_int)
  N <- nrow(X_all_int)

  rho <- n / N

  beta_res <- compute_residuals(Y_lab, X_lab_int, beta_lab, family)

  gamma_res     <- compute_residuals(pred_lab, X_lab_int, beta_lab, family)
  gamma_res_all <- compute_residuals(pred_all, X_all_int, beta_lab, family)

  residuals_beta <- beta_res$residuals
  D1 <- beta_res$hessian

  residuals_gamma <- gamma_res$residuals

  score_beta <- X_lab_int * residuals_beta

  score_gamma_all <- X_all_int * gamma_res_all$residuals
  score_gamma_all_int <- cbind(1, score_gamma_all)

  score_gamma <- X_lab_int * residuals_gamma
  score_gamma_int <- cbind(1, score_gamma)

  weights <- matrix(NA_real_,
                    nrow = ncol(score_beta), ncol = ncol(score_gamma_int))

  for (i in seq_len(ncol(score_beta))) {

    s <- score_beta[, i]

    weights[i, ] <- stats::glm(s ~ 0 + score_gamma_int)$coefficients
  }

  diff <- matrix(colMeans(score_gamma_int) - colMeans(score_gamma_all_int),
                 ncol(score_gamma_all_int), 1)

  D1_inv <- solve(D1)

  augmentation <- -D1_inv %*% weights %*% diff

  theta_hat <- as.vector(beta_lab - augmentation)

  C11 <- crossprod(score_beta) / n
  C12 <- crossprod(score_beta, score_gamma) / n
  C22 <- crossprod(score_gamma) / n

  C22_inv <- solve(C22)
  D1C12 <- D1_inv %*% C12
  D1C12C22 <- D1C12 %*% C22_inv

  theta_hat_cov <- (D1_inv %*% C11 %*% D1_inv -
                      (1 - rho) * D1C12C22 %*% t(D1C12)) / n

  list(theta_hat = theta_hat, theta_hat_cov = theta_hat_cov)
}

#--- NEW CHEN-CHEN MIN-MSE (SEPARATE AUXILIARY DESIGN/FAMILY) ------------------

#' Chen-Chen (min-MSE) estimator with separate auxiliary design/family           ### THIS IS NOT YET IMPEMENTED
#'
#' This is your "new CC" core estimator: target model uses `X_*_int` and
#' auxiliary prediction model uses `X_*_int_cc` and `family_cc`.
#'
#' @param Y_lab Labeled outcomes, length n.
#' @param pred_lab_cc Labeled predictions to use for CC auxiliary model, length n.
#' @param X_lab_int Target labeled design matrix (with intercept), n x p.
#' @param X_all_int Target all-data design matrix (with intercept), N x p.
#' @param X_lab_int_cc Auxiliary labeled design matrix (with intercept), n x q.
#' @param X_all_int_cc Auxiliary all-data design matrix (with intercept), N x q.
#' @param beta_lab Target coefs on labeled, length p.
#' @param gamma_lab Auxiliary coefs on labeled, length q.
#' @param gamma_all Auxiliary coefs on all, length q.
#' @param family Target family (character).
#' @param family_cc Auxiliary family (character).
#'
#' @return A list with `theta_hat` and `theta_hat_cov`.
#'
#' @keywords internal
#' @noRd
compute_min_mse_est_new <- function(Y_lab, pred_lab_cc,
                                    X_lab_int, X_all_int,
                                    X_lab_int_cc, X_all_int_cc,
                                    beta_lab, gamma_lab, gamma_all,
                                    family, family_cc) {

  n <- nrow(X_lab_int)
  N <- nrow(X_all_int)
  rho <- n / N

  beta_res <- compute_residuals(Y_lab, X_lab_int, beta_lab, family)

  gamma_res <- compute_residuals(pred_lab_cc, X_lab_int_cc, gamma_all, family_cc)

  residuals_beta <- beta_res$residuals
  D1 <- beta_res$hessian

  residuals_gamma <- gamma_res$residuals
  D2 <- gamma_res$hessian

  score_beta <- X_lab_int * residuals_beta
  score_gamma <- X_lab_int_cc * residuals_gamma

  C11 <- crossprod(score_beta) / n
  C12 <- crossprod(score_beta, score_gamma) / n
  C22 <- crossprod(score_gamma) / n

  D1_inv <- solve(D1)
  C22_inv <- solve(C22)

  D1C12 <- D1_inv %*% C12
  D1C12C22 <- D1C12 %*% C22_inv

  theta_hat <- as.vector(beta_lab - D1C12C22 %*% D2 %*% (gamma_lab - gamma_all))

  theta_hat_cov <- (D1_inv %*% C11 %*% D1_inv -
                      (1 - rho) * D1C12C22 %*% t(D1C12)) / n

  list(theta_hat = theta_hat, theta_hat_cov = theta_hat_cov)
}
