#' Efron bootstrap PPI for binary outcomes
#'
#' @description
#' `efron_boot_logistic()` implements the PPboot1 algorithm from Efron (2025)
#' for binary responses. The prediction rule is calibrated on the labeled data,
#' bootstrap outcomes are generated from calibrated probabilities, and a
#' downstream logistic regression statistic is evaluated on each bootstrap
#' sample.
#'
#' @details
#' The calibration model is a binomial GLM using a raw polynomial in
#' `logit(f_l)` as the calibration design. The default
#' `calibration_degree = 1` matches the main PPboot1 calibration model. The
#' default downstream statistic is the coefficient vector from
#' `stats::glm.fit(x = X, y = y, family = stats::binomial())`.
#'
#' Predicted probabilities in `f_l` and `f_u` are clipped internally to avoid
#' infinite logits.
#'
#' @param X_l Numeric design matrix for labeled observations. Rows must
#'   correspond to `Y_l` and `f_l`.
#' @param Y_l Numeric vector or one-column matrix of observed binary outcomes
#'   for labeled observations. Outcomes must be coded as 0 and 1.
#' @param f_l Numeric vector or one-column matrix of predicted probabilities
#'   for labeled observations.
#' @param X_u Numeric design matrix for unlabeled observations. Rows must
#'   correspond to `f_u`.
#' @param f_u Numeric vector or one-column matrix of predicted probabilities
#'   for unlabeled observations.
#' @param nboot Integer number of bootstrap replicates. Default is 1000.
#' @param calibration_degree Integer degree of the raw polynomial in
#'   `logit(f)` used by the calibration model. Default is 1.
#' @param conditional Logical. If `TRUE`, use conditional parametric resampling
#'   from the calibrated probabilities. If `FALSE`, resample labeled outcomes
#'   nonparametrically in the first bootstrap level.
#' @param statistic Optional function of the form `function(X, y)` returning a
#'   numeric vector. Defaults to logistic regression coefficients.
#' @param seed Optional integer random seed.
#' @param return_boot Logical. If `TRUE`, return bootstrap replicate matrices.
#' @param ... Additional arguments accepted for compatibility with the `ipd()`
#'   wrapper. They are currently unused.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{est}{Named numeric vector of bootstrap estimates computed on
#'   generated unlabeled outcomes.}
#'   \item{se}{Named numeric vector of bootstrap standard errors computed from
#'   generated unlabeled outcomes.}
#'   \item{labeled_est}{Named numeric vector of bootstrap estimates computed on
#'   generated labeled outcomes.}
#'   \item{labeled_se}{Named numeric vector of bootstrap standard errors
#'   computed from generated labeled outcomes.}
#'   \item{classical_est}{Named numeric vector of the downstream statistic
#'   computed from the observed labeled data.}
#'   \item{calibration_fit}{The initial binomial calibration fit returned by
#'   `stats::glm.fit()`.}
#'   \item{calibration_degree}{Integer calibration polynomial degree.}
#'   \item{conditional}{Logical value indicating whether conditional resampling
#'   was used.}
#'   \item{nboot}{Integer number of bootstrap replicates.}
#'   \item{boot_unlabeled}{Bootstrap replicate matrix for unlabeled outcomes
#'   when `return_boot = TRUE`; otherwise `NULL`.}
#'   \item{boot_labeled}{Bootstrap replicate matrix for labeled outcomes when
#'   `return_boot = TRUE`; otherwise `NULL`.}
#' }
#'
#' @examples
#' set.seed(1)
#' n_l <- 80
#' n_u <- 120
#' x_l <- rnorm(n_l)
#' x_u <- rnorm(n_u)
#'
#' X_l <- cbind("(Intercept)" = 1, X1 = x_l)
#' X_u <- cbind("(Intercept)" = 1, X1 = x_u)
#' f_l <- plogis(-0.1 + 0.9 * x_l + rnorm(n_l, sd = 0.3))
#' f_u <- plogis(-0.1 + 0.9 * x_u + rnorm(n_u, sd = 0.3))
#' Y_l <- rbinom(n_l, size = 1, prob = plogis(-0.3 + 1.1 * x_l))
#'
#' fit <- efron_boot_logistic(
#'   X_l = X_l, Y_l = Y_l, f_l = f_l,
#'   X_u = X_u, f_u = f_u,
#'   nboot = 25, seed = 123
#' )
#'
#' fit$est
#' fit$se
#'
#' @export
efron_boot_logistic <- function(
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    nboot = 1000,
    calibration_degree = 1,
    conditional = TRUE,
    statistic = NULL,
    seed = NULL,
    return_boot = FALSE,
    ...
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(statistic)) {
    statistic <- .default_logistic_statistic
  }

  Y_l <- as.numeric(Y_l)

  if (!all(Y_l %in% c(0, 1))) {
    stop("`efron_boot_logistic()` requires binary outcomes coded as 0/1.",
         call. = FALSE)
  }

  n_l <- length(Y_l)
  n_u <- length(f_u)

  if (n_l != nrow(X_l)) {
    stop("Length of `Y_l` must equal `nrow(X_l)`.", call. = FALSE)
  }

  if (length(f_l) != n_l) {
    stop("Length of `f_l` must equal length of `Y_l`.", call. = FALSE)
  }

  if (length(f_u) != n_u || n_u != nrow(X_u)) {
    stop("Length of `f_u` must equal `nrow(X_u)`.", call. = FALSE)
  }

  L_l <- .make_logit_calibration_design(f_l, degree = calibration_degree)
  L_u <- .make_logit_calibration_design(f_u, degree = calibration_degree)

  # Initial calibration model: y_l ~ logit(f_l)
  cal0 <- .fit_binomial_calibration(L_l, Y_l)
  pi_l_hat <- cal0$fitted

  # Determine statistic dimension and names.
  stat0 <- statistic(X_l, Y_l)
  p <- length(stat0)

  boot_u <- matrix(NA_real_, nrow = nboot, ncol = p)
  boot_l <- matrix(NA_real_, nrow = nboot, ncol = p)

  colnames(boot_u) <- colnames(boot_l) <- names(stat0)

  for (b in seq_len(nboot)) {
    # Step 1: bootstrap labeled outcomes.
    if (conditional) {
      y_l_star <- stats::rbinom(n_l, size = 1, prob = pi_l_hat)
    } else {
      y_l_star <- sample(Y_l, size = n_l, replace = TRUE)
    }

    # Step 2: refit calibration model on labeled bootstrap data.
    cal_star <- .fit_binomial_calibration(L_l, y_l_star)

    # Step 3: map calibration relationship to unlabeled f-values.
    pi_u_star <- .predict_binomial_calibration(L_u, cal_star$beta)
    pi_l_star <- .predict_binomial_calibration(L_l, cal_star$beta)

    # Step 4: generate bootstrap outcomes.
    y_u_star <- stats::rbinom(n_u, size = 1, prob = pi_u_star)
    y_l2_star <- stats::rbinom(n_l, size = 1, prob = pi_l_star)

    # Step 5: compute downstream statistic.
    boot_u[b, ] <- statistic(X_u, y_u_star)
    boot_l[b, ] <- statistic(X_l, y_l2_star)
  }

  est <- colMeans(boot_u, na.rm = TRUE)
  se <- apply(boot_u, 2, stats::sd, na.rm = TRUE)

  labeled_est <- colMeans(boot_l, na.rm = TRUE)
  labeled_se <- apply(boot_l, 2, stats::sd, na.rm = TRUE)

  classical_est <- statistic(X_l, Y_l)

  out <- list(
    est = est,
    se = se,
    labeled_est = labeled_est,
    labeled_se = labeled_se,
    classical_est = classical_est,
    calibration_fit = cal0$fit,
    calibration_degree = calibration_degree,
    conditional = conditional,
    nboot = nboot,
    boot_unlabeled = if (return_boot) boot_u else NULL,
    boot_labeled = if (return_boot) boot_l else NULL
  )

  out
}


