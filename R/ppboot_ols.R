#' Efron bootstrap PPI for continuous outcomes
#'
#' @description
#' `efron_boot_ols()` implements a PPboot2-inspired bootstrap procedure for
#' continuous responses. The helper calibrates the relationship between the
#' observed labeled outcomes and their predictions, simulates bootstrap
#' outcomes for the unlabeled sample, and evaluates an OLS-style downstream
#' statistic on each bootstrap sample.
#'
#' @details
#' The calibration model estimates the conditional mean of `Y_l` given `f_l`
#' with a raw polynomial in `f_l`. When `degree = NULL`, the polynomial degree is
#' selected by AIC over degrees 1 through `max_degree`. A smooth conditional
#' scale estimate is then fit as a function of `f_l`, standardized residuals are
#' sampled with replacement, and bootstrap outcomes are generated for both the
#' unlabeled and labeled design matrices.
#'
#' @param X_l Numeric design matrix for labeled observations. Rows must
#'   correspond to `Y_l` and `f_l`.
#' @param Y_l Numeric vector or one-column matrix of observed continuous
#'   outcomes for labeled observations.
#' @param f_l Numeric vector or one-column matrix of predicted outcomes for
#'   labeled observations.
#' @param X_u Numeric design matrix for unlabeled observations. Rows must
#'   correspond to `f_u`.
#' @param f_u Numeric vector or one-column matrix of predicted outcomes for
#'   unlabeled observations.
#' @param nboot Integer number of bootstrap replicates. Default is 1000.
#' @param degree Optional integer polynomial degree for the conditional mean
#'   calibration model. If `NULL`, the degree is selected by AIC.
#' @param max_degree Maximum polynomial degree considered when `degree = NULL`.
#'   Default is 5.
#' @param conditional Logical. If `TRUE`, use model-based residual resampling.
#'   If `FALSE`, nonparametrically resample labeled rows before refitting the
#'   calibration model.
#' @param statistic Optional function of the form `function(X, y)` returning a
#'   numeric vector. Defaults to OLS coefficients from `stats::lm.fit()`.
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
#'   \item{calibration_model}{List containing the fitted mean and scale
#'   calibration objects and prediction functions.}
#'   \item{degree}{Integer polynomial degree used in the mean calibration
#'   model.}
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
#' f_l <- 0.5 + 1.2 * x_l + rnorm(n_l, sd = 0.5)
#' f_u <- 0.5 + 1.2 * x_u + rnorm(n_u, sd = 0.5)
#' Y_l <- 0.2 + 0.8 * x_l + 0.7 * f_l + rnorm(n_l)
#'
#' fit <- efron_boot_ols(
#'   X_l = X_l, Y_l = Y_l, f_l = f_l,
#'   X_u = X_u, f_u = f_u,
#'   nboot = 25, seed = 123
#' )
#'
#' fit$est
#' fit$se
#'
#' @export
efron_boot_ols <- function(
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    nboot = 1000,
    degree = NULL,
    max_degree = 5,
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
    statistic <- .default_ols_statistic
  }

  Y_l <- as.numeric(Y_l)
  f_l <- as.numeric(f_l)
  f_u <- as.numeric(f_u)

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

  model0 <- .fit_ppboot2_model(
    y = Y_l,
    f = f_l,
    degree = degree,
    max_degree = max_degree
  )

  mu_l <- model0$predict_mu(f_l)
  sd_l <- model0$predict_sd(f_l)

  stat0 <- statistic(X_l, Y_l)
  p <- length(stat0)

  boot_u <- matrix(NA_real_, nrow = nboot, ncol = p)
  boot_l <- matrix(NA_real_, nrow = nboot, ncol = p)

  colnames(boot_u) <- colnames(boot_l) <- names(stat0)

  for (b in seq_len(nboot)) {
    if (conditional) {
      eps_l <- sample(model0$eps, size = n_l, replace = TRUE)
      y_l_star <- mu_l + sd_l * eps_l
      f_l_star <- f_l
    } else {
      idx <- sample(seq_len(n_l), size = n_l, replace = TRUE)
      y_l_star <- Y_l[idx]
      f_l_star <- f_l[idx]
    }

    model_star <- .fit_ppboot2_model(
      y = y_l_star,
      f = f_l_star,
      degree = model0$degree,
      max_degree = max_degree
    )

    mu_u_star <- model_star$predict_mu(f_u)
    sd_u_star <- model_star$predict_sd(f_u)

    mu_l_star <- model_star$predict_mu(f_l)
    sd_l_star <- model_star$predict_sd(f_l)

    eps_u <- sample(model_star$eps, size = n_u, replace = TRUE)
    eps_l2 <- sample(model_star$eps, size = n_l, replace = TRUE)

    y_u_star <- mu_u_star + sd_u_star * eps_u
    y_l2_star <- mu_l_star + sd_l_star * eps_l2

    boot_u[b, ] <- statistic(X_u, y_u_star)
    boot_l[b, ] <- statistic(X_l, y_l2_star)
  }

  est <- colMeans(boot_u, na.rm = TRUE)
  se <- apply(boot_u, 2, stats::sd, na.rm = TRUE)

  labeled_est <- colMeans(boot_l, na.rm = TRUE)
  labeled_se <- apply(boot_l, 2, stats::sd, na.rm = TRUE)

  classical_est <- statistic(X_l, Y_l)

  list(
    est = est,
    se = se,
    labeled_est = labeled_est,
    labeled_se = labeled_se,
    classical_est = classical_est,
    calibration_model = model0,
    degree = model0$degree,
    conditional = conditional,
    nboot = nboot,
    boot_unlabeled = if (return_boot) boot_u else NULL,
    boot_labeled = if (return_boot) boot_l else NULL
  )
}
