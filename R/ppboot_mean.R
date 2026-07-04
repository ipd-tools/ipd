#' Efron bootstrap PPI for mean estimation
#'
#' @description
#' `efron_boot_mean()` estimates a marginal mean using the Efron PPboot
#' bootstrap helpers. It dispatches to the binary PPboot1 helper when the
#' outcome is binary and the predictions are valid probabilities; otherwise it
#' dispatches to the continuous PPboot2-inspired helper.
#'
#' @details
#' The design matrix arguments are accepted for compatibility with the
#' matrix-based `ipd()` helper interface. Mean estimation uses an intercept-only
#' downstream statistic, so `X_l` and `X_u` are not otherwise used.
#'
#' When `binary = NULL`, the function treats the problem as binary if all
#' non-missing labeled outcomes are in `{0, 1}` and all labeled and unlabeled
#' predictions are strictly between 0 and 1.
#'
#' @param X_l Numeric design matrix for labeled observations. Accepted for
#'   interface compatibility; not used for mean estimation.
#' @param Y_l Numeric vector or one-column matrix of observed labeled outcomes.
#' @param f_l Numeric vector or one-column matrix of predictions for labeled
#'   observations.
#' @param X_u Numeric design matrix for unlabeled observations. Accepted for
#'   interface compatibility; not used for mean estimation.
#' @param f_u Numeric vector or one-column matrix of predictions for unlabeled
#'   observations.
#' @param nboot Integer number of bootstrap replicates. Default is 1000.
#' @param binary Optional logical value indicating whether to use the binary
#'   PPboot helper. If `NULL`, the function infers this from `Y_l`, `f_l`, and
#'   `f_u`.
#' @param seed Optional integer random seed.
#' @param return_boot Logical. If `TRUE`, return bootstrap replicate matrices
#'   from the selected PPboot helper.
#' @param ... Additional arguments passed to `efron_boot_logistic()` when
#'   `binary = TRUE` or to `efron_boot_ols()` when `binary = FALSE`.
#'
#' @return A list with the same components returned by
#'   `efron_boot_logistic()` or `efron_boot_ols()`. The `est`, `se`,
#'   `labeled_est`, `labeled_se`, and `classical_est` components are
#'   length-one named numeric vectors representing the marginal mean.
#'
#' @examples
#' set.seed(1)
#' n_l <- 80
#' n_u <- 120
#' f_l <- rnorm(n_l, mean = 1)
#' f_u <- rnorm(n_u, mean = 1)
#' Y_l <- 0.5 + 0.8 * f_l + rnorm(n_l)
#'
#' fit <- efron_boot_mean(
#'   X_l = matrix(1, n_l, 1), Y_l = Y_l, f_l = f_l,
#'   X_u = matrix(1, n_u, 1), f_u = f_u,
#'   nboot = 25, seed = 123
#' )
#'
#' fit$est
#' fit$se
#'
#' @export
efron_boot_mean <- function(
    X_l,
    Y_l,
    f_l,
    X_u,
    f_u,
    nboot = 1000,
    binary = NULL,
    seed = NULL,
    return_boot = FALSE,
    ...
) {
  if (is.null(binary)) {
    binary <- all(stats::na.omit(Y_l) %in% c(0, 1)) &&
      all(f_l > 0 & f_l < 1) &&
      all(f_u > 0 & f_u < 1)
  }

  mean_stat <- function(X, y) {
    c("(Intercept)" = mean(y))
  }

  if (binary) {
    out <- efron_boot_logistic(
      X_l = matrix(1, nrow = length(Y_l), ncol = 1,
                   dimnames = list(NULL, "(Intercept)")),
      Y_l = Y_l,
      f_l = f_l,
      X_u = matrix(1, nrow = length(f_u), ncol = 1,
                   dimnames = list(NULL, "(Intercept)")),
      f_u = f_u,
      nboot = nboot,
      statistic = mean_stat,
      seed = seed,
      return_boot = return_boot,
      ...
    )
  } else {
    out <- efron_boot_ols(
      X_l = matrix(1, nrow = length(Y_l), ncol = 1,
                   dimnames = list(NULL, "(Intercept)")),
      Y_l = Y_l,
      f_l = f_l,
      X_u = matrix(1, nrow = length(f_u), ncol = 1,
                   dimnames = list(NULL, "(Intercept)")),
      f_u = f_u,
      nboot = nboot,
      statistic = mean_stat,
      seed = seed,
      return_boot = return_boot,
      ...
    )
  }

  out
}
