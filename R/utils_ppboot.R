#' Clip probabilities away from zero and one
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small positive value used as the lower and upper clipping margin.
#'
#' @return Numeric vector with all values in `[eps, 1 - eps]`.
#'
#' @keywords internal
#' @noRd
.clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

#' Build a logit-polynomial calibration design
#'
#' @param f Numeric vector of predicted probabilities.
#' @param degree Integer degree of the raw polynomial in `logit(f)`.
#'
#' @return Numeric design matrix with an intercept column and one column per
#'   polynomial degree.
#'
#' @keywords internal
#' @noRd
.make_logit_calibration_design <- function(f, degree = 1) {
  f <- .clip_prob(f)
  eta <- qlogis(f)

  if (degree < 1) {
    stop("`degree` must be >= 1.", call. = FALSE)
  }

  out <- cbind("(Intercept)" = 1, sapply(seq_len(degree), function(j) eta^j))
  colnames(out) <- c("(Intercept)", paste0("logit_f_", seq_len(degree)))
  out
}

#' Fit a binomial calibration model
#'
#' @param L Numeric calibration design matrix.
#' @param y Numeric binary outcome vector coded as 0 and 1.
#'
#' @return A list with `beta`, the fitted calibration coefficients; `fitted`,
#'   the clipped fitted probabilities; and `fit`, the raw `glm.fit` object.
#'
#' @keywords internal
#' @noRd
.fit_binomial_calibration <- function(L, y) {
  fit <- suppressWarnings(
    stats::glm.fit(x = L, y = y, family = stats::binomial())
  )

  beta <- stats::coef(fit)
  beta[is.na(beta)] <- 0

  list(
    beta = beta,
    fitted = .clip_prob(stats::fitted(fit)),
    fit = fit
  )
}

#' Predict from a binomial calibration model
#'
#' @param L Numeric calibration design matrix.
#' @param beta Numeric calibration coefficient vector.
#'
#' @return Numeric vector of clipped calibrated probabilities.
#'
#' @keywords internal
#' @noRd
.predict_binomial_calibration <- function(L, beta) {
  .clip_prob(as.numeric(stats::plogis(L %*% beta)))
}

#' Default downstream logistic statistic
#'
#' @param X Numeric design matrix.
#' @param y Numeric binary outcome vector coded as 0 and 1.
#'
#' @return Numeric vector of logistic regression coefficients.
#'
#' @keywords internal
#' @noRd
.default_logistic_statistic <- function(X, y) {
  fit <- suppressWarnings(
    stats::glm.fit(x = X, y = y, family = stats::binomial())
  )

  beta <- stats::coef(fit)
  beta[is.na(beta)] <- NA_real_
  beta
}

#' Default downstream OLS statistic
#'
#' @param X Numeric design matrix.
#' @param y Numeric continuous outcome vector.
#'
#' @return Numeric vector of OLS coefficients.
#'
#' @keywords internal
#' @noRd
.default_ols_statistic <- function(X, y) {
  fit <- stats::lm.fit(x = X, y = y)
  beta <- stats::coef(fit)
  beta[is.na(beta)] <- NA_real_
  beta
}

#--- PPBOOT2-INSPIRED: CONTINUOUS OUTCOME / OLS TARGET -------------------------
# This is intentionally a little more conservative in naming. Efron's PPboot2
# uses a polynomial mean model m(f), a local scale model s(f), and a fitted
# standardized Gamma residual distribution. The implementation below follows
# the same mean/scale/residual-bootstrap structure, but uses empirical
# standardized residual resampling rather than fitting the Gamma residual family.

#' Fit the continuous PPboot calibration model
#'
#' @param y Numeric continuous outcome vector.
#' @param f Numeric prediction vector.
#' @param max_degree Maximum polynomial degree considered when `degree = NULL`.
#' @param degree Optional polynomial degree for the conditional mean model.
#'
#' @return A list containing the selected `degree`, fitted `mean_fit` and
#'   `scale_fit` objects, standardized residuals `eps`, and prediction
#'   functions `predict_mu()` and `predict_sd()`.
#'
#' @keywords internal
#' @noRd
.fit_ppboot2_model <- function(y, f, max_degree = 5, degree = NULL) {
  dat <- data.frame(y = as.numeric(y), f = as.numeric(f))

  if (is.null(degree)) {
    degrees <- seq_len(max_degree)

    aic <- vapply(degrees, function(d) {
      fit <- stats::lm(y ~ stats::poly(f, degree = d, raw = TRUE), data = dat)
      stats::AIC(fit)
    }, numeric(1))

    degree <- degrees[which.min(aic)]
  }

  mean_fit <- stats::lm(
    y ~ stats::poly(f, degree = degree, raw = TRUE),
    data = dat
  )

  mu <- as.numeric(stats::fitted(mean_fit))
  resid <- y - mu

  # Smooth local SD as a function of f.
  # Efron suggests a smoothed version of ordered adjacent differences
  # 0.886 * |y_(i+1) - y_(i)|. This implements that idea with loess.
  ord <- order(f)
  f_ord <- f[ord]
  y_ord <- y[ord]

  if (length(y_ord) >= 5) {
    local_sd_raw <- 0.886 * abs(diff(y_ord))
    local_f <- stats::filter(f_ord, rep(1 / 2, 2), sides = 1)[-1]
    local_dat <- data.frame(
      local_sd = pmax(as.numeric(local_sd_raw), .Machine$double.eps),
      f = as.numeric(local_f)
    )

    scale_fit <- try(
      stats::loess(log(local_sd) ~ f, data = local_dat, span = 0.75),
      silent = TRUE
    )
  } else {
    scale_fit <- NULL
  }

  pred_sd <- function(new_f) {
    if (!inherits(scale_fit, "try-error") && !is.null(scale_fit)) {
      s <- exp(as.numeric(stats::predict(scale_fit, newdata = data.frame(f = new_f))))
      s[!is.finite(s)] <- stats::sd(resid, na.rm = TRUE)
    } else {
      s <- rep(stats::sd(resid, na.rm = TRUE), length(new_f))
    }

    pmax(s, .Machine$double.eps)
  }

  sigma <- pred_sd(f)

  eps <- resid / sigma
  eps <- eps[is.finite(eps)]
  eps <- eps - mean(eps)
  eps <- eps / stats::sd(eps)

  list(
    degree = degree,
    mean_fit = mean_fit,
    scale_fit = scale_fit,
    eps = eps,
    predict_mu = function(new_f) {
      as.numeric(stats::predict(mean_fit, newdata = data.frame(f = new_f)))
    },
    predict_sd = pred_sd
  )
}
