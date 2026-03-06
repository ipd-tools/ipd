#--- PDC OLS -------------------------------------------------------------------

#' PDC OLS
#'
#' @description
#' Helper function for PDC OLS estimation.
#'
#' @details
#' Prediction de-correlated inference: A safe approach for post-prediction
#' inference (Gan et al., 2024) \doi{10.1111/anzs.12429}
#'
#' @param X_l (matrix): n x p matrix of covariates in the labeled data.
#' @param Y_l (vector): n-vector of labeled outcomes.
#' @param f_l (vector): n-vector of predictions in the labeled data.
#' @param X_u (matrix): N x p matrix of covariates in the unlabeled data.
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#' @param intercept (Logical): Do the design matrices include intercept
#' columns? Default is \code{TRUE}.
#'
#' @return (list): A list containing the following:
#' \describe{
#'    \item{est}{(vector): vector of PDC OLS regression coefficient estimates.}
#'    \item{se}{(vector): vector of standard errors of the coefficients.}
#' }
#'
#' @examples
#'
#' dat <- simdat(model = "ols")
#'
#' form <- Y - f ~ X1
#'
#' X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])
#'
#' Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>
#'   matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>
#'   matrix(ncol = 1)
#'
#' X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>
#'   matrix(ncol = 1)
#'
#' pdc_ols(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE)
#'
#' @import stats
#'
#' @export
pdc_ols <- function(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE) {
  if (!intercept) stop("ipd currently supplies design matrices with intercept; set intercept=TRUE.", call. = FALSE)

  X_a <- rbind(X_l, X_u)

  beta_lab <- .fit_glm_from_matrix(Y_l, X_l, family = "gaussian")

  pdc <- compute_pdc_est(
    Y_lab = as.numeric(Y_l),
    pred_lab = as.numeric(f_l),
    pred_all = c(as.numeric(f_l), as.numeric(f_u)),
    X_lab_int = X_l,
    X_all_int = X_a,
    beta_lab = beta_lab,
    family = "gaussian"
  )

  list(est = pdc$theta_hat, se = sqrt(diag(pdc$theta_hat_cov)))
}
