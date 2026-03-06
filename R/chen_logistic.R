#--- CHEN & CHEN LOGISTIC ------------------------------------------------------

#' Chen & Chen Logistic
#'
#' @description
#' Helper function for Chen & Chen logistic regression estimation.
#'
#' @details
#' Another look at statistical inference with machine learning-imputed data
#' (Gronsbell et al., 2026) \doi{10.48550/arXiv.2411.19908}
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
#' @param intercept (Logical): Do the design matrices include intercept
#' columns? Default is \code{TRUE}.
#'
#' @return (list): A list containing the following:
#' \describe{
#'    \item{est}{(vector): vector of Chen & Chen logistic regression coefficient
#'    estimates.}
#'    \item{se}{(vector): vector of standard errors of the coefficients.}
#' }
#'
#' @examples
#'
#' dat <- simdat(model = "logistic")
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
#' chen_logistic(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE)
#'
#' @import stats
#'
#' @export
chen_logistic <- function(X_l, Y_l, f_l, X_u, f_u, intercept = TRUE) {
  if (!intercept) stop("ipd currently supplies design matrices with intercept; set intercept=TRUE.", call. = FALSE)

  X_a <- rbind(X_l, X_u)

  beta_lab  <- .fit_glm_from_matrix(Y_l, X_l, family = "binomial")
  gamma_lab <- .fit_glm_from_matrix(f_l, X_l, family = "binomial")
  gamma_all <- .fit_glm_from_matrix(c(f_l, f_u), X_a, family = "binomial")

  cc <- compute_min_mse_est(
    Y_lab = as.numeric(Y_l),
    pred_lab = as.numeric(f_l),
    X_lab_int = X_l,
    X_all_int = X_a,
    beta_lab = beta_lab,
    gamma_lab = gamma_lab,
    gamma_all = gamma_all,
    family = "binomial"
  )

  list(est = cc$theta_hat, se = sqrt(diag(cc$theta_hat_cov)))
}
