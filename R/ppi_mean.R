#===============================================================================
# PPI MEAN ESTIMATION
#===============================================================================

#--- PPI MEAN ESTIMATION -------------------------------------------------------

#' PPI Mean Estimation
#'
#' @description
#' Helper function for PPI mean estimation
#'
#' @details
#' Prediction Powered Inference (Angelopoulos et al., 2023)
#' <https://www.science.org/doi/10.1126/science.adi6000>
#'
#' @param Y_l (vector): n-vector of labeled outcomes.
#'
#' @param f_l (vector): n-vector of predictions in the labeled data.
#'
#' @param f_u (vector): N-vector of predictions in the unlabeled data.
#'
#' @param alpha (scalar): type I error rate for hypothesis testing - values in
#' (0, 1); defaults to 0.05.
#'
#' @param alternative (string): Alternative hypothesis. Must be one of
#' \code{"two-sided"}, \code{"less"}, or \code{"greater"}.
#'
#' @return tuple: Lower and upper bounds of the prediction-powered confidence
#' interval for the mean.
#'
#' @examples
#'
#' dat <- simdat(model = "mean")
#'
#' form <- Y - f ~ 1
#'
#' Y_l <- dat[dat$set_label == "labeled",   all.vars(form)[1]] |> matrix(ncol = 1)
#'
#' f_l <- dat[dat$set_label == "labeled",   all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |> matrix(ncol = 1)
#'
#' ppi_mean(Y_l, f_l, f_u)
#'
#' @import stats
#'
#' @export

ppi_mean <- function(Y_l, f_l, f_u, alpha = 0.05, alternative = "two-sided") {

  #- Compute Dimensions of Inputs

  n <- ifelse(is.null(dim(Y_l)), length(Y_l), nrow(Y_l))

  N <- ifelse(is.null(dim(f_u)), length(f_u), nrow(f_u))

  p <- if (length(dim(f_l)) > 1) dim(f_l)[2] else 1

  est <- ppi_plusplus_mean_est(Y_l, f_l, f_u, lhat = 1)

  imputed_std <- sd(f_u) * sqrt((N - 1) / N) / sqrt(N)

  rectifier_std <- sd(Y_l - f_l) * sqrt((n - 1) / n) / sqrt(n)

  return(zconfint_generic(

    est, sqrt(imputed_std^2 + rectifier_std^2), alpha, alternative))
}

#=== END =======================================================================
