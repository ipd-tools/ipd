#===============================================================================
# METHODS
#===============================================================================

#=== STANDARD METHODS ==========================================================

#--- PRINT.IPD -----------------------------------------------------------------

#' Print IPD Fit
#'
#' Prints a brief summary of the IPD method/model combination.
#'
#' @param x An object of class \code{ipd}.
#'
#' @param ... Additional arguments to be passed to the print function.
#'
#' @return The input \code{x}, invisibly.
#'
#' @examples
#'
#' #-- Generate Example Data
#'
#' set.seed(2023)
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- Fit IPD
#'
#' fit <- ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'   data = dat, label = "set_label")
#'
#' #-- Print Output
#'
#' print(fit)
#'
#' @export

print.ipd <- function(x, ...) {

  if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

  cat("\nCall:\n", deparse(x$formula), "\n\n")

  cat("Coefficients:\n")

  print(x$coefficients)

  invisible(x)
}

#--- SUMMARY.IPD ---------------------------------------------------------------

#' Summarize IPD Fit
#'
#' Produces a summary of the IPD method/model combination.
#'
#' @param object An object of class \code{ipd}.
#'
#' @param ... Additional arguments to be passed to the summary function.
#'
#' @return A list containing:
#'
#' \describe{
#'
#'   \item{coefficients}{Model coefficients and related statistics.}
#'
#'   \item{performance}{Performance metrics of the model fit.}
#'
#'   \item{...}{Additional summary information.}
#' }
#'
#' @examples
#'
#' #-- Generate Example Data
#'
#' set.seed(2023)
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- Fit IPD
#'
#' fit <- ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'   data = dat, label = "set_label")
#'
#' #-- Summarize Output
#'
#' summ_fit <- summary(fit)
#'
#' summ_fit
#'
#' @export

summary.ipd <- function(object, ...) {

  if (!inherits(object, "ipd")) stop("Object is not of class 'ipd'")

  coef_table <- data.frame(

    Estimate = object$coefficients,

    Std.Error = object$se,

    `Lower CI` = object$ci[, 1],

    `Upper CI` = object$ci[, 2]
  )

  result <- list(

    call = object$formula,

    coefficients = coef_table,

    method = object$method,

    model = object$model,

    intercept = object$intercept
  )

  class(result) <- "summary.ipd"

  return(result)
}

#--- PRINT.SUMMARY.IPD ---------------------------------------------------------

#' Print Summary of IPD Fit
#'
#' Prints a detailed summary of the IPD method/model combination.
#'
#' @param x An object of class \code{summary.ipd}.
#'
#' @param ... Additional arguments to be passed to the print function.
#'
#' @examples
#'
#' #-- Generate Example Data
#'
#' set.seed(2023)
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- Fit IPD
#'
#' fit <- ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'   data = dat, label = "set_label")
#'
#' #-- Summarize Output
#'
#' summ_fit <- summary(fit)
#'
#' print(summ_fit)
#'
#' @return The input \code{x}, invisibly.
#'
#' @export

print.summary.ipd <- function(x, ...) {

  if (!inherits(x, "summary.ipd")) stop("Object is not of class 'summary.ipd'")

  cat("\nCall:\n", deparse(x$call), "\n\n")

  cat("Method:", x$method, "\n")

  cat("Model:", x$model, "\n")

  cat("Intercept:", ifelse(x$intercept, "Yes", "No"), "\n\n")

  cat("Coefficients:\n")

  printCoefmat(x$coefficients, P.values = FALSE, has.Pvalue = FALSE)

  invisible(x)
}

#=== BROOM TIDIER METHODS ======================================================

#--- TIDY.IPD ------------------------------------------------------------------

#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy an IPD Fit
#'
#' Tidies the IPD method/model fit into a data frame.
#'
#' @param x An object of class \code{ipd}.
#'
#' @param ... Additional arguments to be passed to the tidy function.
#'
#' @return A tidy data frame of the model's coefficients.
#'
#' @examples
#'
#' #-- Generate Example Data
#'
#' set.seed(2023)
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- Fit IPD
#'
#' fit <- ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'   data = dat, label = "set_label")
#'
#' #-- Tidy Output
#'
#' tidy(fit)
#'
#' @export

tidy.ipd <- function(x, ...) {

  if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

  result <- data.frame(

    term = names(x$coefficients),

    estimate = x$coefficients,

    std.error = x$se,

    conf.low = x$ci[, 1],

    conf.high = x$ci[, 2])

  return(result)
}

#--- GLANCE.IPD ----------------------------------------------------------------

#' @importFrom generics glance
#' @export
generics::glance

#' Glance at an IPD Fit
#'
#' Glances at the IPD method/model fit, returning a one-row summary.
#'
#' @param x An object of class \code{ipd}.
#'
#' @param ... Additional arguments to be passed to the glance function.
#'
#' @return A one-row data frame summarizing the IPD method/model fit.
#'
#' @examples
#'
#' #-- Generate Example Data
#'
#' set.seed(2023)
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- Fit IPD
#'
#' fit <- ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'   data = dat, label = "set_label")
#'
#' #-- Glance Output
#'
#' glance(fit)
#'
#' @export

glance.ipd <- function(x, ...) {

  if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

  glance_df <- data.frame(

    method = x$method,

    model = x$model,

    include_intercept = x$intercept,

    nobs_labeled = nrow(x$data_l),

    nobs_unlabeled = nrow(x$data_u),

    call = deparse(x$formula)
  )

  return(glance_df)
}

#--- AUGMENT.IPD ---------------------------------------------------------------

#' @importFrom generics augment
#' @export
generics::augment

#' Augment Data from an IPD Fit
#'
#' Augments the data used for an IPD method/model fit with additional
#' information about each observation.
#'
#' @param x An object of class \code{ipd}.
#'
#' @param data The \code{data.frame} used to fit the model. Default is
#' \code{x$data}.
#'
#' @param ... Additional arguments to be passed to the augment function.
#'
#' @return A \code{data.frame} containing the original data along with fitted
#' values and residuals.
#'
#' @examples
#'
#' #-- Generate Example Data
#'
#' set.seed(2023)
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- Fit IPD
#'
#' fit <- ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'   data = dat, label = "set_label")
#'
#' #-- Augment Data
#'
#' augmented_df <- augment(fit)
#'
#' head(augmented_df)
#'
#' @export

augment.ipd <- function(x, data = x$data_u, ...) {

  if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

  data_aug <- data

  model_matrix <- model.matrix(x$formula, data)

  fitted_values <- model_matrix %*% x$coefficients

  data_aug$.fitted <- fitted_values

  data_aug$.resid <- data_aug[, all.vars(x$formula)[1]] - fitted_values

  return(data_aug)
}

#=== END =======================================================================
