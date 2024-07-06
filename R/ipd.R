#===============================================================================
#
#  FILE:    ipd.R
#
#  PURPOSE: Defines main IPD wrapper function
#
#  UPDATED: 2024.07.05
#
#  NOTES:   1. Add logistic, Poisson, and multiclass models to all methods
#
#===============================================================================

#=== IPD WRAPPER FUNCTION ======================================================

#' Valid and Efficient Inference on Predicted Data (IPD)
#'
#' @description The main wrapper function to conduct IPD using various methods
#' and models, and returns a list of fitted model components.
#'
#' @param formula An object of class \code{formula}: a symbolic description of
#' the model to be fitted. Must be of the form \code{Y - f ~ X}, where \code{Y}
#' is the name of the column corresponding to the observed outcome in the
#' labeled data, \code{f} is the name of the column corresponding to the
#' predicted outcome in both labeled and unlabeled data, and \code{X}
#' corresponds to the features of interest (i.e., \code{X = X1 + ... + Xp}).
#'
#' @param method The method to be used for fitting the model. Must be one of
#' \code{"postpi_analytic"}, \code{"postpi_boot"}, \code{"ppi"},
#' \code{"popinf"}, or \code{"ppi_plusplus"}.
#'
#' @param model The type of model to be fitted. Must be one of \code{"mean"},
#' \code{"quantile"}, \code{"ols"}, or \code{"logistic"}.
#'
#' @param data A \code{data.frame} containing the variables in the model,
#' either a stacked data frame with a specific column identifying the labeled
#' versus unlabeled observations (\code{label}), or only the labeled data
#' set. Must contain columns for the observed outcomes (\code{Y}), the
#' predicted outcomes (\code{f}), and the features (\code{X}) needed to specify
#' the \code{formula}.
#'
#' @param label A \code{string}, \code{int}, or \code{logical} specifying the
#' column in the data that distinguishes between the labeled and unlabeled
#' observations. See the \code{Details} section for more information. If NULL,
#' \code{unlabeled_data} must be specified.
#'
#' @param unlabeled_data (optional) A \code{data.frame} of unlabeled data. If
#' NULL, \code{label} must be specified. Specifying both the \code{label} and
#' \code{unlabeled_data} arguments will result in an error message. If
#' specified, must contain columns for the predicted outcomes (\code{f}), and
#' the features (\code{X}) needed to specify the \code{formula}.
#'
#' @param seed (optional) An \code{integer} seed for random number generation.
#'
#' @param intercept \code{Logical}. Should an intercept be included in the
#' model? Default is \code{TRUE}.
#'
#' @param alpha The significance level for confidence intervals. Default is
#' \code{0.05}.
#'
#' @param alternative A string specifying the alternative hypothesis. Must be
#' one of \code{"two-sided"}, \code{"less"}, or \code{"greater"}.
#'
#' @param ... Additional arguments to be passed to the fitting function. See
#' the \code{Details} section for more information.
#'
#' @returns a summary of model output.
#'
#' @details
#'
#' \strong{1. Formula:}
#'
#' The \code{ipd} function uses one formula argument that specifies both the
#' calibrating model (e.g., PostPI "relationship model", PPI "rectifier" model)
#' and the inferential model. These separate models will be created internally
#' based on the specific \code{method} called.
#'
#' \strong{2. Data:}
#'
#' The data can be specified in two ways:
#'
#' \enumerate{
#'    \item Single data argument (\code{data}) containing a stacked
#'    \code{data.frame} and a label identifier (\code{label}).
#'    \item Two data arguments, one for the labeled data (\code{data}) and one
#'    for the unlabeled data (\code{unlabeled_data}).
#' }
#'
#' For option (1), provide one data argument (\code{data}) which contains a
#' stacked \code{data.frame} with both the unlabeled and labeled data and a
#' \code{label} argument that specify the column that identifies the labeled
#' versus the unlabeled observations in the stacked \code{data.frame}
#'
#' NOTE: Labeled data identifiers can be:
#'
#' \describe{
#'    \item{String}{"l", "lab", "label", "labeled", "labelled", "tst", "test",
#'    "true"}
#'    \item{Logical}{TRUE}
#'    \item{Factor}{Non-reference category (i.e., binary 1)}
#' }
#'
#' Unlabeled data identifiers can be:
#'
#' \describe{
#'    \item{String}{"u", "unlab", "unlabeled", "unlabelled", "val",
#'    "validation", "false"}
#'    \item{Logical}{FALSE}
#'    \item{Factor}{Non-reference category (i.e., binary 0)}
#' }
#'
#' For option (2), provide separate data arguments for the labeled data set
#' (\code{data}) and the unlabeled data set (\code{unlabeled_data}). If the
#' second argument is provided, the function ignores the label identifier and
#' assumes the data provided is stacked.
#'
#' \strong{3. Method:}
#'
#' Use the \code{method} argument to specify the fitting method:
#'
#' \describe{
#'    \item{"postpi"}{Wang et al. (2020) Post-Prediction Inference (PostPI)}
#'    \item{"ppi"}{Angelopoulos et al. (2023) Prediction-Powered Inference
#'    (PPI)}
#'    \item{"ppi_plusplus"}{Angelopoulos et al. (2023) PPI++}
#'    \item{"popinf"}{Miao et al. (2023) Assumption-Lean and Data-Adaptive
#'    Post-Prediction Inference (POP-Inf)}
#' }
#'
#' \strong{4. Model:}
#'
#' Use the \code{model} argument to specify the type of model:
#'
#' \describe{
#'    \item{"mean"}{Mean value of the outcome}
#'    \item{"quantile"}{\code{q}th quantile of the outcome}
#'    \item{"ols"}{Linear regression}
#'    \item{"logistic"}{Logistic regression}
#'    \item{"poisson"}{Poisson regression}
#' }
#'
#' The \code{ipd} wrapper function will concatenate the \code{method} and
#' \code{model} arguments to identify the required helper function, following
#' the naming convention "method_model".
#'
#' \strong{5. Auxiliary Arguments:}
#'
#' The wrapper function will take method-specific auxiliary arguments (e.g.,
#' \code{q} for the quantile estimation models) and pass them to the helper
#' function through the "..." with specified defaults for simplicity.
#'
#' \strong{6. Other Arguments:}
#'
#' All other arguments that relate to all methods (e.g., alpha, ci.type), or
#' other method-specific arguments, will have defaults.
#'
#' @return A list containing the fitted model components:
#'
#' \describe{
#'   \item{coefficients}{Estimated coefficients of the model}
#'   \item{se}{Standard errors of the estimated coefficients}
#'   \item{ci}{Confidence intervals for the estimated coefficients}
#'   \item{formula}{The formula used to fit the IPD model.}
#'   \item{data}{The data frame used for model fitting.}
#'   \item{method}{The method used for model fitting.}
#'   \item{model}{The type of model fitted.}
#'   \item{intercept}{Logical. Indicates if an intercept was included in the
#'   model.}
#'   \item{fit}{Fitted model object containing estimated coefficients, standard
#'   errors, confidence intervals, and additional method-specific output.}
#'   \item{...}{Additional output specific to the method used.}
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
#' #-- PostPI Analytic Correction (Wang et al., 2020)
#'
#' ipd(formula, method = "postpi_analytic", model = "ols",
#'
#'     data = dat, label = "set")
#'
#' #-- PostPI Bootstrap Correction (Wang et al., 2020)
#'
#' nboot <- 200
#'
#' ipd(formula, method = "postpi_boot", model = "ols",
#'
#'     data = dat, label = "set", nboot = nboot)
#'
#' #-- PPI (Angelopoulos et al., 2023)
#'
#' ipd(formula, method = "ppi", model = "ols",
#'
#'     data = dat, label = "set")
#'
#' #-- PPI++ (Angelopoulos et al., 2023)
#'
#' ipd(formula, method = "ppi_plusplus", model = "ols",
#'
#'     data = dat, label = "set")
#'
#' #-- POP-Inf (Miao et al., 2023)
#'
#' ipd(formula, method = "popinf", model = "ols",
#'
#'     data = dat, label = "set")
#'
#' @import stats
#'
#' @export

ipd <- function(formula, method, model, data,

  label = NULL, unlabeled_data = NULL, seed = NULL,

  intercept = T, alpha = 0.05, alternative = "two-sided", ...) {

  #--- CHECKS & ASSERTIONS -----------------------------------------------------

  #-- CHECK FOR DATA

  if (missing(data)) data <- environment(formula)

  #-- CHECK IF BOTH 'label' AND 'unlabeled_data' ARE UNSPECIFIED

  if(is.null(label) & is.null(unlabeled_data)) {

    stop(paste("at least one of 'label' and 'unlabeled_data' must be",

      "specified.\nSee the help('ipd') documentation for more information."))
  }

  #-- CHECK IF BOTH 'label' AND 'unlabeled_data' ARE SPECIFIED

  if(!is.null(label) & !is.null(unlabeled_data)) {

    stop(paste("specify only one of 'label' and 'unlabeled_data' argument.",

      "\nSee the help('ipd') documentation for more information."))
  }

  #-- CHECK IF SPECIFIED 'label' EXISTS IN DATA

  if (!is.null(label)) {

    if (!exists(label, where = data)) {

      stop(paste(label, "does not exist in the data set.\nSee the",

        "help('ipd') documentation for more information."))
    }
  }

  #-- CHECK FOR VALID METHOD

  if (!(method %in% c("postpi_analytic", "postpi_boot", "ppi", "popinf",

    "ppi_plusplus"))) {

    stop(paste("'method' must be one of c('postpi_analytic', 'postpi_boot',",

      "'ppi', 'popinf', 'ppi_plusplus').\nSee the 'Details' section of the",

      "documentation for more information."))
  }

  #-- CHECK FOR VALID MODEL

  if (!(model %in% c("mean", "quantile", "ols", "logistic", "poisson"))) {

    stop(paste("'model' must be one of c('mean', 'quantile', 'ols',",

      "'logistic', 'poisson').\nSee the 'Details' for more information."))
  }

  #--- SET SEED ----------------------------------------------------------------

  if (!is.null(seed)) set.seed(seed)

  #--- PREPARE DATA ------------------------------------------------------------

  #-- IF STACKED DATA ARE PROVIDED

  if(!is.null(label) & is.null(unlabeled_data)) {

    #---- DEFINE VALID label IDENTIFERS

    ##--- CHECK ONE OF EACH labeled AND unlabeled IDENTIFIERS EXIST

    valid_labeled_df_id <- c("l", "lab", "label", "labeled", "labelled",

      "tst", "test", "true", 1, TRUE)

    valid_unlabeled_df_id <- c("u", "unlab", "unlabeled", "unlabelled",

      "val", "validation", "false", 0, FALSE)

    if(!((sum(unique(data[[label]]) %in% valid_labeled_df_id) == 1) &

      (sum(unique(data[[label]]) %in% valid_unlabeled_df_id) == 1))) {

      stop(paste(label,

        "must have one valid identifier for labeled and unlabeled data set",

        "each. See the 'Details' section of the documentation for more",

        "information."))
    }

    data_l <- data[data[[label]] %in% valid_labeled_df_id, ]

    data_u <- data[data[[label]] %in% valid_unlabeled_df_id, ]
  }

  if(is.null(label) & !is.null(unlabeled_data)) {

    #- IF UNSTACKED DATA ARE PROVIDED

    data_l <- data

    data_u <- unlabeled_data
  }

  #-- LABELED DATA

  if (intercept) {

    X_l <- model.matrix(formula, data = data_l)

  } else {

    X_l <- model.matrix(formula - 1, data = data_l)
  }

  Y_l <- data_l[ , all.vars(formula)[1]] |> matrix(ncol = 1)

  f_l <- data_l[ , all.vars(formula)[2]] |> matrix(ncol = 1)

  #-- UNLABELED DATA

  if (intercept) {

    X_u <- model.matrix(formula, data = data_u)

  } else {

    X_u <- model.matrix(formula - 1, data = data_u)
  }

  f_u <- data_u[ , all.vars(formula)[2]] |> matrix(ncol = 1)

  #--- METHOD ------------------------------------------------------------------

  func <- get(paste(method, model, sep = "_"))

  fit <- func(X_l, Y_l, f_l, X_u, f_u)

  cat(colnames(X_u))

  names(fit$est) <- colnames(X_u)

  ci  <- zconfint_generic(fit$est, fit$se, alpha, alternative)

  #--- RETURN ------------------------------------------------------------------

  result <- list(

    coefficients = fit$est,
    se = fit$se,
    ci = ci,
    formula = formula,
    data = data,
    method = method,
    model = model,
    intercept = intercept,
    fit = fit,
    ...
  )

  class(result) <- c("ipd")

  return(result)
}

#=== END =======================================================================
