#===============================================================================
#
#  PROGRAM: ipd.R [This is the development version of ipd.R]
#
#  AUTHORS: Awan Afiaz (aafiaz@uw.edu)
#           Kentaro Hoffman (khoffm3@uw.edu)
#           Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: A complete user-friendly wrapper function with necessary defaults
#           that calls each individual Inference-on-predicted-data (IPD) methods
#           (postpi, ppi, ppi++, popinf, etc) with a specified model
#           argument and returns a summary of outputs.
#
#  INPUTS:  1. Formula: Use one formula argument that specifies the "rectifier"
#           model formula.
#
#           NOTE: As the PostPI methods are the only ones to deviate from
#                 the "rectifier" framework, we will create the relationship
#                 model formula internally for when calling the PostPI-based
#                 methods. Please see further information in the function
#                 details found in help tab ("?ipd()").
#
#           2. Data: The data can be specified in two ways -
#
#              a) Single data argument with stacked data and a label identifier:
#                 Take one data argument where both the unlabeled and labeled
#                 data are "stacked" on top of each other and specify the
#                 string/logical/factor column that identifies the rows that are
#                 labeled and unlabeled data respectively.
#
#                 NOTE: We will consider labeled data to be defined as
#                       string  = labeled (US), labelled (UK), TRUE (character)
#                       logical = TRUE
#                       factor  = non-reference category (i.e. binary 1)
#
#              b) Two data arguments, one for labeled and one for unlabeled.
#
#                 NOTE: If the 2nd data argument is provided the function does
#                       not check for (or ignores the label) identifier.
#                       Otherwise, if the 2nd data argument is unspecified, the
#                       function assumes the data provided is stacked.
#
#           3. Method: Use a "method" argument to signify which specific
#              method to use. Note that, going by the author-specified
#              convention, we denote the methods as following:
#
#              a) postpi       = Wang et al. post prediction inference method
#              b) ppi          = Angelopoulos et al. prediction-powered inference
#                                method (ppi)
#              c) ppi_plusplus = Angelopoulos et al. PPI++ method
#              d) popinf       = Miao et al. Assumption-lean and Data-adaptive
#                                Post-Prediction Inference
#
#           4. Model: Use a model "model" argument to signify what model
#              (e.g., ols, logistic, cox, etc). The wrapper will concatenate
#              the Method and Model arguments to identify the required helper
#              function, all of which will have the same naming convention:
#              "method_model.R".
#
#           5. Auxiliary arguments: The wrapper will take method-specific
#              auxiliary arguments and pass them to the helper function through
#              the triple dots "..." with specified defaults for simplicity.
#
#           6. Other arguments: All other arguments that relate to all the
#              methods (e.g., alpha, ci.type),
#              or all other method-specific arguments and will have defaults.
#
#           NOTES:
#              i) Propose all of the other methods take matrix arguments, to
#                 minimize additional parsing after wrapper function assertions
#                 and data pre-processing.
#
#  Updated: 2024-01-21
#
#===============================================================================

#=== MAIN IPD FUNCTION =========================================================

#' Valid and Efficient Inference on Predicted Data (IPD) using State-of-the-Art
#' Methods
#'
#' @description A complete wrapper function with necessary defaults that calls
#' on each individual post-prediction inference based methods (postpi, ppi,
#' ppi++, popinf, etc) along with a specified model argument and
#' returns a summary of outputs.
#'
#' @details  Note that, going by the author-specified conventions, we denote
#' the methods as following:
#' (a) postpi       = Wang et al. (2020) post prediction inference (PostPI) method
#' (b) ppi          = Angelopoulos et al. (2023) prediction-powered inference
#'                    (PPI) method
#' (c) ppi_plusplus = Angelopoulos et al. (2023) PPI++ method
#' (d) popinf       = Miao et al. (2023) Assumption-lean and Data-adaptive
#'                    Post-Prediction Inference method
#'
#' Since the PostPI methods (requiring a relationship model) are
#' the only ones to deviate from the 'rectifier' model framework (used in PPI,
#' PPI++, and POP-Inf), we will create the relationship model formula internally
#' when calling the PostPI methods under the assumption that the relationship
#' model arguments can be recovered from the rectifier model.
#'
#' @param formula (formula): an argument of the form Y - Yhat ~ X, where Y is the
#' name of the column corresponding to the observed outcome in the labeled data,
#' Yhat is the name of the column corresponding to the predicted outcome in both
#' labeled and unlabeled data, and X corresponds to the features of interest
#' (i.e., X = X1 + X2 + ... + Xp).
#'
#' @param method (string): Method to be used. The available methods include
#' 'postpi' (PostPI), 'ppi' (PPI), 'ppi_plusplus' (PPI++), and 'popinf' (POP-Inf).
#'
#' @param model (string): Type of the regression model to be fit or
#' estimand to be calculated.
#' This can take on values: 'ols' (for continuous  variable), 'logistic'
#' (for binary variable), and estimands: 'mean' and 'quantile'.
#' Note that, not all method-model combinations are currently available. For newer
#' method-model additions, please use the development version from the github
#' repo: awanafiaz/IPD.
#'
#' @param data (data.frame): either a stacked data frame with a specific column identifying the
#' labeled and unlabeled data sets, or only the labeled data set.
#' The stacked data frame (or the labeled data set) consists of observed outcomes
#' for labeled data set (Y), the predicted outcomes (Yhat) for both labeled and
#' unlabeled data sets, and the features (X).
#'
#' @param label_index (string/int/logical, optional): The column name for the
#' indexing variable that identifies the labeled and unlabeled data sets. This
#' column must have the following pairs of values for labeled and unlabeled data:
#' string: 'labeled' and 'unlabeled', or integer: 1 (labeled) and 0 (unlabeled), or
#' logical: TRUE (labeled) and FALSE (unlabled)).
#'
#' @param df_unlabeled (data.frame, optional): a data frame containing only the
#' unlabeled data set. Specify this argument ONLY if the data provided is not
#' already stacked. This data frame consists of the predicted outcomes (Yhat)
#' for the unlabeled data and the features (X). Specifying both 'label_index'
#' and 'df_unlabeled' arguments will result in an error message.
#'
#' @param seed (int, optional): an integer value provided as seed.
#'
#' @param alpha (float, optional): specified level of significance of the the
#' confidence interval; defaults to 0.05.
#'
#' @param alternative (string, optional): Specify the alternative hypothesis,
#' must be one of "two.sided", "greater" or "less"; defaults to "two.sided".
#'
#' @param ... further arguments passed to or from specific methods. See details.
#'
#' @returns a summary of model output.
#'
#' @examples
#' # example code
#'
#'
#' @import stats
#'
#' @export

ipd <- function(formula, method, model, data,

  label_index = NULL, df_unlabeled = NULL, seed = NULL,

  alpha = 0.05, alternative = "two-sided", ...) {

  #--- 1. CHECKS & ASSERTIONS --------------------------------------------------

  #-- A1. CHECK FOR DATA

  if (missing(data)) data <- environment(formula)

  #-- A2. CHECK IF BOTH 'label_index' AND 'df_unlabeled' args are UNSPECIFIED

  if(is.null(label_index) & is.null(df_unlabeled)) {

    stop(paste("at least one of 'label_index' and 'df_unlabeled' must be specified.",

      "\nSee the help('ipd') documentation for more information."))
  }

  #-- A3. CHECK IF BOTH 'label_index' AND 'df_unlabeled' args are SPECIFIED

  if(!is.null(label_index) & !is.null(df_unlabeled)) {

    stop(paste("specify only one of 'label_index' and 'df_unlabeled' argument.",

      "\nSee the help('ipd') documentation for more information."))
  }

  #-- A4. CHECK IF SPECIFIED 'label_index' exists in data

  if (!is.null(label_index)) {

    if (!exists(label_index, where = data)) {

      stop(paste(label_index, "does not exist in the data set.",

        "\nSee the help('ipd') documentation for more information."))
    }
  }

  #-- B. CHECK FOR VALID METHOD

  if (!(method %in% c("postpi_analytic", "postpi_boot", "postpi_mi",

    "ppi", "popinf", "ppi_plusplus"))) {

    stop(paste("'method' must be one of",

      "c('postpi_analytic', 'postpi_boot', 'postpi_mi', 'ppi', 'popinf').",

      "See the 'Details' section of the documentation for more information."))
  }

  #-- C. CHECK FOR VALID MODEL

  if (!(model %in% c("mean", "quantile", "ols", "logistic", "multiclass"))) {

    stop(paste("'model' must be one of",

      "c('mean', 'quantile', 'ols', 'logistic', 'multiclass').",

      "See the 'Details' section of the documentation for more information."))
  }

  #--- D. CHECK FOR VALID METHOD-MODEL COMBINATION

  # if (!exists(paste(method, model, sep = "_"), inherits = FALSE))  {
  #
  #   stop(paste0("This method/model combination (",
  #
  #     paste(method, model, sep = "_"), ") has not yet been implemented. ",
  #
  #     "See the 'Details' section of the documentation for more information."))
  # }

  #--- E. SET SEED

  if (!is.null(seed)) set.seed(seed)

  #--- PREPARE DATA ------------------------------------------------------------

  #---- DEFINE VALID label_index IDENTIFERS

  ##--- CHECK ONE OF EACH labeled AND unlabeled IDENTIFIERS EXIST

  valid_labeled_df_id <- c(

    "tst", "test", 1, TRUE, "lab", "label", "labeled", "labelled")

  valid_unlabeled_df_id <- c(

    "val", "validation", 0, FALSE, "unlab", "unlabeled", "unlabelled")

  if(!((sum(unique(data[[label_index]]) %in% valid_labeled_df_id) == 1) &

    (sum(unique(data[[label_index]]) %in% valid_unlabeled_df_id) == 1))) {

    stop(paste(label_index,

      "must have one valid identifier for labeled and unlabeled data set each.",

      "See the 'Details' section of the documentation for more information."))
  }

  ##-- IF STACKED DATA IS PROVIDED

  if(!is.null(label_index) & is.null(df_unlabeled)) {

    data_l <- data[data[[label_index]] %in% valid_labeled_df_id, ]

    data_u <- data[data[[label_index]] %in% valid_unlabeled_df_id, ]
  }

  if(is.null(label_index) & !is.null(df_unlabeled)) {

    ##-- IF UNSTACKED DATA IS PROVIDED

    data_l <- data

    data_u <- df_unlabeled
  }

  #-- LABELED DATA

  X_l <- model.matrix(formula, data = data_l)

  Y_l <- data_l[ , all.vars(formula)[1]] |> matrix(ncol = 1)

  f_l <- data_l[ , all.vars(formula)[2]] |> matrix(ncol = 1)

  #-- UNLABELED DATA

  X_u <- model.matrix(formula, data = data_u)

  f_u <- data_u[ , all.vars(formula)[2]] |> matrix(ncol = 1)

  #--- METHOD ------------------------------------------------------------------

  func <- get(paste(method, model, sep = "_"))

  fit <- func(X_l, Y_l, f_l, X_u, f_u)

  ci  <- zconfint_generic(fit$est, fit$se, alpha, alternative)

  #--- RETURN ------------------------------------------------------------------

  obj <- list(est = fit$est, se = fit$se, ci = ci, ...)

  class(obj) <- c("ipd", fit$class)

  return(obj)
}

#=== END =======================================================================
