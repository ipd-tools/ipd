#===============================================================================
# WRAPPER FUNCTION
#===============================================================================

#--- MAIN WRAPPER FUNCTION -----------------------------------------------------

#' Inference on Predicted Data (ipd)
#'
#' @description The main wrapper function to conduct ipd using various methods
#' and models, and returns a list of fitted model components.
#'
#' @param formula An object of class \code{formula}: a symbolic description of
#' the model to be fitted. Must be of the form \code{Y - f ~ X}, where \code{Y}
#' is the name of the column corresponding to the observed outcome in the
#' labeled data, \code{f} is the name of the column corresponding to the
#' predicted outcome in both labeled and unlabeled data, and \code{X}
#' corresponds to the features of interest (i.e., \code{X = X1 + ... + Xp}).
#' See \strong{1. Formula} in the \strong{Details} below for more information.
#'
#' @param method The IPD method to be used for fitting the model. Must be one of
#' \code{"postpi_analytic"}, \code{"postpi_boot"}, \code{"ppi"},
#' \code{"ppi_plusplus"}, or \code{"pspa"}.
#' See \strong{3. Method} in the \strong{Details} below for more information.
#'
#' @param model The type of downstream inferential model to be fitted, or the
#' parameter being estimated. Must be one of \code{"mean"},
#' \code{"quantile"}, \code{"ols"}, \code{"logistic"}, or \code{"poisson"}.
#' See \strong{4. Model} in the \strong{Details} below for more information.
#'
#' @param data A \code{data.frame} containing the variables in the model,
#' either a stacked data frame with a specific column identifying the labeled
#' versus unlabeled observations (\code{label}), or only the labeled data
#' set. Must contain columns for the observed outcomes (\code{Y}), the
#' predicted outcomes (\code{f}), and the features (\code{X}) needed to specify
#' the \code{formula}. See \strong{2. Data} in the \strong{Details} below for
#' more information.
#'
#' @param label A \code{string}, \code{int}, or \code{logical} specifying the
#' column in the data that distinguishes between the labeled and unlabeled
#' observations. See the \code{Details} section for more information. If NULL,
#' \code{unlabeled_data} must be specified. See \strong{2. Data} in the
#' \strong{Details} below for more information.
#'
#' @param unlabeled_data (optional) A \code{data.frame} of unlabeled data. If
#' NULL, \code{label} must be specified. Specifying both the \code{label} and
#' \code{unlabeled_data} arguments will result in an error message. If
#' specified, must contain columns for the predicted outcomes (\code{f}), and
#' the features (\code{X}) needed to specify the \code{formula}. See
#' \strong{2. Data} in the \strong{Details} below for more information.
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
#' @param n_t (integer, optional) Size of the dataset used to train the
#' prediction function (necessary for the \code{"postpi_analytic"} and
#' \code{"postpi_boot"} methods if \code{n_t} < \code{nrow(X_l)}.
#' Defaults to \code{Inf}.
#'
#' @param na_action (string, optional) How missing covariate data should be
#' handled. Currently \code{"na.fail"} and \code{"na.omit"} are accommodated.
#' Defaults to \code{"na.fail"}.
#'
#' @param ... Additional arguments to be passed to the fitting function. See
#' the \code{Details} section for more information. See
#' \strong{5. Auxiliary Arguments} and \strong{6. Other Arguments} in the
#' \strong{Details} below for more information.
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
#' \code{label} argument that specifies the column identifying the labeled
#' versus the unlabeled observations in the stacked \code{data.frame} (e.g.,
#' \code{label = "set_label"} if the column "set_label" in the stacked data
#' denotes which set an observation belongs to).
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
#' second argument is provided, the function ignores the \code{label} identifier
#' and assumes the data provided are not stacked.
#'
#' NOTE: Not all columns in \code{data} or \code{unlabeled_data} may be used
#' unless explicitly referenced in the \code{formula} argument or in the
#' \code{label} argument (if the data are passed as one stacked data frame).
#'
#' \strong{3. Method:}
#'
#' Use the \code{method} argument to specify the fitting method:
#'
#' \describe{
#'    \item{"postpi_analytic"}{Wang et al. (2020) Post-Prediction Inference (PostPI) Analytic Correction}
#'    \item{"postpi_boot"}{Wang et al. (2020) Post-Prediction Inference (PostPI) Bootstrap Correction}
#'    \item{"ppi"}{Angelopoulos et al. (2023) Prediction-Powered Inference
#'    (PPI)}
#'    \item{"ppi_plusplus"}{Angelopoulos et al. (2023) PPI++}
#'    \item{"pspa"}{Miao et al. (2023) Assumption-Lean and Data-Adaptive
#'    Post-Prediction Inference (PSPA)}
#' }
#'
#' \strong{4. Model:}
#'
#' Use the \code{model} argument to specify the type of downstream inferential
#' model or parameter to be estimated:
#'
#' \describe{
#'    \item{"mean"}{Mean value of a continuous outcome}
#'    \item{"quantile"}{\code{q}th quantile of a continuous outcome}
#'    \item{"ols"}{Linear regression coefficients for a continuous outcome}
#'    \item{"logistic"}{Logistic regression coefficients for a binary outcome}
#'    \item{"poisson"}{Poisson regression coefficients for a count outcome}
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
#'   \item{formula}{The formula used to fit the ipd model.}
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
#' set.seed(12345)
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
#'     data = dat, label = "set_label")
#'
#' #-- PostPI Bootstrap Correction (Wang et al., 2020)
#'
#' nboot <- 200
#'
#' ipd(formula, method = "postpi_boot", model = "ols",
#'
#'     data = dat, label = "set_label", nboot = nboot)
#'
#' #-- PPI (Angelopoulos et al., 2023)
#'
#' ipd(formula, method = "ppi", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' #-- PPI++ (Angelopoulos et al., 2023)
#'
#' ipd(formula, method = "ppi_plusplus", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' #-- PSPA (Miao et al., 2023)
#'
#' ipd(formula, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' @import stats
#'
#' @export

ipd <- function(formula, method, model, data,

  label = NULL, unlabeled_data = NULL, seed = NULL, intercept = TRUE,

  alpha = 0.05, alternative = "two-sided", n_t = Inf, na_action = "na.fail",

  ...) {

  #--- CHECKS & ASSERTIONS -----------------------------------------------------

  #-- CHECK ARGUMENTS

  if(na_action != "na.fail" & na_action != "na.omit") {

    stop("na_action should be either 'na.fail' or 'na.omit'")
  }

  #-- CHECK FOR DATA

  if (missing(data)) data <- environment(formula)

  #-- CHECK DATA CLASS

  if (!inherits(data, "data.frame")) {

    data <- try(as.data.frame(data), silent = TRUE)

    if (inherits(data, "try-error")) {

      stop(paste("'data' cannot be coerced to a data.frame.",

        "Please provide a valid data.frame."))
    }
  }

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

  if (!(method %in% c("postpi_analytic", "postpi_boot", "ppi", "pspa",

    "ppi_plusplus"))) {

    stop(paste("'method' must be one of c('postpi_analytic', 'postpi_boot',",

      "'ppi', 'pspa', 'ppi_plusplus').\nSee the 'Details' section of the",

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

  #-- DROP UNUSED FACTOR LEVELS IN DATA AND REPORT DROPPED LEVELS

  factor_vars <- names(Filter(is.factor, data))

  #-- IF STACKED DATA ARE PROVIDED

  if (!is.null(label) & is.null(unlabeled_data)) {

    if (!is.null(factor_vars)) {

      dropped_levels <- sapply(data[factor_vars],

        function(x) setdiff(levels(x), levels(droplevels(x))))

      if (any(lengths(dropped_levels) > 0)) {

        message("Dropped unused factor levels in the following variables:\n",

                paste(names(dropped_levels)[lengths(dropped_levels) > 0],

                      ":", sapply(dropped_levels[lengths(dropped_levels) > 0],

                                  paste, collapse = ", "), collapse = "\n"))
      }

      data <- droplevels(data)
    }

    #- DEFINE VALID 'label' IDENTIFERS

    # CHECK ONE OF EACH 'labeled' AND 'unlabeled' IDENTIFIERS EXIST

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

  #-- IF UNSTACKED DATA ARE PROVIDED

  if(is.null(label) & !is.null(unlabeled_data)) {

    data_l <- data

    data_u <- unlabeled_data
  }

  #-- CHECK FOR UNUSED FACTOR LEVELS AFTER SPLITTING

  data_l <- droplevels(data_l)

  data_u <- droplevels(data_u)

  differing_levels <- sapply(factor_vars, function(var) {

    levels_l <- levels(data_l[[var]])

    levels_u <- levels(data_u[[var]])

    if (!identical(levels_l, levels_u)) {

      return(paste("Differing levels in '", var, "': labeled = [",

        paste(levels_l, collapse = ", "), "], unlabeled = [",

        paste(levels_u, collapse = ", "), "]", sep = ""))
    }

    return(NULL)
  })

  differing_levels <- differing_levels[!sapply(differing_levels, is.null)]

  if (length(differing_levels) > 0) {

    stop(paste0("The following variables have differing factor levels ",

      "between the labeled and unlabeled data:\n",

      paste(differing_levels, collapse = "\n")))
  }

  #-- LABELED DATA

  #- CHECK FOR MISSING COVARIATE DATA

  missing_covs_l <- apply(data_l[ , all.vars(formula)[-(1:2)], drop = FALSE], 2,

    function(x) any(is.na(x)))

  if (any(missing_covs_l)) {

    missing_vars_l <- names(missing_covs_l[missing_covs_l])

    stop(paste0(

      "Missing values detected in the following labeled covariate data: ",

      paste(missing_vars_l, collapse = ", "),

      ". Please ensure there are no missing values in these covariates."
    ))
  }

  if (intercept) {

    X_l <- model.matrix(formula,

      model.frame(formula, data = data_l, na.action = na_action))

  } else {

    X_l <- model.matrix(update(formula, . ~ . - 1),

      model.frame(update(formula, . ~ . - 1),

        data = data_l, na.action = na_action))
  }

  Y_l <- data_l[ , all.vars(formula)[1]] |> matrix(ncol = 1)

  f_l <- data_l[ , all.vars(formula)[2]] |> matrix(ncol = 1)

  #-- UNLABELED DATA

  formula_u <- as.formula(paste(all.vars(formula)[2], "~",

    paste(all.vars(formula)[-c(1, 2)], collapse = " + ")))


  #- CHECK FOR MISSING COVARIATE DATA

  missing_covs_u <- apply(data_u[ , all.vars(formula_u)[-1], drop = FALSE], 2,

    function(x) any(is.na(x)))

  if (any(missing_covs_u)) {

    missing_vars_u <- names(missing_covs_u[missing_covs_u])

    stop(paste0(

      "Missing values detected in the following unlabeled covariate data: ",

      paste(missing_vars_u, collapse = ", "),

      ". Please ensure there are no missing values in these covariates."
    ))
  }

  if (intercept) {

    X_u <- model.matrix(formula_u,

      model.frame(formula_u, data = data_u, na.action = na_action))

  } else {

    X_u <- model.matrix(update(formula_u, . ~ . - 1),

      model.frame(update(formula_u, . ~ . - 1),

        data = data_u, na.action = na_action))
  }

  f_u <- data_u[ , all.vars(formula)[2]] |> matrix(ncol = 1)

  #--- METHOD ------------------------------------------------------------------

  func <- get(paste(method, model, sep = "_"))

  if(grepl("postpi", method) && model == "ols") {

    fit <- func(X_l, Y_l, f_l, X_u, f_u, n_t = n_t, ...)

  } else {

    fit <- func(X_l, Y_l, f_l, X_u, f_u, ...)
  }

  names(fit$est) <- colnames(X_u)

  ci  <- zconfint_generic(fit$est, fit$se, alpha, alternative)

  #--- RETURN ------------------------------------------------------------------

  result <- list(

    coefficients = fit$est,
    se = fit$se,
    ci = ci,
    formula = formula,
    data_l = data_l,
    data_u = data_u,
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
