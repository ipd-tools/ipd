#--- IPD WRAPPER FUNCTION -----------------------------------------------------

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
#' @param method The IPD method to be used for fitting the model. Must be one
#' of \code{"chen"}, \code{"postpi_analytic"}, \code{"postpi_boot"},
#' \code{"ppi"}, \code{"ppi_a"}, \code{"ppi_plusplus"}, or \code{"pspa"}.
#' See \strong{3. Method} in the \strong{Details} below for more information.
#'
#' @param model The type of downstream inferential model to be fitted, or the
#' parameter being estimated. Must be one of \code{"mean"}, \code{"quantile"},
#' \code{"ols"}, \code{"logistic"}, or \code{"poisson"}. See \strong{4. Model}
#' in the \strong{Details} below for more information.
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
#' @param intercept \code{Logical}. Should an intercept be included in the
#' model? Default is \code{TRUE}.
#'
#' @param alpha The significance level for confidence intervals. Default is
#' \code{0.05}.
#'
#' @param alternative A string specifying the alternative hypothesis. Must be
#' one of \code{"two-sided"}, \code{"less"}, or \code{"greater"}.
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
#' second argument is provided, the function ignores the \code{label}
#' identifier and assumes the data provided are not stacked.
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
#'    \item{"chen"}{Gronsbell et al. (2025) Chen and Chen Correction}
#'    \item{"postpi_analytic"}{Wang et al. (2020) Post-Prediction Inference
#'    (PostPI) Analytic Correction}
#'    \item{"postpi_boot"}{Wang et al. (2020) Post-Prediction Inference
#'    (PostPI) Bootstrap Correction}
#'    \item{"ppi"}{Angelopoulos et al. (2023) Prediction-Powered Inference
#'    (PPI)}
#'    \item{"ppi_a"}{Gronsbell et al. (2025) PPI "All" Correction}
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
#' @return
#' An S4 object of class \code{IPD} with the following slots:
#' \describe{
#'   \item{\code{coefficients}}{Named \code{\link[base]{numeric}}
#'   vector of estimated parameters.}
#'   \item{\code{se}}{Named \code{\link[base]{numeric}}
#'   vector of standard errors.}
#'   \item{\code{ci}}{A \code{\link[base]{matrix}} of confidence intervals,
#'   with columns \code{lower} and \code{upper}.}
#'   \item{\code{coefTable}}{A \code{\link[base]{data.frame}} summarizing
#'   Estimate, Std. Error, z-value, and Pr(>|z|) (glm-style).}
#'   \item{\code{fit}}{The raw output \code{\link[base]{list}} returned by
#'   the method-specific helper function.}
#'   \item{\code{formula}}{The \code{\link[stats]{formula}} used for fitting
#'   the IPD model.}
#'   \item{\code{data_l}}{The labeled \code{\link[base]{data.frame}} used in
#'   the analysis.}
#'   \item{\code{data_u}}{The unlabeled \code{\link[base]{data.frame}} used
#'   in the analysis.}
#'   \item{\code{method}}{A \code{\link[base]{character}} string indicating
#'   which IPD method was applied.}
#'   \item{\code{model}}{A \code{\link[base]{character}} string indicating
#'   the downstream inferential model.}
#'   \item{\code{intercept}}{A \code{\link[base]{logical}} indicating whether
#'   an intercept was included.}
#' }
#'
#' @examples
#'
#' #-- Generate Example Data
#'
#' dat <- simdat(n = c(300, 300, 300), effect = 1, sigma_Y = 1)
#'
#' head(dat)
#'
#' formula <- Y - f ~ X1
#'
#' #-- Chen and Chen Correction (Gronsbell et al., 2025)
#'
#' ipd(formula,
#'   method = "chen", model = "ols",
#'   data = dat, label = "set_label"
#' )
#'
#' #-- PostPI Analytic Correction (Wang et al., 2020)
#'
#' ipd(formula,
#'   method = "postpi_analytic", model = "ols",
#'   data = dat, label = "set_label"
#' )
#'
#' #-- PostPI Bootstrap Correction (Wang et al., 2020)
#'
#' nboot <- 200
#'
#' ipd(formula,
#'   method = "postpi_boot", model = "ols",
#'   data = dat, label = "set_label", nboot = nboot
#' )
#'
#' #-- PPI (Angelopoulos et al., 2023)
#'
#' ipd(formula,
#'   method = "ppi", model = "ols",
#'   data = dat, label = "set_label"
#' )
#'
#' #-- PPI "All" (Gronsbell et al., 2025)
#'
#' ipd(formula,
#'   method = "ppi_a", model = "ols",
#'   data = dat, label = "set_label"
#' )
#'
#' #-- PPI++ (Angelopoulos et al., 2023)
#'
#' ipd(formula,
#'   method = "ppi_plusplus", model = "ols",
#'   data = dat, label = "set_label"
#' )
#'
#' #-- PSPA (Miao et al., 2023)
#'
#' ipd(formula,
#'   method = "pspa", model = "ols",
#'   data = dat, label = "set_label"
#' )
#'
#' @import stats
#' @import tibble
#'
#' @export

ipd <- function(
    formula,
    method,
    model,
    data,
    label = NULL,
    unlabeled_data = NULL,
    intercept = TRUE,
    alpha = 0.05,
    alternative = "two-sided",
    na_action = "na.fail",
    ...) {

    #- Implemented Methods and Models

    valid_methods <- c("chen", "postpi_analytic", "postpi_boot", "ppi",

        "ppi_a", "ppi_plusplus", "pspa")

    valid_models <- c("mean", "quantile", "ols", "logistic", "poisson")

    #- Identify Factor Variables in the Formula

    all_vars    <- all.vars(formula)
    preds       <- all_vars[-c(1,2)]
    factor_vars <- intersect(preds, names(Filter(is.factor, data)))

    #- Drop Unused Levels if Stacked Data

    if (!is.null(label) && is.null(unlabeled_data)) {

        data <- .drop_unused_levels(data, factor_vars)
    }

    #- Validate and Split

    inp <- .parse_inputs(data, label, unlabeled_data, na_action)

    #- Warn on Mismatched Levels

    .warn_differing_levels(inp$data_l, inp$data_u, factor_vars)

    #- Build Design Matrices

    mats <- .build_design(formula, inp$data_l, inp$data_u, intercept, na_action)

    #- Fit Model Using Method

    method <- match.arg(method, valid_methods)

    model  <- match.arg(model,  valid_models)

    helper <- get(paste(method, model, sep = "_"))

    fit <- helper(mats$X_l, mats$Y_l, mats$f_l, mats$X_u, mats$f_u, ...)

    #- Results

    est <- as.numeric(fit$est)
    se  <- as.numeric(fit$se)
    nm  <- colnames(mats$X_u)

    names(est) <- names(se) <- nm

    ci_mat <- zconfint_generic(est, se, alpha, alternative)

    rownames(ci_mat) <- nm
    colnames(ci_mat) <- c("lower", "upper")

    zval <- est / se
    pval <- 2 * pnorm(-abs(zval))

    coef_tab <- data.frame(
        Estimate     = est,
        `Std. Error` = se,
        `z value`    = zval,
        `Pr(>|z|)`   = pval,
        row.names    = nm,
        check.names  = FALSE
    )

    #- Return

    new("ipd",
        coefficients = est,
        se           = se,
        ci           = ci_mat,
        coefTable    = coef_tab,
        fit          = fit,
        formula      = formula,
        data_l       = inp$data_l,
        data_u       = inp$data_u,
        method       = method,
        model        = model,
        intercept    = intercept
    )
}
