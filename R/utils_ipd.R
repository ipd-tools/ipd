#--- IPD S4 CLASS --------------------------------------------------------------

#' ipd: S4 class for inference on predicted data results
#'
#' @slot coefficients Numeric vector of parameter estimates.
#' @slot se           Numeric vector of standard errors.
#' @slot ci           Numeric matrix of confidence intervals.
#' @slot coefTable    Data frame summarizing \code{Estimate},
#'                    \code{Std. Error}, \code{z value}, and \code{Pr(>|z|)}.
#' @slot fit          The raw list returned by the helper function.
#' @slot formula      The formula used (class "formula").
#' @slot data_l       The labeled data (data.frame).
#' @slot data_u       The unlabeled data (data.frame).
#' @slot method       Character; which IPD method was used.
#' @slot model        Character; which downstream model was fitted.
#' @slot intercept    Logical; was an intercept included?
#'
#' @import methods
#'
#' @exportClass ipd

setClass("ipd",

    representation(
        coefficients = "numeric",
        se           = "numeric",
        ci           = "matrix",
        coefTable    = "data.frame",
        fit          = "list",
        formula      = "formula",
        data_l       = "data.frame",
        data_u       = "data.frame",
        method       = "character",
        model        = "character",
        intercept    = "logical"
    )
)

#--- PARSE INPUTS --------------------------------------------------------------

#' Validate and split input data
#'
#' @keywords internal
#'
#' @param data           data.frame or coercible object
#' @param label          optional column flagging labeled vs. unlabeled rows
#' @param unlabeled_data optional data.frame of unlabeled observations
#' @param na_action      how to handle missing data: "na.fail" or "na.omit"
#'
#' @return A list with components data_l and data_u (data.frames)

.parse_inputs <- function(
    data,
    label = NULL,
    unlabeled_data = NULL,
    na_action = "na.fail") {

    if (!na_action %in% c("na.fail", "na.omit")) {

        stop("na_action must be 'na.fail' or 'na.omit'")
    }

    if (!is.data.frame(data)) {

        data <- try(as.data.frame(data), silent = TRUE)

        if (inherits(data, "try-error")) {

            stop("data cannot be coerced to a data.frame")
        }
    }

    if (is.null(label) && is.null(unlabeled_data)) {

        stop("One of label or unlabeled_data must be specified")
    }

    if (!is.null(label) && !is.null(unlabeled_data)) {

        stop("Specify only one of label or unlabeled_data")
    }

    if (!is.null(label)) {

        if (!label %in% names(data)) {

            stop(sprintf("%s not found in data", label))
        }

        valid_l <- c("l", "lab", "label", "labeled", "labelled",

            "tst", "test", "true", 1, TRUE)

        valid_u <- c("u", "unlab", "unlabeled", "unlabelled",

            "val", "validation", "false", 0, FALSE)

        ids <- unique(data[[label]])

        if (!(sum(ids %in% valid_l) == 1L && sum(ids %in% valid_u) == 1L)) {

            stop(sprintf(
                "%s must flag exactly one labeled and one unlabeled set",
                label
            ))
        }

        data_l <- droplevels(data[data[[label]] %in% valid_l, ])
        data_u <- droplevels(data[data[[label]] %in% valid_u, ])

    } else {

        data_l <- data
        data_u <- unlabeled_data
    }

    list(data_l = data_l, data_u = data_u)
}

#--- DROP FACTOR LEVELS --------------------------------------------------------

#' Drop unused factor levels and report which were removed
#'
#' @keywords internal
#'
#' @param data         data.frame
#' @param factor_vars  character vector of factor column names
#'
#' @return data.frame with unused levels dropped

.drop_unused_levels <- function(
    data,
    factor_vars) {

    dropped <- lapply(

        data[, factor_vars, drop = FALSE],

        function(x) setdiff(levels(x), levels(droplevels(x)))
    )

    to_drop <- names(dropped)[lengths(dropped) > 0]

    if (length(to_drop)) {

        message("Dropped unused factor levels in these variables:")

        for (var in to_drop) {

            lvls <- dropped[[var]]

            message(" - ", var, ": ", paste(lvls, collapse = ", "))
        }

        data <- droplevels(data)
    }

    data
}

#--- WARN DIFFERING LEVELS -----------------------------------------------------

#' Warn on differing factor levels between labeled and unlabeled data
#'
#' @keywords internal
#'
#' @param data_l       labeled data.frame
#' @param data_u       unlabeled data.frame
#' @param factor_vars  character vector of factor column names
#'
#' @return Invisibly returns \code{NULL}.  Messages are printed for each
#' variable whose levels differ.

.warn_differing_levels <- function(
    data_l,
    data_u,
    factor_vars) {

    diffs <- list()

    for (var in factor_vars) {

        lv_l <- levels(data_l[[var]])

        lv_u <- levels(data_u[[var]])

        if (!identical(lv_l, lv_u)) {

            diffs[[var]] <- list(labeled = lv_l, unlabeled = lv_u)
        }
    }

    if (length(diffs)) {

        warning("Differing factor levels between labeled and unlabeled data:")

        for (var in names(diffs)) {

            lv   <- diffs[[var]]

            msg  <- sprintf(" - %s:\n    labeled = [%s]\n    unlabeled = [%s]",

                var,
                paste(lv$labeled, collapse = ", "),
                paste(lv$unlabeled, collapse = ", "))

            warning(msg, call. = FALSE)
        }
    }

    invisible(NULL)
}

#--- BUILD DESIGN --------------------------------------------------------------

#' Build design matrices and outcome vectors
#'
#' @keywords internal
#'
#' @param formula    two-sided formula `Y - f ~ X...`
#' @param data_l     labeled data.frame
#' @param data_u     unlabeled data.frame
#' @param intercept  include intercept?
#' @param na_action  "na.fail" or "na.omit"
#'
#' @return list(X_l, Y_l, f_l, X_u, f_u)

.build_design <- function(
    formula,
    data_l,
    data_u,
    intercept = TRUE,
    na_action = "na.fail") {

    vars <- all.vars(formula)

    preds <- vars[-c(1,2)]

    rhs <- if (length(preds)) paste(preds, collapse = " + ") else "1"

    mx_form <- if (intercept) {

        as.formula(paste("~", rhs))

    } else {

        as.formula(paste("~", rhs, "-1"))
    }

    mf_l <- model.frame(mx_form, data = data_l, na.action = na_action)

    X_l <- model.matrix(mx_form, mf_l)
    Y_l <- data_l[[vars[1]]] |> matrix(ncol = 1)
    f_l <- data_l[[vars[2]]] |> matrix(ncol = 1)

    mf_u <- model.frame(mx_form, data = data_u, na.action = na_action)

    X_u <- model.matrix(mx_form, mf_u)
    f_u <- data_u[[vars[2]]] |> matrix(ncol = 1)

    list(
        X_l = X_l,
        Y_l = Y_l,
        f_l = f_l,
        X_u = X_u,
        f_u = f_u
    )
}
