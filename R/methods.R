#--- S4 SHOW METHOD ------------------------------------------------------------

#' Show an ipd object
#'
#' @description
#' Display a concise summary of an \code{ipd} S4 object, including
#' method, model, formula, and a glm-style coefficient table.
#'
#' @param object An object of S4 class \code{ipd}.
#'
#' @return Invisibly returns \code{object} after printing.
#'
#' @import methods
#'
#' @exportMethod show

setMethod("show", "ipd", function(object) {

    cat("IPD inference summary\n")
    cat("  Method:  ",  object@method,  "\n")
    cat("  Model:   ",  object@model,   "\n")
    cat("  Formula: ",  deparse(object@formula), "\n\n")
    cat("Coefficients:\n")

    printCoefmat(object@coefTable, P.values = TRUE, has.Pvalue = TRUE)

    invisible(object)
})

#--- S3 PRINT.IPD --------------------------------------------------------------

#' Print ipd fit
#'
#' @param x   An object of class \code{ipd}.
#'
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export

print.ipd <- function(x, ...) {

    if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

    show(x)

    invisible(x)
}

#--- S3 SUMMARY.IPD ------------------------------------------------------------

#' Summarize ipd fit
#'
#' @param object An object of class \code{ipd}.
#'
#' @param ...    Ignored.
#'
#' @return An object of class \code{summary.ipd} containing:
#'   \describe{
#'     \item{call}{The model formula.}
#'     \item{coefficients}{A glm-style table of estimates, SE, z, p.}
#'     \item{method}{Which IPD method was used.}
#'     \item{model}{Which downstream model was fitted.}
#'     \item{intercept}{Logical; whether an intercept was included.}
#'   }
#'
#' @export

summary.ipd <- function(object, ...) {

    if (!is(object, "ipd")) stop("Object is not of class 'ipd'")

    coef_tab <- object@coefTable

    summ <- list(
        call         = object@formula,
        coefficients = coef_tab,
        method       = object@method,
        model        = object@model,
        intercept    = object@intercept
    )

    class(summ) <- "summary.ipd"

    summ
}

#--- S3 PRINT.SUMMARY.IPD ------------------------------------------------------

#' Print summary.ipd
#'
#' @param x   An object of class \code{summary.ipd}.
#'
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export

print.summary.ipd <- function(x, ...) {

    if (!inherits(x, "summary.ipd")) {

        stop("Object is not of class 'summary.ipd'")
    }

    cat("\nCall:\n")
    cat(" ", deparse(x$call), "\n\n")
    cat("Method:   ", x$method, "\n")
    cat("Model:    ", x$model, "\n")
    cat("Intercept:", ifelse(x$intercept, "Yes", "No"), "\n\n")
    cat("Coefficients:\n")

    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)

    invisible(x)
}

#--- S3 TIDY.IPD ---------------------------------------------------------------

#' tidy re-exported from generics packages
#'
#' @return A wrapper for the \code{tidy} generic.
#' See \code{\link[generics]{tidy}} for details.
#'
#' @seealso
#' \code{\link[generics]{tidy}}
#'
#' @importFrom generics tidy
#'
#' @examples
#'
#' dat <- simdat()
#'
#' fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' tidy(fit)
#'
#' @export

generics::tidy

#' Tidy an ipd fit
#'
#' @param x   An object of class \code{ipd}.
#'
#' @param ... Ignored.
#'
#' @return A \code{\link[tibble]{tibble}} with columns
#'   \code{term, estimate, std.error, conf.low, conf.high}.
#'
#' @importFrom tibble tibble
#'
#' @examples
#'
#' dat <- simdat()
#'
#' fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' tidy(fit)
#'
#' @export

tidy.ipd <- function(x, ...) {

    if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

    tibble::tibble(
        term      = names(x@coefficients),
        estimate  = x@coefficients,
        std.error = x@se,
        conf.low  = x@ci[ , "lower"],
        conf.high = x@ci[ , "upper"]
    )
}

#--- S3 GLANCE.IPD -------------------------------------------------------------

#' glance re-exported from generics packages
#'
#' @return A wrapper for the \code{glance} generic.
#' See \code{\link[generics]{glance}} for details.
#'
#' @seealso
#' \code{\link[generics]{glance}}
#'
#' @importFrom generics glance
#'
#' @examples
#'
#' dat <- simdat()
#'
#' fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' glance(fit)
#'
#' @export

generics::glance

#' Glance at an ipd fit
#'
#' @param x   An object of class \code{ipd}.
#' @param ... Ignored.
#'
#' @return A one-row \code{\link[tibble]{tibble}} summarizing the fit.
#'
#' @importFrom tibble tibble
#'
#' @examples
#'
#' dat <- simdat()
#'
#' fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' glance(fit)
#'
#' @export

glance.ipd <- function(x, ...) {

    if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

    tibble::tibble(
        method           = x@method,
        model            = x@model,
        intercept        = x@intercept,
        nobs_labeled     = nrow(x@data_l),
        nobs_unlabeled   = nrow(x@data_u),
        call             = deparse(x@formula)
    )
}

#--- S3 AUGMENT.IPD ------------------------------------------------------------

#' augment re-exported from generics packages
#'
#' @return A wrapper for the \code{augment} generic.
#' See \code{\link[generics]{augment}} for details.
#'
#' @seealso
#' \code{\link[generics]{augment}}
#'
#' @importFrom generics augment
#'
#' @examples
#'
#' dat <- simdat()
#'
#' fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' augmented_df <- augment(fit)
#'
#' head(augmented_df)
#'
#' @export

generics::augment

#' Augment data from an ipd fit
#'
#' @param x    An object of class \code{ipd}.
#' @param data A \code{data.frame} to augment; defaults to \code{x@data_u}.
#' @param ...  Ignored.
#'
#' @return The \code{data.frame} with columns \code{.fitted} and \code{.resid}.
#'
#' @examples
#'
#' dat <- simdat()
#'
#' fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",
#'
#'     data = dat, label = "set_label")
#'
#' augmented_df <- augment(fit)
#'
#' head(augmented_df)
#'
#' @export

augment.ipd <- function(x, data = x@data_u, ...) {

    if (!inherits(x, "ipd")) stop("Object is not of class 'ipd'")

    data_aug <- data

    mm <- model.matrix(x@formula, data_aug)

    fv <- as.vector(mm %*% x@coefficients)

    data_aug$.fitted <- fv

    resp <- all.vars(x@formula)[1]

    data_aug$.resid  <- data_aug[[resp]] - fv

    data_aug
}
