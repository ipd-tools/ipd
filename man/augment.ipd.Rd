% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{augment.ipd}
\alias{augment.ipd}
\title{Augment data from an ipd fit}
\usage{
\method{augment}{ipd}(x, data = x@data_u, ...)
}
\arguments{
\item{x}{An object of class \code{ipd}.}

\item{data}{A \code{data.frame} to augment; defaults to \code{x@data_u}.}

\item{...}{Ignored.}
}
\value{
The \code{data.frame} with columns \code{.fitted} and \code{.resid}.
}
\description{
Augment data from an ipd fit
}
\examples{

dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

augmented_df <- augment(fit)

head(augmented_df)

}
