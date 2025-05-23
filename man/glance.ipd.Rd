% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{glance.ipd}
\alias{glance.ipd}
\title{Glance at an ipd fit}
\usage{
\method{glance}{ipd}(x, ...)
}
\arguments{
\item{x}{An object of class \code{ipd}.}

\item{...}{Ignored.}
}
\value{
A one-row \code{\link[tibble]{tibble}} summarizing the fit.
}
\description{
Glance at an ipd fit
}
\examples{

dat <- simdat()

fit <- ipd(Y - f ~ X1, method = "pspa", model = "ols",

    data = dat, label = "set_label")

glance(fit)

}
