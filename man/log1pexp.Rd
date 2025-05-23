% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_ppi.R
\name{log1pexp}
\alias{log1pexp}
\title{Log1p Exponential}
\usage{
log1pexp(x)
}
\arguments{
\item{x}{(vector): A numeric vector of inputs.}
}
\value{
(vector): A numeric vector where each element is the result of
log(1 + exp(x)).
}
\description{
Computes the natural logarithm of 1 plus the exponential of the input,
to handle large inputs.
}
