% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_ppi.R
\name{compute_cdf}
\alias{compute_cdf}
\title{Empirical CDF of the Data}
\usage{
compute_cdf(Y, grid, w = NULL)
}
\arguments{
\item{Y}{(matrix): n x 1 matrix of observed data.}

\item{grid}{(matrix): Grid of values to compute the CDF at.}

\item{w}{(vector, optional): n-vector of sample weights.}
}
\value{
(list): Empirical CDF and its standard deviation at the specified
grid points.
}
\description{
Computes the empirical CDF of the data.
}
