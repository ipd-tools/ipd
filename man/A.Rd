% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_pspa.R
\name{A}
\alias{A}
\title{Calculation of the matrix A based on single dataset}
\usage{
A(
  X,
  Y,
  quant = NA,
  theta,
  method = c("ols", "quantile", "mean", "logistic", "poisson")
)
}
\arguments{
\item{X}{Array or data.frame containing covariates}

\item{Y}{Array or data.frame of outcomes}

\item{quant}{quantile for quantile estimation}

\item{theta}{parameter theta}

\item{method}{indicates the method to be used for M-estimation.
Options include "mean", "quantile", "ols", "logistic", and "poisson".}
}
\value{
matrix A based on single dataset
}
\description{
\code{A} function for the calculation of the matrix A based on single dataset
}
