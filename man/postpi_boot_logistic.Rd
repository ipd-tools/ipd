% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postpi_boot_logistic.R
\name{postpi_boot_logistic}
\alias{postpi_boot_logistic}
\title{PostPI Logistic Regression (Bootstrap Correction)}
\usage{
postpi_boot_logistic(X_l, Y_l, f_l, X_u, f_u, nboot = 100, se_type = "par")
}
\arguments{
\item{X_l}{(matrix): n x p matrix of covariates in the labeled data.}

\item{Y_l}{(vector): n-vector of labeled outcomes.}

\item{f_l}{(vector): n-vector of predictions in the labeled data.}

\item{X_u}{(matrix): N x p matrix of covariates in the unlabeled data.}

\item{f_u}{(vector): N-vector of predictions in the unlabeled data.}

\item{nboot}{(integer): Number of bootstrap samples. Defaults to 100.}

\item{se_type}{(string): Which method to calculate the standard errors.
Options include "par" (parametric) or "npar" (nonparametric).
Defaults to "par".}
}
\value{
A list of outputs: estimate of inference model parameters and
corresponding standard error based on both parametric and non-parametric
bootstrap methods.
}
\description{
Helper function for PostPI logistic regression (bootstrap correction)
}
\details{
Methods for correcting inference based on outcomes predicted by machine
learning (Wang et al., 2020)
\url{https://www.pnas.org/doi/abs/10.1073/pnas.2001238117}
}
\examples{

dat <- simdat(model = "logistic")

form <- Y - f ~ X1

X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>
  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

postpi_boot_logistic(X_l, Y_l, f_l, X_u, f_u, nboot = 200)

}
