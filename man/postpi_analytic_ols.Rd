% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postpi_analytic_ols.R
\name{postpi_analytic_ols}
\alias{postpi_analytic_ols}
\title{PostPI OLS (Analytic Correction)}
\usage{
postpi_analytic_ols(X_l, Y_l, f_l, X_u, f_u, original = FALSE)
}
\arguments{
\item{X_l}{(matrix): n x p matrix of covariates in the labeled data.}

\item{Y_l}{(vector): n-vector of labeled outcomes.}

\item{f_l}{(vector): n-vector of predictions in the labeled data.}

\item{X_u}{(matrix): N x p matrix of covariates in the unlabeled data.}

\item{f_u}{(vector): N-vector of predictions in the unlabeled data.}

\item{original}{(boolean): Logical argument to use original method from
Wang et al. (2020). Defaults to FALSE; TRUE retained for posterity.}
}
\value{
A list of outputs: estimate of the inference model parameters and
corresponding standard error estimate.
}
\description{
Helper function for PostPI OLS estimation (analytic correction)
}
\details{
Methods for correcting inference based on outcomes predicted by machine
learning (Wang et al., 2020)
\url{https://www.pnas.org/doi/abs/10.1073/pnas.2001238117}
}
\examples{

dat <- simdat(model = "ols")

form <- Y - f ~ X1

X_l <- model.matrix(form, data = dat[dat$set_label == "labeled", ])

Y_l <- dat[dat$set_label == "labeled", all.vars(form)[1]] |>
  matrix(ncol = 1)

f_l <- dat[dat$set_label == "labeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

X_u <- model.matrix(form, data = dat[dat$set_label == "unlabeled", ])

f_u <- dat[dat$set_label == "unlabeled", all.vars(form)[2]] |>
  matrix(ncol = 1)

postpi_analytic_ols(X_l, Y_l, f_l, X_u, f_u)

}
