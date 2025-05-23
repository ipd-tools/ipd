% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_ppi.R
\name{zconfint_generic}
\alias{zconfint_generic}
\title{Normal Confidence Intervals}
\usage{
zconfint_generic(mean, std_mean, alpha, alternative)
}
\arguments{
\item{mean}{(float): Estimated normal mean.}

\item{std_mean}{(float): Estimated standard error of the mean.}

\item{alpha}{(float): Significance level in [0,1]}

\item{alternative}{(string): Alternative hypothesis, either 'two-sided',
'larger' or 'smaller'.}
}
\value{
(vector): Lower and upper (1 - alpha) * 100\% confidence limits.
}
\description{
Calculates normal confidence intervals for a given alternative at a given
significance level.
}
