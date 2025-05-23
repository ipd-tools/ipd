% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_ipd.R
\name{.parse_inputs}
\alias{.parse_inputs}
\title{Validate and split input data}
\usage{
.parse_inputs(data, label = NULL, unlabeled_data = NULL, na_action = "na.fail")
}
\arguments{
\item{data}{data.frame or coercible object}

\item{label}{optional column flagging labeled vs. unlabeled rows}

\item{unlabeled_data}{optional data.frame of unlabeled observations}

\item{na_action}{how to handle missing data: "na.fail" or "na.omit"}
}
\value{
A list with components data_l and data_u (data.frames)
}
\description{
Validate and split input data
}
\keyword{internal}
