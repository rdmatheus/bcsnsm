% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/08_regression_methods.R
\name{residuals.bcsnsmreg}
\alias{residuals.bcsnsmreg}
\title{Extract Model Residuals for a Fitted BCS-NSM Regression}
\usage{
\method{residuals}{bcsnsmreg}(
  object,
  type = c("quantile", "marginal", "mahalanobis", "epsilon"),
  ...
)
}
\arguments{
\item{object}{an \code{"bcsnsmreg"} object.}

\item{type}{character; specifies which residual should be extracted. The available arguments are
\code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and \code{"epsilon"}.}

\item{...}{further arguments passed to or from other methods.}
}
\description{
Extract Model Residuals for a Fitted BCS-NSM Regression
}
\author{
Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
}
