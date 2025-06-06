% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_bcs-nsm_distributions.R
\name{mvplot}
\alias{mvplot}
\title{Visualize Multivariate Data}
\usage{
mvplot(y, method = c("kendall", "spearman", "pearson"), ...)
}
\arguments{
\item{y}{a matrix or a data frame of quantitative data with variables in columns.}

\item{method}{a character string indicating which association coefficient
is to be computed for the upper panels of the plot. It is computed
via [cor](stats::cor) function. One of \code{"kendall"} (default),
\code{"pearson"}, or \code{"spearman"} can be used.}

\item{...}{arguments to be passed to or from methods.}
}
\value{
A panel plot displaying histograms on the diagonal, scatterplots in the
lower panels, and Kendall's tau correlation in the upper panels.
}
\description{
This function takes a matrix of quantitative data and generates a panel plot.
The diagonal panels contains histograms for each variable, the lower panels
display pairwise scatter plots, and the upper panels show some association
measure, such as Kendall's tau, Spearman's rho, or Pearson correlation coefficient.
}
\examples{
data <- data.frame(x1 = rnorm(1000),
                   x2 = rgamma(1000, 2, 4),
                   x3 = rbeta(1000, 4, 1),
                   x5 = rlnorm(1000, 3, 0.5))

mvplot(data)
}
\author{
Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
}
