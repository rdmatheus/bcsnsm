% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_bcla.R
\name{lla}
\alias{lla}
\alias{dlla}
\alias{plla}
\alias{qlla}
\alias{rlla}
\title{The Log-Laplace Distribution}
\usage{
dlla(x, mu, sigma, log = FALSE, ...)

plla(q, mu, sigma, lower.tail = TRUE, ...)

qlla(p, mu, sigma, lower.tail = TRUE, ...)

rlla(n, mu, sigma)
}
\arguments{
\item{x, q}{vector of positive quantiles.}

\item{mu}{vector of strictly positive scale parameters.}

\item{sigma}{vector of strictly positive relative dispersion parameters.}

\item{log}{logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.}

\item{...}{further arguments.}

\item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
\code{P[X <= x]}, otherwise, \code{P[X > x]}.}

\item{p}{vector of probabilities.}

\item{n}{number of random values to return.}
}
\value{
\code{dlla} returns the density, \code{plla} gives the distribution
    function, \code{qlla} gives the quantile function, and \code{rlla}
    generates random observations.

    Invalid arguments will result in return value \code{NaN}.

    The length of the result is determined by \code{n} for \code{rlla}, and is the
    maximum of the lengths of the numerical arguments for the other functions.
}
\description{
Density, distribution function, quantile function, and random
    generation for the log-Laplace distribution with parameters \code{mu} and
    \code{sigma}.
}
\details{
A random variable X has a log-Laplace distribution with parameter \code{mu} and
    \code{sigma} if log(X) follows a Laplace distribution with location parameter \code{log(mu)}
    and dispersion parameter \code{sigma}. It can be showed that \code{mu} is the median of X.
}
\examples{
mu <- 5
sigma <- 0.5

# Sample generation
x <- rlla(500, mu, sigma)

# Density
hist(x, prob = TRUE, main = "The Log-Laplace Distribution", col = "white")
curve(dlla(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
legend("topright", "Probability density function", col = 2, lwd = 2, lty = 1)

# Distribution function
plot(ecdf(x), main = "The Log-Laplace Distribution", ylab = "Distribution function")
curve(plla(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
legend("bottomright", c("Emp. distribution function", "Theo. distribution function"),
  col = c(1, 2), lwd = 2, lty = c(1, 1)
)

# Quantile function
plot(seq(0.01, 0.99, 0.001), quantile(x, seq(0.01, 0.99, 0.001)),
  type = "l",
  xlab = "p", ylab = "Quantile function", main = "The Log-Laplace Distribution"
)
curve(qlla(x, mu, sigma), add = TRUE, col = 2, lwd = 2, from = 0, to = 1)
legend("topleft", c("Emp. quantile function", "Theo. quantile function"),
  col = c(1, 2), lwd = 2, lty = c(1, 1)
)
}
\references{
Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
    properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
}
\author{
Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
}
