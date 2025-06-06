# The Box-Cox type II logistic distribution -------------------------------------------------------------------

#' @name bcloii
#'
#' @title The Box-Cox Type II Logistic Distribution
#'
#' @description Density, distribution function, quantile function, and random
#'     generation for the Box-Cox type II logistic distribution with parameters \code{mu},
#'     \code{sigma}, and \code{lambda}.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param lambda vector of real-valued skewness parameters. If \code{lambda = 0}, the Box-Cox
#'     type II logistic distribution reduces to the log-type II logistic distribution with parameters
#'     \code{mu} and \code{sigma} (see \code{\link{lloii}}).
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'     \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... further arguments.
#'
#' @return \code{dbcloii} returns the density, \code{pbcloii} gives the distribution function,
#'     \code{qbcloii} gives the quantile function, and \code{rbcloii} generates random observations.
#'
#'     Invalid arguments will result in return value \code{NaN}, with an warning.
#'
#'     The length of the result is determined by \code{n} for \code{rbcloii}, and is the
#'     maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'     applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
#'     properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' mu <- 8
#' sigma <- 1
#' lambda <- 2
#'
#' # Sample generation
#' x <- rbcloii(10000, mu, sigma, lambda)
#'
#' # Density
#' hist(x, prob = TRUE, main = "The Box-Cox Type II Logistic Distribution", col = "white")
#' curve(dbcloii(x, mu, sigma, lambda), add = TRUE, col = 2, lwd = 2)
#' legend("topright", "Probability density function", col = 2, lwd = 2, lty = 1)
#'
#' # Distribution function
#' plot(ecdf(x), main = "The Box-Cox Type II Logistic Distribution", ylab = "Distribution function")
#' curve(pbcloii(x, mu, sigma, lambda), add = TRUE, col = 2, lwd = 2)
#' legend("bottomright", c("Emp. distribution function", "Theo. distribution function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
#'
#' # Quantile function
#' plot(seq(0.01, 0.99, 0.001), quantile(x, seq(0.01, 0.99, 0.001)),
#'   type = "l",
#'   xlab = "p", ylab = "Quantile function", main = "The Box-Cox Type II Logistic Distribution"
#' )
#' curve(qbcloii(x, mu, sigma, lambda), add = TRUE, col = 2, lwd = 2, from = 0, to = 1)
#' legend("topleft", c("Emp. quantile function", "Theo. quantile function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
NULL

# Density
#' @rdname bcloii
#' @export
dbcloii <- function(x, mu, sigma, lambda, log = FALSE, ...) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma), length(lambda)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0)] <- NaN

  # Positive density index
  id1 <- which(x > 0 & lambda != 0 & !is.nan(pmf))
  id2 <- which(x > 0 & lambda == 0 & !is.nan(pmf))

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  z[id1] <- ((x[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(x[id2] / mu[id2]) / sigma[id2]

  pmf[id1] <- (lambda[id1] - 1) * log(x[id1]) + stats::dlogis(z[id1], log = TRUE) -
    stats::plogis(1 / (sigma[id1] * abs(lambda[id1])), log.p = TRUE) -
    lambda[id1] * log(mu[id1]) - log(sigma[id1])

  pmf[id2] <- stats::dlogis(z[id2], log = TRUE) - log(sigma[id2] * x[id2])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
#' @rdname bcloii
#' @export
pbcloii <- function(q, mu, sigma, lambda, lower.tail = TRUE, ...) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma), length(lambda)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  # Extended Box-Cox transformation
  z <- rep(NaN, length.out = maxl)

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda != 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda == 0)

  z[id1] <- ((q[id1] / mu[id1])^lambda[id1] - 1) / (sigma[id1] * lambda[id1])
  z[id2] <- log(q[id2] / mu[id2]) / sigma[id2]

  id1 <- which(q > 0 & mu > 0 & sigma > 0 & lambda <= 0)
  id2 <- which(q > 0 & mu > 0 & sigma > 0 & lambda > 0)

  cdf <- rep(NaN, length.out = maxl)
  cdf[id1] <- stats::plogis(z[id1]) / stats::plogis(1 / (sigma[id1] * abs(lambda[id1])))
  cdf[id2] <- (stats::plogis(z[id2]) - stats::plogis(-1 / (sigma[id2] * lambda[id2]))) /
    stats::plogis(1 / (sigma[id2] * lambda[id2]))

  cdf[which(q <= 0 & mu > 0 & sigma > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
#' @rdname bcloii
#' @export
qbcloii <- function(p, mu, sigma, lambda, lower.tail = TRUE, ...) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma), length(lambda)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)
  lambda <- rep(lambda, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda <= 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda > 0)

  zp[id1] <- stats::qlogis(p[id1] * stats::plogis(1 / (sigma[id1] * abs(lambda[id1]))))
  zp[id2] <- stats::qlogis(1 - (1 - p[id2]) * stats::plogis(1 / (sigma[id2] * abs(lambda[id2]))))

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda != 0)
  id2 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0 & lambda == 0)
  id3 <- which(p == 0 & mu > 0 & sigma > 0)
  id4 <- which(p == 1 & mu > 0 & sigma > 0)

  qtf[id1] <- exp(log(mu[id1]) + (1 / lambda[id1]) * log1p(sigma[id1] * lambda[id1] * zp[id1]))
  qtf[id2] <- exp(log(mu[id2]) + sigma[id2] * zp[id2])
  qtf[id3] <- 0
  qtf[id4] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
#' @rdname bcloii
#' @export
rbcloii <- function(n, mu, sigma, lambda) {
  u <- stats::runif(n)
  qbcloii(u, mu, sigma, lambda)
}

# BCS class
bcloii <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "bcloii"

  # Name
  out$name <- "Box-Cox Type II Logistic"

  # Number of parameters
  out$npar <- 3

  # Extra parameter
  out$extrap <- FALSE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCT(mu.link = "log"), trace = FALSE), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }

    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")), exp(stats::coef(gamlss_fit, "sigma")), stats::coef(gamlss_fit, "nu"))
    } else {
      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)
      c(stats::median(x), asinh(CV / 1.5) * stats::qlogis(0.75), 0L)
    }

  }

  structure(out, class = "bcs")
}




# The log type II logistic distribution ----------------------------------------------------------------------

#' @name lloii
#'
#' @title The Log-Type II Logistic Distribution
#'
#' @description Density, distribution function, quantile function, and random
#'     generation for the log-type II logistic distribution with parameters \code{mu} and
#'     \code{sigma}.
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of strictly positive scale parameters.
#' @param sigma vector of strictly positive relative dispersion parameters.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'     \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param log logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... further arguments.
#'
#' @details A random variable X has a log-type II logistic distribution with parameter \code{mu} and
#'     \code{sigma} if log(X) follows a type II logistic distribution with location parameter \code{log(mu)}
#'     and dispersion parameter \code{sigma}. It can be showed that \code{mu} is the median of X.
#'
#' @return \code{dlloii} returns the density, \code{plloii} gives the distribution
#'     function, \code{qlloii} gives the quantile function, and \code{rlloii}
#'     generates random observations.
#'
#'     Invalid arguments will result in return value \code{NaN}.
#'
#'     The length of the result is determined by \code{n} for \code{rlloii}, and is the
#'     maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical
#'     properties and parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' mu <- 5
#' sigma <- 0.2
#'
#' # Sample generation
#' x <- rlloii(1000, mu, sigma)
#'
#' # Density
#' hist(x, prob = TRUE, main = "The Log-Type II Logistic Distribution", col = "white")
#' curve(dlloii(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
#' legend("topright", "Probability density function", col = 2, lwd = 2, lty = 1)
#'
#' # Distribution function
#' plot(ecdf(x), main = "The Log-Type II Logistic Distribution", ylab = "Distribution function")
#' curve(plloii(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
#' legend("bottomright", c("Emp. distribution function", "Theo. distribution function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
#'
#' # Quantile function
#' plot(seq(0.01, 0.99, 0.001), quantile(x, seq(0.01, 0.99, 0.001)),
#'   type = "l",
#'   xlab = "p", ylab = "Quantile function", main = "The Log-Type II Logistic Distribution"
#' )
#' curve(qlloii(x, mu, sigma), add = TRUE, col = 2, lwd = 2, from = 0, to = 1)
#' legend("topleft", c("Emp. quantile function", "Theo. quantile function"),
#'   col = c(1, 2), lwd = 2, lty = c(1, 1)
#' )
NULL

# Density
#' @rdname lloii
#' @export
dlloii <- function(x, mu, sigma, log = FALSE, ...) {
  if (is.matrix(x)) d <- ncol(x) else d <- 1L

  maxl <- max(c(length(x), length(mu), length(sigma)))

  x <- rep(x, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  pmf <- rep(-Inf, maxl)

  # NaN index
  pmf[which(mu <= 0 | sigma <= 0)] <- NaN

  # Positive density index
  id <- which(x > 0 & !is.nan(pmf))

  # Transformed variables
  z <- rep(NaN, length.out = maxl)

  z[id] <- log(x[id] / mu[id]) / sigma[id]

  pmf[id] <- stats::dlogis(z[id], log = TRUE) - log(sigma[id] * x[id])

  if (!log) pmf <- exp(pmf)

  if (d > 1L) matrix(pmf, ncol = d) else pmf
}

## Distribution function
#' @rdname lloii
#' @export
plloii <- function(q, mu, sigma, lower.tail = TRUE, ...) {
  if (is.matrix(q)) d <- ncol(q) else d <- 1L

  maxl <- max(c(length(q), length(mu), length(sigma)))

  q <- rep(q, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  # Transformed variables
  z <- rep(NaN, length.out = maxl)

  id <- which(q > 0 & mu > 0 & sigma > 0)

  z[id] <- log(q[id] / mu[id]) / sigma[id]

  cdf <- rep(NaN, length.out = maxl)
  cdf[id] <- stats::plogis(z[id])

  cdf[which(q <= 0 & mu > 0 & sigma > 0)] <- 0

  if (!lower.tail) cdf <- 1 - cdf

  if (d > 1L) matrix(cdf, ncol = d) else cdf
}

## Quantile function
#' @rdname lloii
#' @export
qlloii <- function(p, mu, sigma, lower.tail = TRUE, ...) {
  if (is.matrix(p)) d <- ncol(p) else d <- 1L

  maxl <- max(c(length(p), length(mu), length(sigma)))

  p <- rep(p, length.out = maxl)
  mu <- rep(mu, length.out = maxl)
  sigma <- rep(sigma, length.out = maxl)

  if (!lower.tail) p <- 1 - p

  qtf <- zp <- rep(NaN, length.out = maxl)

  # z_p
  id <- which(p > 0 & p < 1 & mu > 0 & sigma > 0)
  zp[id] <- stats::qlogis(p[id])

  # Quantile function
  id1 <- which(p > 0 & p < 1 & mu > 0 & sigma > 0)
  id2 <- which(p == 0 & mu > 0 & sigma > 0)
  id3 <- which(p == 1 & mu > 0 & sigma > 0)

  qtf[id1] <- exp(log(mu[id1]) + sigma[id1] * zp[id1])
  qtf[id2] <- 0
  qtf[id3] <- Inf

  if (d > 1L) matrix(qtf, ncol = d) else qtf
}

# Random generation
#' @rdname lloii
#' @export
rlloii <- function(n, mu, sigma) {
  exp(stats::rlogis(n, log(mu), sigma))
}

# BCS class
lloii <- function(x) {
  out <- list()

  # Abbreviation
  out$abb <- "lloii"

  # Name
  out$name <- "Log-Type II Logistic"

  # Number of parameters
  out$npar <- 2

  # Extra parameter
  out$extrap <- FALSE

  # Initial values -------------------------------------------------------------
  out$start <- function(x) {

    n <- length(x)

    gamlss_fit <- suppressWarnings(try(gamlss::gamlss(x ~ 1, family = gamlss.dist::BCT(mu.link = "log"),
                                                      trace = FALSE, nu.fix = TRUE, nu.start = 0L), silent = TRUE))

    if (unique(grepl("Error", gamlss_fit))) {
      convergence <- FALSE
    } else {
      convergence <- gamlss_fit$converged
    }


    if (convergence) {
      c(exp(stats::coef(gamlss_fit, "mu")), exp(stats::coef(gamlss_fit, "sigma")))
    } else {
      CV <- 0.75 * (stats::quantile(x, 0.75) - stats::quantile(x, 0.25)) / stats::median(x)
      c(stats::median(x), asinh(CV / 1.5) * stats::qlogis(0.75))
    }


  }

  structure(out, class = "bcs")
}
