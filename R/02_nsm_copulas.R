# Copula information -------------------------------------------------------------------------------
make_copula <- function(copula, delta) {

  switch(copula,

         ## Multivariate normal distribution (mu = (0, ..., 0)') -----------------------------------
         gaussian = {

           ### Marginal specifications
           dPSI <- function(x, log = FALSE) stats::dnorm(x, log = log)
           pPSI <- function(q) stats::pnorm(q)
           qPSI <- function(p) stats::qnorm(p)

           ### Density generating function
           dgf <- function(u, d, log = FALSE) {
             out <- -0.5 * u - (d/2) * log(2 * pi)

             if (!log) out <- exp(out)

             out
           }

           ### Joint density
           dmv <- function(x, Gamma) {

             mvnfast::dmvn(x, rep(0, ncol(Gamma)), Gamma, log = FALSE, ncores = 2)

           }

           ### Random generator
           rmv <- function(n, Gamma) {

             d <- ncol(Gamma)
             A <- Rfast::cholesky(Gamma)
             Z <- matrix(stats::rnorm(n * d), ncol = d)

             Z%*%A
           }

           ### Mahalanobis distance distribution function
           maha <- function(q, d) stats::pchisq(q, d)

         },

         ## Multivariate Student's t distribution (mu = (0, ..., 0)') ------------------------------
         t = {

           ### Marginal specifications
           dPSI <- function(x, log = FALSE) stats::dt(x, delta, log = log)
           pPSI <- function(q) stats::pt(q, delta)
           qPSI <- function(p) stats::qt(p, delta)

           ### Density generating function
           dgf <- function(u, d, log = FALSE){

             out <- lgamma(0.5 * (delta + d)) - 0.5 * (delta + d) * log(1 + u / delta) -
               lgamma(delta / 2) - (d/2) * log(delta * pi)

             if (!log) out <- exp(out)

             out
           }

           ### Joint density
           dmv <- function(x, Gamma) {

             mvnfast::dmvt(x, rep(0, ncol(Gamma)), Gamma, log = FALSE, df = delta, ncores = 2)

           }

           ### Random generator
           rmv <- function(n, Gamma) {

             d <- ncol(Gamma)
             A <- Rfast::cholesky(Gamma)
             Z <- matrix(stats::rnorm(n * d), ncol = d)

             Z%*%A / sqrt(stats::rgamma(n, delta / 2, delta / 2))

           }

           ### Mahalanobis distances distribution function
           maha <- function(q, d) stats::pf(q / d, d, delta)

         },

         ## Multivariate slash distribution (mu = (0, ..., 0)') ------------------------------------
         slash = {

           ### Marginal specifications
           Wsl <- distr::AbscontDistribution(
             d = function(x) dslash(x, nu = delta),
             Symmetry = distr::SphericalSymmetry(0)
           )

           dPSI <- function(x, log = FALSE) dslash(x, nu = delta, log = log)
           pPSI <- function(q) distr::p(Wsl)(q)
           qPSI <- function(p) distr::q(Wsl)(p)

           ### Density generating function
           dgf <- function(u, d, log = FALSE){

             id1 <- which(u == 0)
             id2 <- which(u != 0)

             out <- vector("numeric", length(u))

             out[id1] <- log(2 * delta) - (d/2) * log(2 * pi) - log(2 * delta + d)
             out[id2] <- log(delta) + delta * log(2) + log(ig(delta + d / 2, u[id2] / 2)) -
               (d/2) * log(pi) - (delta + d/2) * log(u[id2])

             if (!log) out <- exp(out)

             out
           }

           ### Joint density
           dmv <- function(x, Gamma) {

               delta <- 2 * delta

               if (is.vector(x))
                 x <- matrix(x, ncol = length(x))

               d <- ncol(Gamma)

               dec <- tryCatch(Rfast::cholesky(Gamma), error = function(e) e)
               tmp <- backsolve(dec, t(x), transpose = TRUE)
               rss <- colSums(tmp^2)

               id1 <- which(rss == 0)
               id2 <- which(rss != 0)

               out <- vector("numeric", length(rss))

               out[id1] <- -sum(log(diag(dec))) + log(delta) - (d/2) * log(2 * pi) - log(delta + d)
               out[id2] <- -sum(log(diag(dec))) + log(delta) + (delta / 2 - 1) * log(2) +
                 lgamma(0.5 * (delta + d)) + stats::pgamma(rss[id2] / 2, 0.5 * (delta + d), scale = 1, log.p = TRUE) -
                 (d/2) * log(pi) - (0.5 * (delta + d)) * log(rss[id2])


               exp(out)

           }

           ### Random generator
           rmv <- function(n, Gamma) {

             d <- ncol(Gamma)
             A <- Rfast::cholesky(Gamma)
             Z <- matrix(stats::rnorm(n * d), ncol = d)

             Z%*%A / sqrt(stats::rbeta(n, delta, 1))

           }

           ### Mahalanobis distances distribution function
           maha <- function(q, d) {
             stats::pchisq(q, d) - 2^delta * gamma(delta + d / 2) * stats::pchisq(q, d + 2 * delta) /
               (q^delta * gamma(d / 2))
           }


         },

         ## Multivariate hyperbolic distribution (mu = (0, ..., 0)') -------------------------------
         hyp = {

           ### Marginal specifications
           Whp <- distr::AbscontDistribution(
             d = function(x) dhyp(x, nu = delta),
             Symmetry = distr::SphericalSymmetry(0)
           )

           dPSI <- function(x, log = FALSE) dhyp(x, nu = delta, log = log)
           pPSI <- function(q) distr::p(Whp)(q)
           qPSI <- function(p) distr::q(Whp)(p)

           ### Density generating function
           dgf <- function(u, d, log = FALSE){

             out <- (d - 1) * log(delta) + log(besselK(delta * sqrt(1 + u), nu = 1 - d/2)) -
               (d/2) * log(2 * pi) - log(besselK(delta, nu = 1)) -
               (d/2 - 1) * log(delta * sqrt(1 + u))

             if (!log) out <- exp(out)

             out
           }

           ### Joint density
           dmv <- function(x, Gamma) {

               if (is.vector(x))
                 x <- matrix(x, ncol = length(x))

               d <- ncol(Gamma)

               dec <- tryCatch(Rfast::cholesky(Gamma), error = function(e) e)
               tmp <- backsolve(dec, t(x), transpose = TRUE)
               rss <- colSums(tmp^2)

               out <- -sum(log(diag(dec))) + (d - 1) * log(delta) + log(besselK(delta * sqrt(1 + rss), nu = 1 - d/2)) -
                 (d/2) * log(2 * pi) - log(besselK(delta, nu = 1)) -
                 (d/2 - 1) * log(delta * sqrt(1 + rss))

               exp(out)
           }

           ### Random generator
           rmv <- function(n, Gamma) {

             d <- ncol(Gamma)
             A <- Rfast::cholesky(Gamma)
             Z <- matrix(stats::rnorm(n * d), ncol = d)

             Z%*%A * sqrt(GIGrvg::rgig(n, 1.0, 1.0, delta^2))

           }

           ### Mahalanobis distances distribution function
           maha <- function(q, d){


             aux_f <- function(x){

               integrand <- function(u){

                 exp((d / 2) * (log(x) - log(2)) - lgamma(d / 2) +
                       log(ghyp::pgig(1 / u, 1, 1, delta^2)) + (0.5 * d - 1) * log(u) -0.5 * u * x)

               }

               stats::integrate(integrand, .Machine$double.eps^(1/2), Inf)$val

             }

             apply(matrix(q), 1, aux_f)

           }


         },

         stop(gettextf("%s copula not recognised", sQuote(copula)), domain = NA))

  list(dPSI = dPSI, pPSI = pPSI, qPSI = qPSI, dgf = dgf, dmv = dmv, rmv = rmv, maha = maha)

}

