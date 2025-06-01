# Additional methods --------------------------------------------------------------------------
#' @name bcsnsm-methods
#' @title Methods for BCS-NSM Multivariate Regression Fits
#' 
#' @description These are S3 methods for extracting information from 
#'     fitted objects of class \code{"bcsnsm"}, which represent multivariate 
#'     regression models based on the BCS-NSM framework. The available methods 
#'     include functions to extract estimated coefficients, standard errors, 
#'     variance-covariance matrices, log-likelihood values, information criteria,
#'     and the original model call.
#' 
#' @param object an object of class \code{"bcsnsm"}, a result of a call to \code{\link{bcsnsm}}.
#' @param parm a character; specifies which set of parameters should be used to extract the
#'     asymptotic covariance matrix of the maximum likelihood estimators (for \code{vcov} function). The options are 
#'     \code{"beta"} (default), \code{"sigma"}, \code{"lambda"}, \code{"nu"}, 
#'     \code{"gamma"}, and \code{"all"}, where all parameters are considered
#' @param formula a model \link{formula} or \link{terms} object or an \code{"sdlrm"} object.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical Akaike information criteria (AIC).
#' @param ... additional argument(s) for methods.
#'
#' @returns
#' \itemize{
#' \item \code{model.frame} returns a \code{data.frame} containing the variables required
#'     by \code{formula} and any additional arguments provided via \code{...}.
#' \item \code{model.matrix} returns a list, where each element consists of the corresponding 
#'     design matrix of that marginal distribution used in the regression structure.
#' \item \code{coef} returns a list, each element of which is a numeric vector with the 
#'     estimated regression coefficients for each margin.
#' \item \code{vcov} returns the asymptotic covariance matrix of the regression coefficients,
#'     based on the \code{parm} argument.
#' \item \code{logLik} returns the log-likelihood value of the fitted model.
#' \item \code{AIC} returns a numeric value representing the generalized information 
#'     criterion, such as AIC or BIC, depending on the value of \code{k}.
#'     The default \code{k = 2} corresponds to AIC, while \code{k = log(n)} gives BIC.
#' }
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#' 
#' @examples
#' # Data set: macronutrients (for description run ?macronutrients)
#' # Consider modeling the animal and plant protein regarding income level
#' 
#' ## Response distribution:
#' mvplot(macronutrients[, 1:2])
#' 
#' ## Animal and plant protein by income level:
#' boxplot(animal_p ~ income, macronutrients)
#' boxplot(plant_p ~ income, macronutrients)
#' 
#' ## Fit the model with Gaussian copula and BCNO margins
#' fit <- bcsnsm(animal_p + plant_p ~ income | income, data = macronutrients)
#' 
#' # Coef
#' coef(fit)
#' 
#' # vcov
#' vcov(fit)
#' vcov(fit, parm = "sigma")
#' vcov(fit, parm = "lambda")
#' vcov(fit, parm = "gamma")
#' vcov(fit, parm = "all")
#' 
#' # Log-likelihood value
#' logLik(fit)
#' 
#' # AIC and BIC
#' AIC(fit)
#' AIC(fit, k = log(fit$nobs))
#' 
#' # Model matrices
#' model.matrix(fit)      ## Model matrices for all margins
#' model.matrix(fit)[[1]] ## Model matrix for the first margin
#' model.matrix(fit)[[2]] ## Model matrix for the second margin

## Model frame
#' @export
#' @rdname bcsnsm-methods
model.frame.bcsnsm <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname bcsnsm-methods
model.matrix.bcsnsm <- function(object, ...) {
  f <- object$formula
  mf <- stats::model.frame(object)
  Map(function(x) stats::model.matrix(object = f, data = mf, rhs = x), 1:object$d)
}

## Model matrix
#' @export
#' @rdname bcsnsm-methods
coef.bcsnsm <- function(object, ...) {
  object$coefficients
}

#  Variance-covariance matrix
#' @rdname bcsnsm-methods
#' @export
vcov.bcsnsm <- function(object, parm = c("beta", "sigma", "lambda", "nu", "gamma", "all"), ...) {
  
  parm <- match.arg(parm, c("beta", "sigma", "lambda", "nu", "gamma", "all"))
  covm <- object$vcov

  y <- object$y
  d <- object$d
  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  ## Beta
  beta_id <- object$optim_pars$par_id$beta
  vcov_beta <- Map(function(id) covm[id, id], beta_id)
  vcov_beta <- stats::setNames(vcov_beta, colnames(y))
  vcov_beta <- stats::setNames(vcov_beta, colnames(y))

  ## Sigma
  vcov_sigma <- covm[object$optim_pars$par_id$sigma, object$optim_pars$par_id$sigma]
  colnames(vcov_sigma) <- rownames(vcov_sigma) <- colnames(y)

  ## Lambda
  vcov_lambda <- as.matrix(covm[object$optim_pars$par_id$lambda, object$optim_pars$par_id$lambda])
  colnames(vcov_lambda) <- rownames(vcov_lambda) <- colnames(y)[lambda_id]

  ## Nu
  vcov_nu <- as.matrix(covm[object$optim_pars$par_id$nu, object$optim_pars$par_id$nu])
  colnames(vcov_nu) <- rownames(vcov_nu) <- colnames(y)[nu_id]

  ## Gamma
  vcov_gamma <- as.matrix(covm[object$optim_pars$par_id$gamma, object$optim_pars$par_id$gamma])

  if (object$association == "uniform") {
    names_gamma <- "gamma"
  } else if (object$association == "non-associative") {
    names_gamma <- NULL
  } else if (object$association == "unstructured") {
    names_gamma <- vector()
    id <- matrix(which(upper.tri(diag(object$d)), arr.ind = TRUE)[order(which(upper.tri(diag(object$d)), arr.ind = TRUE)[, 1]), ], ncol = 2)
    for (i in 1:nrow(id)) {
      names_gamma[i] <- paste0("gamma", id[i, 1], id[i, 2])
    }
  }

  colnames(vcov_gamma) <- rownames(vcov_gamma) <- names_gamma
  
  switch(parm,
         "beta" = vcov_beta,
         "sigma" = vcov_sigma,
         "lambda" = vcov_lambda,
         "nu" = vcov_nu,
         "gamma" = vcov_gamma,
         "all" = covm
  )
}


# Log-likelihood
#' @rdname bcsnsm-methods
#' @export
logLik.bcsnsm <- function(object, ...) {
  structure(object$logLik,
            df = length(object$optim_pars$par) + as.numeric(!is.null(object$eta)),
            class = "logLik")
}

# AIC
#' @export
#' @rdname bcsnsm-methods
AIC.bcsnsm <- function(object, ..., k = 2) {
  npars <- length(object$optim_pars$par) + as.numeric(!is.null(object$eta))
  AIC <- -2 * object$logLik + k * npars

  class(AIC) <- "AIC"
  return(AIC)
}


# Residuals -----------------------------------------------------------------------------------
#' Extract Residuals from BCS-NSM Multivariate Regression Fits
#' 
#' @description
#' Extracts several types of residuals for assessing the goodness-of-fit of BCS-NSM multivariate 
#'     regression models. Available residuals include: (i) marginal quantile residuals for 
#'     evaluating the adequacy of each marginal distribution; (ii) epsilon transformations for 
#'     diagnostic analysis under the copula structure; (iii) squared Mahalanobis distances of the 
#'     epsilon transformations to detect multivariate outliers; and (iv) overall quantile residuals 
#'     based on the distribution of the Mahalanobis distances.
#'
#' @param object an object of class \code{"bcsnsm"}, a result of a call to \code{\link{bcsnsm}}.
#' @param type character; specifies which residual should be extracted. The available options are
#'     \code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and \code{"epsilon"}.
#'     See details below.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' The fit of a BCS-NSM multivariate regression has some assumptions that need to be evaluated
#'     with a diagnostic analysis. Medeiros and Ferrari (2025), using the properties of normal
#'     scale mixture copulas, proposed diagnostic tools that can be used to assess the 
#'     goodness-of-fit in practical applications. 
#' 
#' The simplest is based on obtaining the quantile residuals (Dunn and Smith, 1996) for each 
#'     Box-Cox symmetric marginal distribution, which is useful for assessing goodness-of-fit at 
#'     the margins of the postulated model. These residuals are called marginal quantile residuals,
#'     and they can be accessed with the function \code{residuals.bcsnsm} with the 
#'     argument \code{type = "marginal"}. If the \eqn{j}th BCS marginal distribution is well-fitted
#'     to the model, the marginal quantile residuals are expected to have an approximately standard 
#'     normal distribution.
#' 
#' The other options available are based on the relationship that the BCS-NSM distributions have 
#'     with the normal scale mixtures (NSM) distributions. Let 
#'     \eqn{\bold{Y}_i = (Y_{i1}, \ldots,Y_{id})^\top}, for \eqn{i \in \{1, \ldots, n\}}, be a 
#'     sample of random vectors under the specification of a multivariate BCS-NSM regression, with 
#'     association matrix \eqn{\bold{\Gamma}}.
#' 
#' The BCS-NSM distributions are related to the NSM distributions through epsilon transformations:
#'   
#'   \eqn{
#'     \varepsilon_{ij} = \Psi_{H}^{-1}\left\{F_j(Y_{ij})\right\},
#'   }
#' 
#' where \eqn{\Psi_H} denotes the cumulative distribution function (CDF) of a standard univariate 
#'     NSM distribution and \eqn{F_j} denotes the CDF of the \eqn{j}th marginal BCS distribution. 
#'     The NSM distributions are a class of multivariate distributions whose elements are 
#'     characterized by the CDF \eqn{H}, called the mixing CDF. 
#' 
#' If \eqn{\bold{\varepsilon}_j = (\varepsilon_{1j}, \ldots, \varepsilon_{nj})^\top}, it can be 
#'     shown that the vectors \eqn{\bold{\varepsilon_{1}}, \ldots, \bold{\varepsilon_{d}}} form 
#'     a \eqn{d}-variate random sample of the NSM distribution with location vector 
#'     \eqn{(0, \ldots, 0)}, dispersion matrix \eqn{\bold{\Gamma}} (the same association matrix 
#'     as the corresponding BCS-NSM distribution), and mixing CDF \eqn{H}. 
#'   
#'   
#'   For instance, if a BCS-NSM regression with Gaussian copula is assumed, the epsilon 
#'       transformations are expected to have a multivariate normal distribution approximately. 
#'       On the other hand, if the \eqn{t} copula is assumed, the transformations are expected to 
#'       have a multivariate \eqn{t} distribution, and so on.
#'   
#'   The matrix with the epsilon transformations can be accessed with the function 
#'       \code{residuals.bcsnsm} with the argument \code{type = "epsilon"}. The squared 
#'       Mahalanobis distances of the epsilon transformations, i.e.,
#'   
#'   \eqn{
#'   M_i = \bold{\varepsilon}_i^\top \bold{\Gamma}^{-1} \bold{\varepsilon}_i, \quad i \in \{1, \ldots, n\}
#'   }
#' 
#'     are useful for identifying outliers in the multivariate context. The transformed Mahalanobis 
#'     distances can be accessed with the function \code{residuals.bcsnsm} with the argument 
#'     \code{type = "mahalanobis"}.
#' 
#' Furthermore, the CDF of the distances is known and can be written as (Lange and Sinsheimer, 1993)
#' 
#' \eqn{
#'   F_M(m) = P(M_i \leqslant m) = \dfrac{m^{d/2}}{2^{d/2} \Gamma(d/2)} \int_0^\infty H(1/u) u^{d/2-1} e^{-mu/2} \textrm{d}u, \quad m > 0.
#' }
#' 
#' This result allows us to define a general quantile residual to evaluate the overall 
#'     goodness-of-fit more simply than with the epsilon transformations. Medeiros and 
#'     Ferrari (2025) define the so-called overall quantile residuals as
#' 
#' \eqn{
#'   r_{i}^q = \Psi^{-1}[F_M(M_i)], \quad i \in \{1, \ldots, n\}
#' }
#' 
#' If the model is well-fit to the data, the overall quantile residuals are expected to have an 
#'     approximately standard normal distribution. The overall quantile residuals can be accessed 
#'     with the function \code{residuals.bcsnsm} with the argument \code{type = "quantile"} (default).
#' 
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#' 
#' @returns If \code{type = "quantile"} or \code{type = "mahalanobis"}, the \code{residual.bcsnsm}
#'     function returns a vector with the overall quantile residuals or the transformed Mahalanobis 
#'     (squared) distances, respectively. Conversely, if \code{type = "marginal"} or 
#'     \code{type = "epsilon"}, it returns a matrix with the quantile residuals or the epsilon 
#'     transforms of each marginal fit, respectively.
#'     
#' @references
#'  Dunn, P. K., and Smyth, G. K.(1996). Randomized quantile residuals. 
#'  \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 236–-24
#' 
#'  Ferrari, S. L., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, \bold{101}, 321--344.
#'      
#'  Lange, K., and Sinsheimer, J. S. (1993). Normal/independent distributions and their applications 
#'      in robust regression. \emph{Journal of Computational and Graphical Statistics}, \bold{2}, 175-–198      
#'
#'  Medeiros, R. M. R., and Ferrari, S. L. P. (2025). Multivariate Box-Cox symmetric regression 
#'      models generated by a normal scale mixture copula.
#'  
#'  Zeileis A., and Croissant, Y. (2010). Extended model formulas in
#'  R: multiple parts and multiple responses. \emph{Journal of Statistical Software},
#'  \bold{34}, 1–13.     
#'     
#' @examples
#' # Data set: macronutrients (for description run ?macronutrients)
#' # Consider modeling the animal and plant protein regarding income level
#' 
#' ## Response distribution:
#' mvplot(macronutrients[, 1:2])
#' 
#' ## Animal and plant protein by income level:
#' boxplot(animal_p ~ income, macronutrients)
#' boxplot(plant_p ~ income, macronutrients)
#' 
#' ## Fit the model with Gaussian copula and BCNO margins
#' fit <- bcsnsm(animal_p + plant_p ~ income | income, data = macronutrients)
#' 
#' ## Overall quantile residuals (default)
#' rq <- residuals(fit)
#' qqnorm(rq)
#' qqline(rq)
#' 
#' ## Marginal quantile residuals
#' rm <- residuals(fit, "marginal")
#' qqnorm(rm[, 1]); qqline(rm[, 1])
#' qqnorm(rm[, 2]); qqline(rm[, 2])
#' 
#' ## Transformed mahalanobis distances
#' maha <- residuals(fit, "mahalanobis")
#' plot(maha, type = "h")
#' 
#' ## Epsilon transformations
#' epsilon <- residuals(fit, "epsilon")
#' plot(epsilon, pch = 16)
#' qqnorm(epsilon[, 1]); qqline(epsilon[, 1])
#' qqnorm(epsilon[, 2]); qqline(epsilon[, 2])
#' @export
residuals.bcsnsm <- function(object, type = c("quantile", "marginal", "mahalanobis", "epsilon"), ...) {

  ## raw response residuals and desired type
  type <- match.arg(type, c("quantile", "marginal", "mahalanobis", "epsilon"))

  y <- object$y
  n <- object$nobs
  d <- object$d

  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  copula <- object$copula
  eta <- object$eta

  mcopula <- make_copula(copula, eta / (1 - eta))
  qPSI <- mcopula$qPSI
  maha <- mcopula$maha

  mu <- matrix(object$marginal.pars$mu, ncol = d)
  sigma <- object$marginal.pars$sigma
  lambda <- rep(NA, d)
  lambda[lambda_id] <- object$marginal.pars$lambda
  nu <- rep(NA, d)
  nu[nu_id] <- object$marginal.pars$nu

  EPS <- .Machine$double.eps

  epsilon <- rm <- matrix(NA, n, d)
  for(j in 1:d){

    epsilon[, j] <- qPSI(pmin(pmax(
      get(paste0("p", margins[j]), envir = asNamespace("bcsnsm"))(y[, j], mu = mu[, j],
                                                           sigma = sigma[j], lambda = lambda[j],
                                                           nu = nu[j]),
      EPS), 1 - EPS))


    rm[, j] <- stats::qnorm(pmin(pmax(
      get(paste0("p", margins[j]), envir = asNamespace("bcsnsm"))(as.matrix(y[, j]), mu = mu[, j],
                                                           sigma = sigma[j], lambda = lambda[j],
                                                           nu = nu[j]), EPS), 1 - EPS))

  }

  colnames(epsilon) <- colnames(rm) <- colnames(y)

  association <- if (object$association == "non-associative") "nonassociative" else object$association
  Gamma <- get(association, envir = asNamespace("bcsnsm"))(d)$Gamma(object$gamma)

  # Squared Mahalanobis distance
  mahalanobis <- mahalanobis(epsilon, rep(0L, d), Gamma)

  # Quantile residuals
  rq <- as.numeric(stats::qnorm(pmin(pmax(
    maha(mahalanobis, d), EPS), 1 - EPS)))

  res <- switch(type,
                "quantile" = rq,
                "marginal" = rm,
                "mahalanobis" = mahalanobis,
                "epsilon" = epsilon)

  res

}



# Summary method ------------------------------------------------------------------------------
#' @name summary.bcsnsm
#'
#' @title Summarizing a BCS-NSM Multivariate Regression Fit
#'
#' @description \code{summary} method for class \code{"bcsnsm"}.
#'
#' @param object an object of class \code{"bcsnsm"}, a result of a call to \code{\link{bcsnsm}}.
#' @param x an object of class \code{"bcsnsm.sdlrm"}, a result of a call to \code{summary.bcsnsm}.
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values.
#' @param margin an integer value specifying a marginal for which the summary will
#'     be presented. If missing (default), a summary for all marginal distributions 
#'     is considered.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The function \code{summary.bcsnsm} returns an object of class \code{"summary.bcsnsm"},
#'     which consists of a list with the following components:
#'  \describe{
#'     \item{call}{the original function call, given in \code{object}.}
#'     \item{beta_tab}{a list with summary statistics for the regression coefficients.}
#'     \item{sigma_tab}{summary statistics for the \code{sigma} marginal parameter vector.}
#'     \item{lambda_tab}{summary statistics for the \code{lambda} marginal parameter vector, if any.}
#'     \item{nu_tab}{summary statistics for the \code{nu} marginal parameter vector, if any.}
#'     \item{Gamma}{the fitted association matrix.}
#'     \item{margins}{a character vector with the specified marginal BCS distributions.}
#'     \item{links}{a character vector with the specified link functions.}
#'     \item{association}{structure specified for the association matrix. It can be one of
#'         \code{"non-associative"}, \code{"unstructured"}, or \code{"uniform"}.}
#'     \item{copula,eta}{\code{copula} is a character which informs which normal scale mixture 
#'         distribution was used to generate the NSM copula and \code{eta} is the 
#'         possible extra parameter associated with the copula, parameterized to take 
#'         values at (0, 1).}
#'     \item{residuals}{the overall quantile residuals.}
#'     \item{logLik}{log-likelihood value of the fitted model.}
#'     \item{y}{the response matrix used.}
#'     \item{AIC, BIC}{Akaike and Bayesian information criteria.}
#'  }
#'
#' @export
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @references Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, \bold{101}, 321--344.
#'      
#'  Medeiros, R. M. R., and Ferrari, S. L. P. (2025). Multivariate Box-Cox symmetric regression 
#'      models generated by a normal scale mixture copula.
#'
#' @examples
#' # Data set: macronutrients (for description run ?macronutrients)
#' # Consider modeling the animal and plant protein regarding income level
#' 
#' ## Response distribution:
#' mvplot(macronutrients[, 1:2])
#' 
#' ## Animal and plant protein by income level:
#' boxplot(animal_p ~ income, macronutrients)
#' boxplot(plant_p ~ income, macronutrients)
#' 
#' ## Reference model
#' fit0 <- bcsnsm(animal_p + plant_p ~ income | income, data = macronutrients)
#' summary(fit0)
#' 
#' ## Improved model
#' fit <- bcsnsm(animal_p + plant_p ~ income | income, data = macronutrients,
#'               copula = "t", eta = 0.9, margins = c("bcloii", "bct"))
#' summary(fit)               
#' 
#' ## Summary with a specific margin
#' print(summary(fit), margin = 1)
#' print(summary(fit), margin = 2)
summary.bcsnsm <- function(object, ...)
{

  y <- object$y
  n <- object$nobs
  d <- object$d

  links <- object$links
  margins <- object$margins
  association <- object$association
  gamma <- object$gamma
  copula <- object$copula
  eta <- object$eta

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)
  
  ## Summary for quantile residuals
  res <- stats::residuals(object)

  # Summary for beta
  beta_coef <- object$coefficients
  beta_se <- lapply(stats::vcov(object, "beta"), function(x) sqrt(diag(x)))
  beta_z <- Map(function(x, y) x / y, beta_coef, beta_se)
  beta_pvalue <- Map(function(x) 2 * stats::pnorm(abs(x), lower.tail = FALSE), beta_z)
  beta_tab <- list()
  for (j in 1:d){
    beta_tab[[j]] <- cbind(Est = beta_coef[[j]],
                           `Std. Error` = beta_se[[j]],
                           `z value` = beta_z[[j]],
                            `Pr(>|z|)` = beta_pvalue[[j]])
  }
  names(beta_tab) <- colnames(y)
  
  # Summary for sigma
  sigma_coef <- object$marginal.pars$sigma
  sigma_se <- sqrt(diag(as.matrix(stats::vcov(object, "sigma"))))
  sigma_tab <- matrix(c(sigma_coef, sigma_se), byrow = TRUE, nrow = 2)
  colnames(sigma_tab) <- colnames(y)
  rownames(sigma_tab) <- c(" ", "Std. Error")

  # Summary for lambda
  lambda_tab <- NULL
  if (any(lambda_id)){
    lambda_coef <- lambda_se <- rep(NA, d)

    lambda_coef[lambda_id] <- object$marginal.pars$lambda
    lambda_se[lambda_id] <- sqrt(diag(as.matrix(stats::vcov(object, "lambda"))))

    lambda_tab <- matrix(c(lambda_coef, lambda_se), byrow = TRUE, nrow = 2)
    colnames(lambda_tab) <- colnames(y)
    rownames(lambda_tab) <- c(" ", "Std. Error")
  }

  nu_tab <- NULL
  if (any(nu_id)) {
    nu_coef <- nu_se <- rep(NA, d)

    nu_coef[nu_id] <- object$marginal.pars$nu
    nu_se[nu_id] <- sqrt(diag(as.matrix(stats::vcov(object, "nu"))))

    nu_tab <- matrix(c(nu_coef, nu_se), byrow = TRUE, nrow = 2)
    colnames(nu_tab) <- colnames(y)
    rownames(nu_tab) <- c(" ", "Std. Error")
  }

  # Summary for gamma
  if (object$association == "non-associative"){
      Gamma <- diag(d)
  }else{
      Gamma <- get(tolower(object$association), envir = asNamespace("bcsnsm"))(d)$Gamma(gamma)
  }
  colnames(Gamma) <- rownames(Gamma) <- colnames(y)
  

  out <- list(call = object$call,
              beta_tab = beta_tab,
              sigma_tab = sigma_tab,
              lambda_tab = lambda_tab,
              nu_tab = nu_tab,
              Gamma = Gamma,
              margins = margins,
              links = links,
              association = association,
              copula = copula,
              eta = eta,
              residuals = res,
              logLik = object$logLik,
              y = y,
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = log(n)))

  class(out) <- "summary.bcsnsm"
  out
}

# Print summary
#' @rdname summary.bcsnsm
#' @export
print.summary.bcsnsm <- function(x, digits = max(3, getOption("digits") - 3), 
                                 margin, ...)
{
  
  y <- x$y
  n <- nrow(y)
  d <- ncol(y)
  
  if (!missing(margin)) {
    if (margin > d | margin < 1) 
      stop(paste0("The argument margin = ", margin, " is outside the dimensions of the data: 1:", d))
  }
  
  copula <- x$copula
  eta <- x$eta
  association <- x$association
  margins <- x$margins
  links <- x$links

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")
  gamma_name <- paste(toupper(substr(association, 1, 1)), substr(association, 2, nchar(association)), sep = "")

  cat(crayon::cyan("\nMultivariate BCS-NSM Regression\n"))

  cat(paste0(crayon::cyan("\nCall:\n"), 
             paste(deparse(x$call, width.cutoff = floor(getOption("width") * 0.7)), collapse = "\n")))

  cat(crayon::cyan("\n\nQuantile residuals:\n"))
  res <- x$residuals
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  residuals_tab <- cbind(mean(res), stats::sd(res), skewness, kurtosis)
  colnames(residuals_tab) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(residuals_tab) <- " "
  print(round(residuals_tab, digits = digits))

  cat(crayon::cyan("\nMarginal models:\n"))
  if (missing(margin)) {
    for (j in 1:d) {
      cat("\n", crayon::cyan(colnames(y)[j]), crayon::cyan(" ~ "),
          crayon::cyan(get(x$margins[j], envir = asNamespace("bcsnsm"))()$name),
          crayon::cyan(" Distribution:\n"), sep = "")
      
      
      cat("\nRegression coefficients with", links[j], "link function:\n")
      cmat <- round(x$beta_tab[[j]], digits)
      stats::printCoefmat(cmat, digits = digits)
      
      TAB <- rbind(
        round(x$sigma_tab[, j], digits),
        if (lambda_id[j]) round(x$lambda_tab[, j], digits),
        if (nu_id[j]) round(x$nu_tab[, j], digits)
      )
      
      rownames(TAB) <- c("sigma", if (lambda_id[j]) "lambda", if (nu_id[j]) "nu")
      colnames(TAB) <- c("Est", "Std. Error")
      cat("\nFurther marginal parameters:\n")
      print(TAB, digits = digits)
    
    }  
  } else {
    j <- margin
    cat("\n", crayon::cyan(colnames(y)[j]), crayon::cyan(" ~ "),
        crayon::cyan(get(x$margins[j], envir = asNamespace("bcsnsm"))()$name),
        crayon::cyan(" Distribution:\n"), sep = "")
    
    cat("\nRegression coefficients with", links[j], "link function:\n")
    cmat <- round(x$beta_tab[[j]], digits)
    stats::printCoefmat(cmat, digits = digits)
      
    TAB <- rbind(
      round(x$sigma_tab[, j], digits),
      if (lambda_id[j]) round(x$lambda_tab[, j], digits),
      if (nu_id[j]) round(x$nu_tab[, j], digits)
    )
    
    rownames(TAB) <- c("sigma", if (lambda_id[j]) "lambda", if (nu_id[j]) "nu")
    colnames(TAB) <- c("Est", "Std. Error")
    cat("\nFurther marginal parameters:\n")
    print(TAB, digits = digits)
  }

  Gamma <- x$Gamma
  Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
  cat(crayon::cyan("\n", gamma_name, " association matrix:\n", sep = ""))
  print(round(Gamma, digits), na.print = ".")

  df <- n - x$AIC / 2 - x$logLik
  
  cat(
    "\n---",
    crayon::cyan("\nMargins:"), margins,
    crayon::cyan("\nCopula:"), if (is.null(eta)) copula_name else paste0(copula_name, "(", round(eta, 2), ")"),
    crayon::cyan("\nAssociation:"), gamma_name,
    crayon::cyan("\nlogLik:"), round(x$logLik, digits), "on", df, "df",
    crayon::cyan("\nAIC:"), round(x$AIC, digits), "|",
    crayon::cyan("BIC:"), round(x$BIC, digits), "\n"
  )
  
  invisible(x)
}

globalVariables(c("z", "upper", "lower", "density", "theo", "emp", "grid_x", "grid_y", "prob"))

# Plot ----------------------------------------------------------------------------------------
#' Diagnostic Plots for a Multivariate BCS-NSM Regression Fit
#' 
#' This function provides tools for evaluating the goodness-of-fit of a BCS-NSM
#'     regression to a multivariate positive data. 
#'
#' @param x an object of class \code{"bcsnsm"}, a result of a call to \code{\link{bcsnsm}}.
#' @param type character; specifies which graphical should be produced. The available options
#'     are \code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and
#'     \code{"epsilon"}. Se details below.
#' @param levels levels for contours plots. Used only when \code{type = "epsilon"}.
#' @param panel a vector of the form \code{c(nr, nc)} with the number of rows and columns, where
#'     the figures will be drawn in an \code{nr-}by\code{-nc} array on the device. Used
#'     only when \code{type = "marginal"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#' 
#' @details
#' The argument \code{type} specifies what type of diagnostic plot should be produced. The available
#'     options are \code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and
#'     \code{"epsilon"}. When \code{type = "quantile"}, the function returns the normal probability
#'     plot with a 95\% confidence region of the overall quantile residuals, where the confidence
#'     region is computed according to Fox (2015, Sec. 3.1.3). 
#'     
#'     When \code{type = "marginal"}, it 
#'     produces the normal probability plots of the quantile residuals resulting from the fit at
#'     each marginal distribution. The argument \code{panel} is a vector \code{c(nr, nc)},
#'     indicating the number of rows and columns in which these graphs will be arranged on the 
#'     device. 
#'     
#'     When \code{type = "mahalanobis"}, a plot of the transformed Mahalanobis distances 
#'     for each observational unit is produced. Finally, when \code{type = "epsilon"}, a panel
#'     graphical summary of the epsilon transformations is returned, where the bottom panels show 
#'     the scatterplots with the fitted contour lines for
#'     each pair at \code{levels}; the diagonal panels show the histograms with the fitted marginal
#'     densities; and the top panels show the upper triangular elements of the fitted association matrix. 
#'     
#'     Details about the diagnostic tools used by the \code{plot.bcsnsm} function can be found in 
#'     \code{\link{residuals.bcsnsm}} and in Medeiros and Ferrari (2025).
#' 
#' @return \code{plot} method for \code{"\link{bcsnsm}"} objects returns four types
#'     of diagnostic plots.
#'
#' @references Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, \bold{101}, 321--344.
#'      
#'  Fox, J. (2015). \emph{Applied Regression Analysis and Generalized Linear Models}. 
#'      Sage Publications.
#'      
#'  Medeiros, R. M. R., and Ferrari, S. L. P. (2025). Multivariate Box-Cox symmetric regression 
#'      models generated by a normal scale mixture copula.
#'
#' @export
#' 
#' @examples
#' # Data set: macronutrients (for description run ?macronutrients)
#' # Consider modeling the animal and plant protein regarding income level
#' 
#' ## Response distribution:
#' mvplot(macronutrients[, 1:2])
#' 
#' ## Animal and plant protein by income level:
#' boxplot(animal_p ~ income, macronutrients)
#' boxplot(plant_p ~ income, macronutrients)
#' 
#' ## Reference model
#' fit <- bcsnsm(animal_p + plant_p ~ income | income, data = macronutrients)
#' 
#' ## Overall quantile residuals (default)
#' plot(fit)
#' 
#' ## Marginal quantile residuals
#' plot(fit, type = "marginal")
#' plot(fit, type = "marginal", panel = c(2, 1)) # Change the layout
#' 
#' ## Mahalanobis distances
#' plot(fit, type = "mahalanobis")
#' 
#' ## The epsilon's transformations
#' plot(fit, type = "epsilon")
plot.bcsnsm <- function(x, type = c("quantile", "marginal", "mahalanobis", "epsilon"),
                          levels = c(1e-1, 1e-2, 1e-3, 1e-4), panel = NULL, ...) {

  type <- match.arg(type, c("quantile", "marginal", "mahalanobis", "epsilon"))

  y <- x$y
  n <- x$nobs
  d <- x$d

  copula <- x$copula
  eta <- x$eta
  association <- if (x$association == "non-associative") "nonassociative" else x$association
  Gamma <- get(association, envir = asNamespace("bcsnsm"))(d)$Gamma(x$gamma)
  mcopula <- make_copula(copula, eta / (1 - eta))

  dPSI <- mcopula$dPSI
  dmv <- mcopula$dmv

  # Information for the normal probability plot with a confidence region
  Pi <- (1:n - 0.5) / n
  zi <- stats::qnorm(Pi)

  if (is.null(panel)) panel <- c(ceiling(d / 4), 4)

  ### Quantile -------------------------------------------------------------------------------------
  if (type == "quantile") {

    rq <- sort(stats::residuals(x))
    mu_hat <- stats::median(rq)
    sigma_hat <- stats::IQR(rq) / 1.349

    se <- sigma_hat * sqrt(Pi * (1 - Pi) / n) / stats::dnorm(zi)
    rq_hat <- mu_hat + sigma_hat * zi

    positions <- data.frame(x = c(zi, rev(zi)),
                            y = c(rq_hat - stats::qnorm(1 - 0.05 / 2) * se,
                                  rev(rq_hat + stats::qnorm(1 - 0.05 / 2) * se)))

    ggplot2::ggplot(data.frame(zi = zi, rq = rq)) +
      ggplot2::geom_polygon(data = positions, ggplot2::aes(x = x, y = y), fill = "#cceff1") +
      ggplot2::geom_point(ggplot2::aes(x = zi, y = rq), pch = "+", cex = 4) +
      ggplot2::xlab("Normal quantiles") + ggplot2::ylab("Residual quantiles") +
      ggplot2::geom_qq_line(ggplot2::aes(sample = rq), col = "#56B1F7", lwd = 1)

    ### Marginal -----------------------------------------------------------------------------------
  } else if (type == "marginal") {

    rm <- stats::residuals(x, "marginal")

    res <- data.frame(emp = c(apply(rm, 2, sort)),
                      theo = rep(stats::qnorm(stats::ppoints(n)), d),
                      marg = factor(rep(colnames(y), each = n), levels = colnames(y)))

    lower_f <- function(res, z) {
      robust_sd <- stats::IQR(res) / 1.349
      robust_line <- stats::median(res) + z * robust_sd
      robust_se <-  robust_sd * sqrt(stats::pnorm(z) * (1 - stats::pnorm(z)) / n) / stats::dnorm(z)
      robust_line - stats::qnorm(1 - 0.05/2) * robust_se
    }

    upper_f <- function(res, z) {
      robust_sd <- stats::IQR(res) / 1.349
      robust_line <- stats::median(res) + z * robust_sd
      robust_se <-  robust_sd * sqrt(stats::pnorm(z) * (1 - stats::pnorm(z)) / n) / stats::dnorm(z)
      robust_line + stats::qnorm(1 - 0.05/2) * robust_se
    }

    band <- data.frame(
      lower = c(apply(rm, 2, lower_f, z = seq(min(res$theo), max(res$theo), length.out = 500))),
      upper = c(apply(rm, 2, upper_f, z = seq(min(res$theo), max(res$theo), length.out = 500))),
      marg = factor(rep(colnames(y), each = 500), levels = colnames(y)),
      z = rep(seq(min(res$theo), max(res$theo), length.out = 500), d)
    )

    ggplot2::ggplot(res) +
      ggplot2::geom_ribbon(ggplot2::aes(x = z, ymax = upper, ymin = lower),  data = band, fill = "#cceff1",  color = "#cceff1",
                           show.legend = FALSE) +
      ggplot2::geom_point(ggplot2::aes(x = theo, y = emp), pch = "+", cex = 2) +
      ggplot2::geom_qq_line(ggplot2::aes(sample = emp), col = "#56B1F7", lwd = 1) +
      ggplot2::facet_wrap(~ marg, nrow = panel[1], ncol = panel[2]) +
      ggplot2::labs(x = "Theoretical quantiles",
                    y = "Sample quantiles")

    ### Mahalanobis ----------------------------------------------------------------------------------
  } else if (type == "mahalanobis") {

    epsilon <- stats::residuals(x, type = "epsilon")
    mahalanobis <- mahalanobis(epsilon, rep(0L, d), Gamma)

    ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x = 1:n, y = rep(0, n),
                                         xend = 1:n, yend = mahalanobis)) +
      ggplot2::labs(x = "Index observations", y = "Transformed Mahalanobis distances")

    ### Epsilon --------------------------------------------------------------------------------------
  } else if (type == "epsilon") {

    epsilon <- stats::residuals(x, "epsilon")

    diag_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)

      ggplot2::ggplot(data, mapping) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                colour = 1, fill = "white",
                                bins = ceiling(1 + 3.33 * log(n))) +
        ggplot2::geom_function(fun = dPSI, col = "#00AFBB")


    }

    ### Upper plots
    upper_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)
      y <- GGally::eval_data_col(data, mapping$y)

      i <- which(colSums(epsilon - x) == 0)
      j <- which(colSums(epsilon - y) == 0)

      corr <- Gamma[i, j]

      colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
      fill <- colFn(1000)[findInterval(corr, seq(-1, 1, length = 1000))]

      GGally::ggally_text(
        label = as.character(round(corr, 4)),
        mapping = ggplot2::aes(),
        xP = 0.5,
        yP = 0.5,
        cex = 2.5,
        color = 'black', ...) +
        ggplot2::theme_void() +
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = fill))
    }

    ### Lower plots
    lower_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)
      y <- GGally::eval_data_col(data, mapping$y)

      i <- which(colSums(epsilon - x) == 0)
      j <- which(colSums(epsilon - y) == 0)

      Gamma_aux <- matrix(c(1, Gamma[i, j], Gamma[i, j], 1), 2, 2)


      grid <- expand.grid(seq(min(x) - 10, max(x) + 10, length.out = 200),
                          seq(min(y) - 10, max(y) + 10, length.out = 200))

      data_aux <- data.frame(grid_x = grid[, 1],
                             grid_y = grid[, 2],
                             prob = dmv(cbind(grid[, 1], grid[, 2]),
                                        Gamma = Gamma_aux))


      ggplot2::ggplot() +
        ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_contour(data = data_aux, mapping =  ggplot2::aes(x = grid_x, y = grid_y, z = prob),
                              col = "#00AFBB", breaks = levels) +
        ggplot2::labs(x = "", y = "")


    }

    colFn <- grDevices::colorRampPalette(c("brown1", "white", "dodgerblue"), interpolate ='spline')
    lab <- ggplot2::ggplot(data.frame(x = stats::runif(200, -1, 1), y = stats::runif(200, -1, 1),
                                      z = stats::runif(200, -1, 1)), ggplot2::aes(x, y, colour = z)) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient2("Fitted association \nparameter",
                                      low = colFn(200)[1], high = colFn(200)[200]) +
      ggplot2::theme(legend.title.align = 0.5,
                     legend.position = "top",
                     legend.key.height = ggplot2::unit(0.3, 'cm'),
                     legend.key.width = ggplot2::unit(1.5, 'cm'))

    GGally::ggpairs(as.data.frame(epsilon),
                    upper = list(continuous = GGally::wrap(upper_func)),
                    lower = list(continuous = GGally::wrap(lower_func)),
                    diag = list(continuous = GGally::wrap(diag_func)),
                    legend = GGally::grab_legend(lab),
                    progress = FALSE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::theme(legend.position = "top")

  }


}
