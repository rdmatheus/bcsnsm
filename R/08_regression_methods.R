#' @name bcsnsmreg-methods
#' @title Methods for \code{"bcsnsmreg"} objects
#' @param x,object an object of class \code{"bcsnsmreg"}.
#' @param digits a non-null value for \code{digits} specifies the minimum number of
#'     digits to be printed in values.
#' @param margin an integer which specifies the marginal distribution.
#' @param formula a symbolic description of the model to be fitted to each marginal distribution.
#'     For instance, \code{formula = c(y1 + y2 + y3 ~ x1 + x2 | x1 + x3 | x2 + x3)} fits a
#'     3-variate BCS-NSM regression that models \code{y1} as a function of \code{x1} and \code{x2};
#'     \code{y2} as a function of \code{x1} and \code{x3}; and \code{y3} as a function of \code{x2}
#'     and \code{x3}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical Akaike information criteria (AIC).
#' @param ... additional argument(s) for methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL

## Model frame
#' @export
#' @rdname bcsnsmreg-methods
model.frame.bcsnsmreg <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname bcsnsmreg-methods
model.matrix.bcsnsmreg <- function(object, margin, ...) {
  f <- object$formula
  mf <- stats::model.frame(object)
  stats::model.matrix(f, data = mf, rhs = margin)
}

#  Variance-covariance matrix
#' @rdname bcsnsmreg-methods
#' @param parm a character; specifies which parameters are to be given the asymptotic covariance
#'     matrix of the maximum likelihood estimators (for \code{vcov} function) or the confidence
#'     intervals (for \code{confint} function). The options are \code{"all"} (default), where
#'     all parameters are considered, \code{"beta"}, \code{"sigma"}, \code{"lambda"}, \code{"nu"},
#'     and \code{"gamma"}.
#'
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
vcov.bcsnsmreg <- function(object, parm = c("all", "beta", "sigma", "lambda", "nu", "gamma"), ...) {
  
  parm <- match.arg(parm, c("all", "beta", "sigma", "lambda", "nu", "gamma"))
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
         "all" = covm,
         "beta" = vcov_beta,
         "sigma" = vcov_sigma,
         "lambda" = vcov_lambda,
         "nu" = vcov_nu,
         "gamma" = vcov_gamma
  )
}

# Confidence intervals
#' @rdname bcsnsmreg-methods
#' @param level the confidence level required.
#' @export
confint.bcsnsmreg <- function(object, parm = c("full", "beta", "sigma", "lambda", "nu", "gamma"), level = 0.95, ...) {

  par <- match.arg(par, c("full", "beta", "sigma", "lambda", "nu", "gamma"))
  est <- object$optim_pars$par
  variance <- diag(object$vcov)

  y <- object$y
  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)
  par_id <- object$optim_pars$par_id

  ## Beta
  beta_id <- object$optim_pars$par_id$beta
  var_beta <- Map(function(id) variance[id], beta_id)
  var_beta <- stats::setNames(variance, colnames(y))

  ## Sigma
  est_sigma <- est[par_id$sigma]
  var_sigma <- variance[par_id$sigma]
  names(est_sigma) <- names(var_sigma) <- colnames(y)

  ## Lambda
  est_lambda <- est[par_id$lambda]
  var_lambda <- variance[par_id$lambda]
  names(est_lambda) <- names(var_lambda) <- colnames(y)[lambda_id]

  ## Nu
  est_nu <- est[par_id$nu]
  var_nu <- variance[par_id$nu]
  names(est_nu) <- names(var_nu) <- colnames(y)[nu_id]

  ## Gamma
  est_gamma <- est[par_id$gamma]
  var_gamma <- variance[par_id$gamma]

  if (object$association == "uniform") {
    names_gamma <- "gamma"
  } else if (object$association == "unstructured") {
    names_gamma <- vector()
    id <- matrix(which(upper.tri(diag(object$d)), arr.ind = TRUE)[order(which(upper.tri(diag(object$d)), arr.ind = TRUE)[, 1]), ], ncol = 2)
    for (i in 1:nrow(id)) {
      names_gamma[i] <- paste0("gamma", id[i, 1], id[i, 2])
    }
  }

  names(est_gamma) <- names(var_gamma) <- names_gamma

  qtl <- stats::qnorm(1 - (1 - level) / 2)

  switch(par,
         "full" = {

           # TO FIX
           # out <- matrix(c(est - qtl * sqrt(variance),
           #                 est + qtl * sqrt(variance)), ncol = 2)
           #
           # colnames(out) <- paste0(c(100 * (1 - level)/2,
           #                           100 * (1 - (1 - level)/2)), "%")

         },

         "beta" = {
           # TO FIX
           # out <- matrix(c(est_beta - qtl * sqrt(var_beta),
           #                 est_beta + qtl * sqrt(var_beta)), ncol = 2)
           #
           # rownames(out) <- names(est_beta)
           # colnames(out) <- paste0(c(100 * (1 - level)/2,
           #                           100 * (1 - (1 - level)/2)), "%")
         },

         "sigma" = {
           out <- matrix(c(est_sigma - qtl * sqrt(var_sigma),
                           est_sigma + qtl * sqrt(var_sigma)), ncol = 2)

           rownames(out) <- names(est_sigma)
           colnames(out) <- paste0(c(100 * (1 - level)/2,
                                     100 * (1 - (1 - level)/2)), "%")
         },

         "lambda" = {
           out <- matrix(c(est_lambda - qtl * sqrt(var_lambda),
                           est_lambda + qtl * sqrt(var_lambda)), ncol = 2)

           rownames(out) <- names(est_lambda)
           colnames(out) <- paste0(c(100 * (1 - level)/2,
                                     100 * (1 - (1 - level)/2)), "%")
         },

         "nu" = {
           out <- matrix(c(est_nu - qtl * sqrt(var_nu),
                           est_nu + qtl * sqrt(var_nu)), ncol = 2)

           rownames(out) <- names(est_nu)
           colnames(out) <- paste0(c(100 * (1 - level)/2,
                                     100 * (1 - (1 - level)/2)), "%")
         },

         "gamma" = {
           out <- matrix(c(est_gamma - qtl * sqrt(var_gamma),
                           est_gamma + qtl * sqrt(var_gamma)), ncol = 2)

           rownames(out) <- names(est_gamma)
           colnames(out) <- paste0(c(100 * (1 - level)/2,
                                     100 * (1 - (1 - level)/2)), "%")
         }
  )

  out

}

# Log-likelihood
#' @rdname bcsnsmreg-methods
#' @export
logLik.bcsnsmreg <- function(object, ...) {
  structure(object$logLik,
            df = length(object$optim_pars$par) + as.numeric(!is.null(object$eta)),
            class = "logLik")
}

# AIC
#' @export
#' @rdname bcsnsmreg-methods
AIC.bcsnsmreg <- function(object, ..., k = 2) {
  npars <- length(object$optim_pars$par) + as.numeric(!is.null(object$eta))
  AIC <- -2 * object$logLik + k * npars

  class(AIC) <- "AIC"
  return(AIC)
}

# Residuals
#' Extract Residuals for an \code{bcsnsmreg} Object
#' 
#' It provides the overall and marginal quantile residuals, transformed Mahalanobis distances,
#'      and the epsilon transformations, \eqn{\varepsilon_{ij} = \Psi_H^{-1}\{\widehat{F}_{j}(y_{ij})\}}.
#'
#' @param object an \code{"bcsnsmreg"} object.
#' @param type character; specifies which residual should be extracted. The available arguments are
#'     \code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and \code{"epsilon"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
residuals.bcsnsmreg <- function(object, type = c("quantile", "marginal", "mahalanobis", "epsilon"), ...) {

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

# Print
#' @rdname bcsnsmreg-methods
#' @export
print.bcsnsmreg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

  y <- x$y
  n <- x$nobs
  d <- x$d

  margins <- x$margins
  association <- x$association
  gamma <- x$gamma
  copula <- x$copula
  eta <- x$eta

  par_id <- x$optim_pars$par_id
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  lambda <- rep(NA, d)
  names(lambda) <- colnames(y)
  lambda[lambda_id] <- x$marginal.pars$lambda

  nu <- rep(NA, d)
  names(nu) <- colnames(y)
  nu[nu_id] <- x$marginal.pars$nu

  # Names
  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")
  gamma_name <- paste(toupper(substr(association, 1, 1)), substr(association, 2, nchar(association)), sep = "")

  links <- x$links
  beta <- x$coefficients

  cat(crayon::cyan(
    "\nMultivariate BCS-NSM Regression with",
    if (is.null(eta)) copula_name else paste0(copula_name, "(", round(eta, 2), ")"),
    "Copula\n\n"
  ))

  cat(crayon::cyan("\nCall:\n"), deparse(x$call, width.cutoff = floor(getOption("width") * 0.7)), "", sep = "\n")

  if (x$optim_pars$convergence > 0) {
    cat("\nmodel did not converge\n")
  } else {

    cat(crayon::cyan("\nRegression coefficients:\n\n"))

    X <- lapply(1:d, function(j) stats::model.matrix(x$formula, data = stats::model.frame(x), rhs = j))
    X_names <- unique(unlist(lapply(X, function(x) colnames(x))))

    BETA <- matrix(NA, length(X_names), d)
    colnames(BETA) <- colnames(y)
    rownames(BETA) <- X_names

    for (j in 1:d) {
      BETA[names(beta[[j]]), j] <- round(beta[[j]], digits)
    }

    print(BETA, na.print = "-")


    tab <- cbind(sigma = round(x$marginal.pars$sigma, digits),
                 if(any(lambda_id)) round(lambda, digits) else NULL,
                 if(any(nu_id)) round(nu, digits) else NULL)

    colnames(tab) <- c("sigma", if(any(lambda_id)) "lambda", if(any(nu_id)) "nu")

    cat(crayon::cyan("\nFurther marginal parameters:\n\n"))
    print(tab, na.print = "-")

    if (length(gamma) > 0) {

      cat(crayon::cyan("\n", gamma_name, " association matrix:\n\n", sep = ""))

      if (x$association == "non-associative"){
        Gamma <- diag(d)
      }else{
        Gamma <- get(tolower(x$association), envir = asNamespace("bcsnsm"))(d)$Gamma(round(gamma, digits))
      }

      Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
      colnames(Gamma) <- rownames(Gamma) <- colnames(y)
      print(Gamma, na.print = ".")
    }

    cat(
      "\n---",
      crayon::cyan("\nMargins:"), margins,
      crayon::cyan("\nlogLik:"), round(x$logLik, digits), "|",
      crayon::cyan("AIC:"), round(stats::AIC(x), digits), "|",
      crayon::cyan("BIC:"), round(stats::AIC(x, k = log(n)), digits), "\n"
    )
  }

  invisible(x)
}

# Summary
#' @rdname bcsnsmreg-methods
#' @export
summary.bcsnsmreg <- function(object, ...)
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
  gamma_id <- length(gamma) > 0

  ## Summary for quantile residuals
  res <- stats::residuals(object)
  skewness <- mean((res - mean(res))^3) / (stats::sd(res)^3)
  kurtosis <- mean((res - mean(res))^4) / (stats::sd(res)^4)
  residuals_tab <- cbind(mean(res), stats::sd(res), skewness, kurtosis)
  colnames(residuals_tab) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(residuals_tab) <- " "

  # Summary for mu
  beta_coef <- object$coefficients
  beta_se <- lapply(stats::vcov(object, "beta"), function(x) sqrt(diag(x)))
  beta_z <- Map(function(x, y) x / y, beta_coef, beta_se)
  beta_pvalue <- Map(function(x) 2 * stats::pnorm(abs(x), lower.tail = FALSE), beta_z)

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
  gamma_tab <- NULL
  if (gamma_id) {

    if (object$association == "non-associative"){
      Gamma <- diag(d)
    }else{
      Gamma <- get(tolower(object$association), envir = asNamespace("bcsnsm"))(d)$Gamma(gamma)
    }

    Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
    colnames(Gamma) <- rownames(Gamma) <- colnames(y)

  }

  out <- list(call = object$call,
              residuals = residuals_tab,
              links = links,
              beta_coef = beta_coef,
              beta_se = beta_se,
              beta_z = beta_z,
              beta_pvalue = beta_pvalue,
              sigma_tab = sigma_tab,
              lambda_tab = lambda_tab,
              nu_tab = nu_tab,
              Gamma = Gamma,
              margins = margins,
              association = association,
              copula = copula,
              eta = eta,
              logLik = object$logLik,
              n = n,
              d = d,
              y = y,
              gamma_id = gamma_id,
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = log(n)))

  class(out) <- "summary.bcsnsmreg"
  out
}

# Print summary
#' @rdname bcsnsmreg-methods
#' @export
print.summary.bcsnsmreg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

  y <- x$y
  n <- x$n
  d <- x$d
  copula <- x$copula
  eta <- x$eta
  association <- x$association
  margins <- x$margins
  links <- x$links

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")
  gamma_name <- paste(toupper(substr(association, 1, 1)), substr(association, 2, nchar(association)), sep = "")

  cat(crayon::cyan(
    "\nMultivariate BCS-NSM Regression with",
    if (is.null(eta)) copula_name else paste0(copula_name, "(", round(eta, 2), ")"),
    "Copula\n\n"
  ))

  cat(crayon::cyan("\nCall:\n"), deparse(x$call, width.cutoff = floor(getOption("width") * 0.7)), "", sep = "\n")

  cat(crayon::cyan("\nQuantile residuals:\n\n"))
  print(x$residuals, digits = digits)

  for (j in 1:d) {

    cat("\n", crayon::cyan(colnames(y)[j]), crayon::cyan(" ~ "),
        crayon::cyan(get(x$margins[j], envir = asNamespace("bcsnsm"))()$name),
        crayon::cyan(" Distribution:\n"), sep = "")


    cat("\nRegression coefficients with", links[j], "link function:\n")

    cmat <- round(cbind(Est = x$beta_coef[[j]],
                        `Std. Error` = x$beta_se[[j]],
                        `z value` = x$beta_z[[j]],
                        `Pr(>|z|)` = x$beta_pvalue[[j]]), digits)
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

  if (x$gamma_id) {

    cat(crayon::cyan("\n", gamma_name, " association matrix:\n\n", sep = ""))
    print(round(x$Gamma, digits), na.print = ".")

  }

  cat(
    "\n---",
    crayon::cyan("\nlogLik:"), round(x$logLik, digits), "|",
    crayon::cyan("AIC:"), round(x$AIC, digits), "|",
    crayon::cyan("BIC:"), round(x$BIC, digits), "\n"
  )

  invisible(x)
}



globalVariables(c("z", "upper", "lower", "density", "theo", "emp", "grid_x", "grid_y", "prob"))
# Plot
#' Plot Diagnostic for an \code{bcsnsmreg} Object
#' 
#' The graphs provided are tools for evaluating the goodness-of-fit of a BCS-NSM
#'     regression to a multivariate positive data. 
#'
#' @param x an object of class \code{bcsnsmreg}.
#' @param type character; specifies which graphical should be produced. The available options
#'     are \code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and
#'     \code{"epsilon"}. Se Details.
#' @param levels levels for contours plots. Used only when \code{type = "epsilon"}.
#' @param panel a vector of the form \code{c(nr, nc)} with the number of rows and columns, where
#'     the figures will be drawn in an \code{nr-}by\code{-nc} array on the device. Used
#'     only when \code{type = "marginal"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#' 
#' @details
#' the argument \code{type} specifies what type of diagnostic plot should be produced. The available
#'     options are \code{"quantile"} (default), \code{"marginal"}, \code{"mahalanobis"}, and
#'     \code{"epsilon"}. When \code{type = "quantile"}, the function returns the normal probability
#'     plot with a 95\% confidence region of the overall quantile residuals, where the confidence
#'     region is computed according to Fox (2015, Sec. 3.1.3). When \code{type = "marginal"}, it 
#'     produces the normal probability plots of the quantile residuals resulting from the fit at
#'     each marginal distribution. The argument \code{panel} is a vector \code{c(nr, nc)},
#'     indicating the number of rows and columns in which these graphs will be arranged on the 
#'     device. When \code{type = "mahalanobis"}, a plot of the transformed Mahalanobis distances 
#'     for each observational unit is produced. Finally, when \code{type = "epsilon"}, a panel
#'     graphical summary of the transformations \eqn{\varepsilon_{ij} = \Psi_H^{-1}\{\widehat{F}_{j}(y_{ij})\}}
#'     is returned, where the bottom panels show the scatterplots with the fitted contour lines for
#'     each pair at \code{levels}; the diagonal panels show the histograms with the fitted marginal
#'     densities; and the top panels show the upper triangular elements of the fitted association matrix. 
#' 
#'
#' @export
#'
plot.bcsnsmreg <- function(x, type = c("quantile", "marginal", "mahalanobis", "epsilon"),
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
