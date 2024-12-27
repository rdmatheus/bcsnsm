#' @name bcsnsm-methods
#' @title Methods for \code{"bcsnsm"} objects
#' @param x,object an object of class \code{bcsnsm}.
#' @param digits number of digits in print methods.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param ... additional argument(s) for methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL

#  Variance-covariance matrix
#' @rdname bcsnsm-methods
#' @param parm character; specifies which submatrix of the asymptotic covariance matrix of the
#'     maximum likelihood estimators should be returned for the \code{vcov} function or
#'     which parameters are to be given confidence intervals for the \code{confint} function. The
#'     options are \code{"all"} (default), \code{"mu"}, \code{"sigma"}, \code{"lambda"}, \code{"nu"}, and \code{"gamma"}.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
vcov.bcsnsm <- function(object, parm = c("all", "mu", "sigma", "lambda", "nu", "gamma"), ...) {

  parm <- match.arg(parm, c("all", "mu", "sigma", "lambda", "nu", "gamma"))
  covm <- object$vcov

  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  d <- object$d

  names_mu <- paste0("mu", 1:d)
  names_sigma <- paste0("sigma", 1:d)
  names_lambda <- if (any(lambda_id)) paste0("lambda", (1:d)[lambda_id]) else NULL
  names_nu <- if (any(nu_id)) paste0("nu", (1:d)[nu_id]) else NULL

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
  
  colnames(covm) <- rownames(covm) <- c(names_mu, names_sigma, names_lambda, names_nu, names_gamma)

  par_id <- object$optim_pars$par_id

  switch(parm,
         "all" = covm,
         "mu" = covm[par_id$mu, par_id$mu],
         "sigma" = covm[par_id$sigma, par_id$sigma],
         "lambda" = covm[par_id$lambda, par_id$lambda],
         "nu" = covm[par_id$nu, par_id$nu],
         "gamma" = covm[par_id$gamma, par_id$gamma]
  )
}

# Confidence intervals
#' @rdname bcsnsm-methods
#' @param level the confidence level required. 
#' @export
confint.bcsnsm <- function(object, parm = c("all", "mu", "sigma", "lambda", "nu", "gamma"),
                           level = 0.95, digits = max(3, getOption("digits") - 3), ...){

  parm <- match.arg(parm, c("all", "mu", "sigma", "lambda", "nu", "gamma"))

  y <- object$y
  d <- object$d

  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)
  par_id <- object$optim_pars$par_id

  est <- object$optim_pars$par
  se <- sqrt(diag(object$vcov))

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

  names(est) <- names(se) <- c(paste0(paste0("mu", 1:d), " (", colnames(y), ")"),
                               paste0(paste0("sigma", 1:d), " (", colnames(y), ")"),
                               paste0(paste0("lambda", 1:d), " (", colnames(y), ")")[lambda_id],
                               paste0(paste0("nu", 1:d), " (", colnames(y), ")")[nu_id],
                               names_gamma)

  qtl <- stats::qnorm(1 - (1 - level) / 2)

  out <- round(matrix(c(est - qtl * se, est + qtl * se), ncol = 2), digits)

  colnames(out) <- paste0(c(100 * (1 - level)/2, 100 * (1 - (1 - level)/2)), "%")
  rownames(out) <- names(est)

  switch(parm,
         "all" = out,

         "mu" = out[par_id$mu, ],

         "sigma" = out[par_id$sigma, ],

         "lambda" = out[par_id$lambda, ],

         "nu" = out[par_id$nu, ],

         "gamma" = out[par_id$gamma, ]
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
  AIC

}

# Residuals
#' @name residuals.bcsnsm
#' @title Extract Model Residuals for a BerG Regression
#'
#' @param object an \code{"bcsnsm"} object.
#' @param type character; specifies which residual should be extracted. The available arguments are
#'     \code{"mahalanobis"} (default), and \code{"epsilon"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
residuals.bcsnsm <- function(object, type = c("mahalanobis", "epsilon"), ...) {

  type <- match.arg(type, c("mahalanobis", "epsilon"))

  y <- object$y

  n <- object$nobs
  d <- object$d
  epsilon <- matrix(NA, n, d)

  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  copula <- object$copula
  eta <- object$eta

  mcopula <- make_copula(copula, eta / (1 - eta))
  qPSI <- mcopula$qPSI

  mu <- object$mu
  sigma <- object$sigma
  lambda <- rep(NA, d)
  lambda[lambda_id] <- object$lambda
  nu <- rep(NA, d)
  nu[nu_id] <- object$nu


  for (j in 1:d) {

    epsilon[, j] <- qPSI(get(paste0("p", margins[j]), envir = asNamespace("bcsnsm"))(q = y[, j],
                                                      mu = mu[j],
                                                      sigma = sigma[j],
                                                      lambda = lambda[j],
                                                      nu = nu[j]))

  }

  epsilon
  colnames(epsilon) <- colnames(y)
  
  association <- if (object$association == "non-associative") "nonassociative" else object$association
  Gamma <- get(association, envir = asNamespace("bcsnsm"))(d)$Gamma(object$gamma)

  # Squared Mahalanobis distance
  mahalanobis <- Rfast::mahala(epsilon, rep(0L, d), Gamma)

  # Out
  res <- switch(type,
                "mahalanobis" = as.numeric(mahalanobis),
                "epsilon" = epsilon
  )

  res
}

# Print
#' @rdname bcsnsm-methods
#' @export
print.bcsnsm <- function(x, digits = max(3, getOption("digits") - 3), ...){

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
  lambda[lambda_id] <- x$lambda

  nu <- rep(NA, d)
  names(nu) <- colnames(y)
  nu[nu_id] <- x$nu

  # Names
  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")
  gamma_name <- paste(toupper(substr(association, 1, 1)), substr(association, 2, nchar(association)), sep = "")

  cat(crayon::cyan(
    "\nMultivariate BCS-NSM Distribution with",
    if (is.null(eta)) copula_name else paste0(copula_name, "(", round(eta, 2), ")"),
    "Copula\n\n"
  ))


  cat(crayon::cyan("\nCall:\n"), deparse(x$call, width.cutoff = floor(getOption("width") * 0.7)), "", sep = "\n")

  if (x$optim_pars$convergence > 0) {

    cat("\nmodel did not converge\n")

  } else {

    cat(crayon::cyan("\nMarginal fit:\n\n"))

    # Parameters
    mu <- x$mu
    sigma <- x$sigma

    comma_lambda <- comma_nu <- rep(",", d)
    if (any(lambda_id)) lambda <- round(lambda, digits)
    if (any(nu_id)) nu <- round(nu, digits)

    lambda[!lambda_id] <- comma_lambda[!lambda_id] <- ""
    nu[!nu_id] <- comma_nu[!nu_id] <- ""

    for (j in 1:d) {

      cat(colnames(y)[j], " ~ ", toupper(margins)[j], "(", round(mu[j], digits), ",",
          round(sigma[j], digits), comma_lambda[j], lambda[j], comma_nu[j], nu[j], ")\n",
          sep = ""
      )

    }

    if (length(gamma) > 0) {

      cat(crayon::cyan("\n", sub("(.)", "\\U\\1", x$association, perl = TRUE),
                       " association matrix:\n\n", sep = ""))

      Gamma <- get(x$association, envir = asNamespace("bcsnsm"))(d)$Gamma(round(x$gamma, digits))
      Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
      colnames(Gamma) <- rownames(Gamma) <- colnames(y)
      print(Gamma, na.print = ".", digits = digits)

    }

    cat(
      "\n---",
      crayon::cyan("\nlogLik:"), round(x$logLik, digits), "|",
      crayon::cyan("AIC:"), round(stats::AIC(x), digits), "|",
      crayon::cyan("BIC:"), round(stats::AIC(x, k = log(n)), digits), "\n"
    )
  }

  invisible(x)

}

# Summary
#' @rdname bcsnsm-methods
#' @export
summary.bcsnsm <- function(object, ...){

  y <- object$y

  n <- object$nobs
  d <- object$d

  margins <- object$margins
  association <- object$association
  gamma <- object$gamma
  copula <- object$copula
  eta <- object$eta

  # Parameters
  mu <- object$mu
  sigma <- object$sigma

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  lambda <- rep(NA, d)
  names(lambda) <- colnames(y)
  lambda[lambda_id] <- object$lambda

  nu <- rep(NA, d)
  names(nu) <- colnames(y)
  nu[nu_id] <- object$nu

  # Names
  y_names <- if (is.null(colnames(y))) paste0("y", 1:d) else colnames(y)
  Margins <- toupper(margins)
  copula_name <- paste(toupper(substr(copula, 1, 1)), substr(copula, 2, nchar(copula)), sep = "")

  # Summary for mu
  se_mu <- sqrt(diag(stats::vcov(object, "mu")))
  TAB_mu <- cbind(`Estimate` = c(mu), `Std. error` = se_mu)
  rownames(TAB_mu) <- colnames(y)

  # Summary for sigma
  se_sigma <- sqrt(diag(stats::vcov(object, "sigma")))
  TAB_sigma <- cbind(`Estimate` = c(sigma), `Std. error` = se_mu)
  rownames(TAB_sigma) <- colnames(y)

  # Summary for lambda
  TAB_lambda <- NULL
  if (any(lambda_id)) {
    se_lambda <- sqrt(diag(stats::vcov(object, "lambda")))
    TAB_lambda <- cbind(`Estimate` = c(lambda[lambda_id]), `Std. error` = se_lambda)
    rownames(TAB_lambda) <- colnames(y)[lambda_id]
  }

  # Summary for nu
  TAB_nu <- NULL
  if (any(nu_id)) {
    se_nu <- sqrt(diag(stats::vcov(object, "nu")))
    TAB_nu <- cbind(`Estimate` = nu[nu_id], `Std. error` = se_nu)
    rownames(TAB_nu) <- colnames(y)[nu_id]
  }

  # Summary for gamma
  if (length(gamma) > 0) {

    if (object$association == "non-associative"){
      Gamma <- diag(d)
    }else{
      Gamma <- get(tolower(object$association), envir = asNamespace("bcsnsm"))(d)$Gamma(gamma)
    }

    Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
    colnames(Gamma) <- rownames(Gamma) <- colnames(y)

  } else {
    Gamma <- NULL
  }

  npar <- length(object$optim_pars$par)
  out <- list(
    mu = TAB_mu,
    sigma = TAB_sigma,
    lambda = TAB_lambda,
    nu = TAB_nu,
    Gamma = Gamma,
    margins = margins,
    association = association,
    copula = copula,
    y = y,
    d = d,
    eta = eta,
    logLik = object$logLik,
    AIC = stats::AIC(object),
    BIC = stats::AIC(object, k = log(n)),
    call = object$call
  )

  class(out) <- "summary.bcsnsm"
  out
}

# Print summary
#' @rdname bcsnsm-methods
#' @export
print.summary.bcsnsm <- function(x, digits = max(3, getOption("digits") - 3), ...){

  d <- x$d
  eta <- x$eta
  margins <- x$margins

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  # Names
  y_names <- if (is.null(colnames(x$y))) paste0("y", 1:d) else colnames(x$y)
  Margins <- toupper(x$margins)
  copula_name <- paste(toupper(substr(x$copula, 1, 1)), substr(x$copula, 2, nchar(x$copula)), sep = "")

  # Print
  cat(crayon::cyan(
    "\nMultivariate Box-Cox Distribution with",
    if (is.null(eta)) copula_name else paste0(copula_name, "(", round(eta, 2), ")"),
    "Copula\n\n"
  ))

  cat(crayon::cyan("Call:\n\n"))
  print(x$call)

  cat(crayon::cyan("\nMarginal fit:\n\n"))

  nu <- lambda <- matrix(NA, d, 2)
  colnames(nu) <- colnames(lambda) <- colnames(x$lambda)
  lambda[lambda_id, ] <- x$lambda
  nu[nu_id, ] <- x$nu

  for (j in 1:d) {

    cat(y_names[j], " ~ ", get(x$margins[j], envir = asNamespace("bcsnsm"))()$name, " Distribution:\n", sep = "")

    TAB <- round(rbind(
      x$mu[j, ],
      x$sigma[j, ],
      if (lambda_id[j]) lambda[j, ],
      if (nu_id[j]) nu[j, ]
    ), digits)

    rownames(TAB) <- c("mu", "sigma", if (lambda_id[j]) "lambda", if (nu_id[j]) "nu")
    print(TAB)
    cat("\n\n")
  }

  if (!is.null(x$Gamma)) {
    cat(crayon::cyan(sub("(.)", "\\U\\1", x$association, perl = TRUE),
                     " association matrix:\n\n", sep = ""))
    print(round(x$Gamma, digits), na.print = ".")
    cat("\n")
  }

  cat(
    "---",
    crayon::cyan("\nLoglik:"), round(x$logLik, digits), "|",
    crayon::cyan("AIC:"), round(x$AIC, digits), "|",
    crayon::cyan("BIC:"), round(x$BIC, digits)
  )

  invisible(x)

}



globalVariables(c("theo", "emp", "marg", "grid_x", "grid_y", "prob", "density"))
# Plot
#' Visualization of the fit of the BCS-NSM distributions
#'
#'
#'
#' @param x an object of class \code{bcsnsm}.
#' @param type character; specifies which graphical should be produced. The available options
#'     are \code{"response"} (default), \code{"margins"}, and \code{"epsilon"}.
#' @param outliers logical; used only when \code{type = "response"}. If \code{TRUE},
#'     possible outliers are highlighted in red.
#' @param alpha criterion according to the squared Mahalanobis distances that identifies a point as
#'     a possible outlier.
#' @param levels levels for contours plots.
#' @param panel A vector of the form \code{c(nr, nc)} with the number of rows and columns, where
#'     the figures will be drawn in an \code{nr-}by\code{-nc} array on the device.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
plot.bcsnsm <- function(x, type = c("response", "marginal", "epsilon"),
                       outliers = FALSE, alpha = 0.01,
                       panel = NULL,
                       levels = c(1, 1e-1, 1e-2, 1e-3, 1e-4), ...) {

  type <- match.arg(type, c("response", "marginal", "epsilon"))

  y <- x$y
  epsilon <- stats::residuals(x, "epsilon")

  n <- x$nobs
  d <- x$d
  margins <- x$margins
  copula <- x$copula
  eta <- x$eta

  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  mu <- x$mu
  sigma <- x$sigma
  lambda <- rep(NA, d)
  lambda[lambda_id] <- x$lambda
  nu <- rep(NA, d)
  nu[nu_id] <- x$nu

  gamma <- x$gamma
  digits <- 4L
  if (is.null(panel)) panel <- c(ceiling(d / 4), 4)

  association <- if (x$association == "non-associative") "nonassociative" else x$association
  Gamma <- get(association, envir = asNamespace("bcsnsm"))(d)$Gamma(x$gamma)

  mcopula <- make_copula(copula, eta / (1 - eta))

  dmv <- mcopula$dmv
  maha <- mcopula$maha
  dPSI <- mcopula$dPSI
  pPSI <- mcopula$pPSI
  qPSI <- mcopula$qPSI

  md <- Rfast::mahala(epsilon, rep(0L, d), Gamma)
  id_md <- maha(md, d) > 1 - alpha

  y_aux <- y
  y_names <- if (is.null(colnames(y))) paste0("y", 1:d) else colnames(y)

  ## Response contour plot -----------------------------------------------------------------------
  if (type == "response") {

    ### Diagonal plots
    diag_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)
      y <- GGally::eval_data_col(data, mapping$y)

      xid <- which(colSums(y_aux - x) == 0)

      dBCS <- function(x) {
        get(paste0("d", margins[xid]), envir = asNamespace("bcsnsm"))(x,
                                       mu = mu[xid],
                                       sigma = sigma[xid],
                                       lambda = lambda[xid],
                                       nu = nu[xid])
      }

      ggplot2::ggplot(data, mapping) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                colour = 1, fill = "white",
                                bins = ceiling(1 + 3.33 * log(n))) +
        ggplot2::geom_function(fun = dBCS, col = "#00AFBB")


    }

    ### Upper plots
    upper_func <- function(data, mapping, ...){

      x <- GGally::eval_data_col(data, mapping$x)
      y <- GGally::eval_data_col(data, mapping$y)

      i <- which(colSums(y_aux - x) == 0)
      j <- which(colSums(y_aux - y) == 0)

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

      i <- which(colSums(y_aux - x) == 0)
      j <- which(colSums(y_aux - y) == 0)

      Gamma_aux <- matrix(c(1, Gamma[i, j], Gamma[i, j], 1), 2, 2)

      mu_aux <- mu[c(i, j)]
      sigma_aux <- sigma[c(i, j)]
      lambda_aux <- lambda[c(i, j)]
      nu_aux <- nu[c(i, j)]

      grid <- expand.grid(seq(min(x) - 10, max(x) + 10, length.out = 200),
                          seq(min(y) - 10, max(y) + 10, length.out = 200))

      data_aux <- data.frame(grid_x = grid[, 1],
                             grid_y = grid[, 2],
                             prob = dbcsnsm(cbind(grid[, 1], grid[, 2]),
                                            mu_aux, sigma_aux,
                                            lambda_aux, nu_aux, Gamma_aux,
                                            eta = eta, copula = copula,
                                            margins = margins[c(i, j)]))

      data_aux2 <- data.frame(x = x[id_md], y = y[id_md])


      out <- ggplot2::ggplot() +
        ggplot2::geom_point(data = data, mapping = ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_contour(data = data_aux, mapping =  ggplot2::aes(x = grid_x, y = grid_y, z = prob),
                              col = "#00AFBB", breaks = levels) +
        ggplot2::labs(x = "", y = "")

      if (outliers) {
        out + ggplot2::geom_point(data = data_aux2, mapping = ggplot2::aes(x = x, y = y), col = "red")
      } else {
        out
      }



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

    GGally::ggpairs(as.data.frame(y), #ggplot2::aes(colour = gender),
                    upper = list(continuous = GGally::wrap(upper_func)),
                    lower = list(continuous = GGally::wrap(lower_func)),
                    diag = list(continuous = GGally::wrap(diag_func)),
                    legend = GGally::grab_legend(lab),
                    progress = FALSE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::theme(legend.position = "top")


  } else if (type == "marginal") {

    ## Marginal distributions ------------------------------------------------------------------------
    theo_CDF <- emp_CDF <- matrix(NA, n, ncol(y))
    ks <- vector()
    for(j in 1:d){

      theo_CDF[ , j] <- get(paste0("p", margins[j]), envir = asNamespace("bcsnsm"))(q = y[, j],
                                                     mu = mu[j],
                                                     sigma = sigma[j],
                                                     lambda = lambda[j],
                                                     nu = nu[j])

      ks[j] <- round(suppressWarnings(stats::ks.test(theo_CDF[, j], "punif")$p.value), 2)

      emp_CDF[, j] <- stats::ecdf(y[, j])(y[, j])

      id <- order(theo_CDF[, j])

      theo_CDF[, j] <- theo_CDF[id, j]
      emp_CDF[, j] <- emp_CDF[id, j]

    }

    dKS <- sfsmisc::KSd(n)
    aux_data <- data.frame(theo = c(theo_CDF),
                           emp = c(emp_CDF),
                           lower = c(emp_CDF - dKS),
                           upper = c(emp_CDF + dKS),
                           marg = factor(rep(colnames(y), each = n),
                                         levels = colnames(y)),
                           ks = rep(paste0("KS p-value: ", ks), each = n))

    positions <- data.frame(x = c(rbind(-1, theo_CDF, 2, apply(theo_CDF, 2, rev))),
                            y = c(rbind(-1 - dKS, theo_CDF - dKS, 2 + dKS, apply(theo_CDF, 2, rev) + dKS)),
                            marg = factor(rep(colnames(y), each = 2 * n + 2),
                                          levels = colnames(y)))

    ggplot2::ggplot() +
      ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, group = marg), positions, fill = "#cceff1") +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::geom_point(ggplot2::aes(x = theo, y = emp, group = marg), aux_data, size = 1) +
      ggplot2::geom_abline(intercept = 0, slope = 1, col = "#00AFBB", lty = 1) +
      ggplot2::geom_text(ggplot2::aes(x = 0.8, y = 0, label = ks, group = marg), aux_data, size = 3) +
      ggplot2::facet_wrap(~ marg, nrow = panel[1], ncol = panel[2]) +
      ggplot2::labs(x = "Fitted distribution function",
                    y = "Empirical distribution function", size = 1)


  } else if (type == "epsilon") {
    ## Epsilon contour plot ----------------------------------------------------------------------


    ### Diagonal plots
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

      mu_aux <- mu[c(i, j)]
      sigma_aux <- sigma[c(i, j)]
      lambda_aux <- lambda[c(i, j)]
      nu_aux <- nu[c(i, j)]

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

    GGally::ggpairs(as.data.frame(epsilon), #ggplot2::aes(colour = gender),
                    upper = list(continuous = GGally::wrap(upper_func)),
                    lower = list(continuous = GGally::wrap(lower_func)),
                    diag = list(continuous = GGally::wrap(diag_func)),
                    legend = GGally::grab_legend(lab),
                    progress = FALSE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::theme(legend.position = "top")


  }

}
