#' @title Multivariate BCS-NSM Regression
#' 
#' @description
#' Fits a multivariate regression model for positive continuous data using the BCS-NSM class.
#'     The marginal distributions are chosen from the Box-Cox symmetric (BCS) family 
#'     (see \code{\link{bcs}}), and the dependence structure is modeled via a normal scale mixture
#'     (NSM) copula.
#'
#' @param formula a symbolic description of the model to be fitted to each marginal distribution.
#'     For instance, \code{formula = c(y1 + y2 + y3 ~ x1 + x2 | x1 + x3 | x2 + x3)} fits a
#'     3-variate BCS-NSM regression that models \code{y1} as a function of \code{x1} and \code{x2};
#'     \code{y2} as a function of \code{x1} and \code{x3}; and \code{y3} as a function of \code{x2}
#'     and \code{x3}. See details below.
#' @param data an optional data frame containing the variables in the formula. By default the
#'     variables are taken from environment(formula).
#' @param subset an optional vector specifying a subset of index observations to be used in the fitting
#'     process. (See additional details about how this argument interacts with data-dependent bases
#'     in the ‘Details’ section of the \code{\link[stats]{model.frame}} documentation.)
#' @param na.action a function which indicates what should happen when the data contain \code{NAs}.
#'     The default is set by the \code{na.action} setting of \link[base]{options}, and is
#'     \code{\link[stats]{na.fail}} if that is unset. The ‘factory-fresh’ default is
#'     \code{\link[stats]{na.omit}}. Another possible value is \code{NULL}, no action.
#'     Value \code{\link[stats]{na.exclude}} can be useful.
#' @param margins a character or a character vector; specifies the marginal BCS distributions. If all
#'     BCS margins are the same, it is sufficient to enter only one character. A table with the
#'     current available BCS distributions can be seen in \code{\link{bcs}}.
#' @param links a character or a character vector; specifies the link functions for the \code{mu}
#'     regression structure. The links \code{"log"} (default) and \code{"identity"} are currently
#'     available. If all links are the same, it is sufficient to enter only one character.
#' @param association one of \code{"unstructured"} (default), \code{"uniform"}, or
#'     \code{"nonassociative"}, which specify the association matrix of the model.
#' @param copula character; informs which normal scale mixture distribution
#'     should be used to generate the NSM copula. Currently,
#'     the copulas available are: Gaussian (\code{"gaussian"}),
#'     Student's \emph{t} (\code{"t"}), slash (\code{"slash"}), and hyperbolic (\code{"hyp"}).
#' @param eta possible extra parameter induced by the copula, parameterized to take values at
#'     \code{(0, 1)}. For instance, under the \emph{t} copula, \code{eta} specifies
#'     \code{eta} / (1 - \code{eta}) degrees of freedom. To estimate the value of \code{eta}
#'     via profile log-likelihood, use the function \code{\link{choose_copula}}.
#' @param control  a list of optimization control arguments specified via \code{\link{control_fit}}.
#' @param ... further optimization control arguments passed to \code{\link{control_fit}}.
#'
#' @return The \code{bcsnsm} function returns an object of class "\code{bcsnsm}",
#'     which consists of a list with the following components:
#' \describe{
#'   \item{coefficients}{a list with the estimated regression coefficients for each marginal distribution.}
#'   \item{fitted.values}{the fitted medians for each marginal distribution.}
#'   \item{marginal.pars}{a list with the fitted values for the marginal parameters \code{mu}, \code{sigma},
#'       \code{lambda}, and \code{nu}.}
#'   \item{margins}{a character vector with the specified marginal BCS distributions.}
#'   \item{links}{a character vector with the specified link functions.}
#'   \item{copula, eta}{\code{"copula"} is a character which informs which normal scale mixture distribution
#'       was used to generate the NSM copula and \code{"eta"} is the possible extra parameter associated with
#'       the copula, parameterized to take values at (0, 1).}
#'   \item{gamma}{the estimated parameters of the association matrix, if any.}
#'   \item{association}{structure of the association matrix. It can be one of \code{"non-associative"},
#'       \code{"unstructured"}, or \code{"uniform"}.}
#'   \item{logLik}{log-likelihood of the fitted model.}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood estimator of the model parameters vector.
#'       Specifically, the inverse of the observed information matrix, obtained via numeric Hessian matrix.}
#'   \item{optim_params}{control optimization parameters used by \code{\link{control_fit}}.}
#'   \item{y}{the response matrix used.}
#'   \item{nobs,d}{the number of observations in the sample and the dimension of the response variable, respectively.}
#'   \item{call}{ the function call.}
#'   \item{formula}{the formula used to specify the model in \code{bcsnsm}.}
#'   \item{terms}{a list with the terms for the "\code{mu}" submodel for each marginal.}
#'  }
#'  
#' @details 
#' Let \eqn{\bold{Y}_1, \ldots, \bold{Y}_n} be \eqn{n} independent positive continuous random 
#'     vectors, where each \eqn{\bold{Y}_i = (Y_{i1}, \ldots, Y_{id})^\top} follows a BCS-NSM 
#'     distribution with marginal parameters \eqn{\bold{\mu}_i = (\mu_{i1}, \ldots, \mu_{id})^\top},
#'     \eqn{\bold{\sigma} = (\sigma_1, \ldots, \sigma_d)^\top}, \eqn{\bold{\lambda} = 
#'     (\lambda_1, \ldots, \lambda_d)^\top}; marginal density generating functions 
#'     \eqn{\bold{r}(u) = (r_1(u), \ldots, r_d(u))^\top}, \eqn{u \geqslant 0}, and NSM copula 
#'     generated by the mixing CDF \eqn{H} with association matrix \eqn{\bold{\Gamma}}. 
#'     The class of the multivariate BCS-NSM regression models is defined as
#'     
#'     \eqn{
#'     \begin{array}{rl}
#'     \bold{Y}_i  & \sim \textrm{BCS--NSM}_d(\bold{\mu}_i, \bold{\sigma}, \bold{\lambda}, \bold{\Gamma}; \bold{r}, H),\\
#'     \bold{\mu}_i & = \left(g_1^{-1}(\bold{x}_{i1}^\top \bold{\beta}_1), \ldots, g_d^{-1}(\bold{x}_{id}^\top \bold{\beta}_d)\right)^\top, 
#'     \end{array}
#'     }
#'     
#'     where \eqn{\bold{\beta}_j} is a \eqn{b_j \times 1} vector of unknown regression coefficients; 
#'     \eqn{\bold{x}_{ij}} is a \eqn{b_j \times 1} vector with the observations on \eqn{b_j} 
#'     explanatory variables, which are assumed to be fixed and known; and \eqn{g_j : \mathbb{R}_+ 
#'     \longrightarrow \mathbb{R}} is a strictly monotone and at least twice differentiable link function.
#'     
#'  The \code{formula} argument uses functionalities inherited from \code{Formula}
#'      package (Zeileis and Croissant, 2010). A basic formula for specifying a 3-variate BCS-NSM
#'      regression model is, for instance, \code{y1 + y2 + y3 ~ x1 + x2 + x3}, which specifies the
#'      same regression structure for the three marginal distributions as a function of the
#'      explanatory variables \code{x1, x2} and \code{x3}. However, using the Formula package
#'      allows different regression structures to be specified for each marginal distribution using
#'      the "\code{|}" operator. For instance, \code{formula = c(y1 + y2 + y3 ~ x1 + x2 | x1 + x3 |
#'      x2 + x3)} fits a 3-variate BCS-NSM regression that models \code{y1} as a function of
#'      \code{x1} and \code{x2}; \code{y2} as a function of \code{x1} and \code{x3}; and \code{y3}
#'      as a function of \code{x2} and \code{x3}.
#'  
#'     The class of the BCS-NSM multivariate regression models  includes models with various 
#'     possibilities of marginal and dependence settings. For instance, each margin can follow a 
#'     different BCS distribution. This property provides the flexibility to accommodate marginal 
#'     distributions with different behaviors, such as different types of skewness and tail-heaviness.
#'
#'     The NSM copula is characterized by the mixing cumulative distribution function (CDF) \eqn{H}. 
#'     The \code{bcsnsm} function provides the Gaussian, \eqn{t}, slash, and hyperbolic copulas, 
#'     specified by the \code{copula} argument. Note that some NSM copulas also have an extra 
#'     parameter, say \eqn{\delta}, which is assumed to be known and fixed in the BCS-NSM 
#'     multivariate regression. This parameter is generally positive, but we consider the 
#'     parameterization \eqn{\eta = \delta / (1 + \delta) \in (0, 1)} to make it easier to choose 
#'     its value in practice.
#'     
#'     The association matrix \eqn{\Gamma} plays an important role in describing the dependence.
#'     Let \eqn{\gamma_{ij}} be the element of the \emph{i}th row and \emph{j}th column of the 
#'     association matrix \eqn{\bold{\Gamma}}, where \eqn{\gamma_{ii} = 1}, for 
#'     \emph{i, j = 1, ..., d}. Thus, Kendall's tau of the \eqn{i}th and \eqn{j}th margins is given by
#'
#'     \deqn{\tau_{ij} = \dfrac{2 \; \textrm{asin}(\gamma_{ij})}{\pi},}
#'
#'     and it is invariant under the copula. Since Kendall’s tau is a rank correlation coefficient, 
#'     it measures any monotonous association between two random variables. 
#'       
#'  Parameter estimation is performed by directly maximizing the log-likelihood function through
#'      the \code{\link[stats]{optim}} function. The \code{\link{choose_copula}} function can be used
#'      to select the possible extra parameter (\code{eta}) of the NSM copula.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @examples
#' \donttest{#' # Data set: macronutrients (for description run ?macronutrients)
#' # Consider modeling the animal and plant protein regarding income level
#' 
#' ## Response distribution:
#' mvplot(macronutrients[, 1:2])
#' 
#' ## Animal and plant protein by income level:
#' boxplot(animal_p ~ income, macronutrients)
#' boxplot(plant_p ~ income, macronutrients)
#' 
#' ## Step 1: Reference model
#' fit0 <- bcsnsm(animal_p + plant_p ~ income | income, data = macronutrients)
#' fit0
#' 
#' ## Marginal quantile residuals of the reference model
#' plot(fit0, "marginal", panel = c(2, 3))
#' 
#' ## Step 2: Improve the fit on margins
#' fit_gaussian <- bcsnsm(animal_p + plant_p ~ income | income, data = macronutrients,
#'                        margins = c("bcloii", "bct"))
#' 
#' ## Marginal quantile residuals of the improved model
#' plot(fit_gaussian, "marginal", panel = c(2, 3))
#' 
#' ## Step 3: Choose a copula
#' 
#' ### Choosing the extra parameter of the t copula (it can be slow)
#' fit_t <- choose_copula(fit_gaussian, copula = "t")
#' plot(fit_t)
#' 
#' ## Comparing with AIC
#' AIC(fit_gaussian)
#' AIC(fit_t[["0.9"]])
#' 
#' ## Overall quantile residuals
#' plot(fit_gaussian)
#' plot(fit_t[["0.9"]])
#' 
#' ## The epsilon's transformations
#' plot(fit_gaussian, "epsilon")
#' plot(fit_t[["0.9"]], "epsilon")}
#' @references
#'  Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, \bold{101}, 321--344.
#'      
#'  Medeiros, R. M. R., and Ferrari, S. L. P. (2025). Multivariate Box-Cox symmetric regression 
#'      models generated by a normal scale mixture copula.
#'
#'  Zeileis A., and Croissant, Y. (2010). Extended model formulas in
#'  R: multiple parts and multiple responses. \emph{Journal of Statistical Software},
#'  \bold{34}, 1–-13.
#'
#'
#' @export
bcsnsm <- function(formula, data, subset, na.action,
                   margins = "bcno",
                   links = "log",
                   association = c("unstructured", "uniform", "nonassociative"),
                   copula = c("gaussian", "t", "slash", "hyp"),
                   eta = NULL,
                   control = control_fit(...), ...)
{

  ## Function call ----
  cl <- match.call()

  ## Model frame ----
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) data <- environment(formula[[1]])

  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE

  f <- Formula::as.Formula(formula)
  if (length(f)[1L] > 2L) {
    f <- Formula::Formula(formula(f, lhs = 1))
    warning("formula must not have more than one LHS part")
  }

  y <- as.matrix(Formula::model.part(f, data = data, lhs = 1L), rownames.force = TRUE)
  n <- dim(y)[1L]
  d <- dim(y)[2L]

  if (any(y < 0) )
    stop("there are responses with negative values")

  if (n < 1)
    stop("empty model")

  simple_formula <- length(f)[2L] < 2L
  if (simple_formula) {

    resp <- as.character(f[2L])
    x <- paste(rep(as.character(f[3L]), d), collapse = " | ")
    f <- Formula::as.Formula(paste(resp, "~", x))

  }

  mf[[1L]] <- as.name("model.frame")
  mf$formula <- f
  mf <- eval(mf, parent.frame())

  ## Extract terms and model matrices ----
  mt <- stats::terms(f, data = data)
  mtX <- lapply(1:d, function(j) stats::terms(f, data = mf, rhs = j))

  X <- lapply(1:d, function(j) stats::model.matrix(f, data = mf, rhs = j))

  ## Marginal information ----
  margins <- as.vector(matrix(margins, 1, d))
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)
  links <- as.vector(matrix(links, 1, d))

  ## Copula information ----
  copula <- match.arg(copula, c("gaussian", "t", "slash", "hyp"))
  if (!is.null(eta)) {
    if (eta <= 0 | eta >= 1)
      stop("eta must be a value on (0, 1)")
  }

  mcopula <- make_copula(copula, eta / (1 - eta))
  qPSI <- mcopula$qPSI
  dgf <- mcopula$dgf

  association <- match.arg(association, c("unstructured", "uniform", "nonassociative"))
  association <- get(association, envir = asNamespace("bcsnsm"))(d)

  ## Optimization ----

  ### Specific optim parameters ----
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  trace <- control$trace
  inits <- control$inits

  control$method <- control$hessian <- control$inits <- NULL

  ### Initial values ----
  if (is.null(inits)) {

    beta <- vector("list", d)
    sigma <- rep(NA, d)
    lambda <- rep(NA, d)
    nu <- rep(NA, d)

    for (j in 1:d) {

      if (margins[j] %in% c("bcloi", "lloi", "bcpe", "lpe", "bchp", "lhp", "bcla", "lla")){
        faux <- gamlss.dist::BCPE
      } else if (margins[j] %in% c("bcno", "lno", "bcloii", "lloii")){
        faux <- gamlss.dist::BCCG
      } else{
        faux <- gamlss.dist::BCT
      }

      fj <- stats::as.formula(paste(colnames(y)[j], "~",
                                    as.character(stats::formula(stats::terms(f, lhs = 0, rhs = j)))[2L]))
      gamlss_fit <- switch (links[j],

                            log = {
                              suppressWarnings(gamlss::gamlss(fj, data = data,
                                                              family = faux(mu.link = "log",sigma.link = "identity",
                                                                            nu.link = "identity"),
                                                              trace = FALSE))
                            },

                            identity = {
                              suppressWarnings(gamlss::gamlss(fj,
                                                              family = faux(mu.link = "identity", sigma.link = "identity",
                                                                            nu.link = "identity"),
                                                              trace = FALSE))
                            }
      )

      beta[[j]] <- stats::coef(gamlss_fit, "mu")
      sigma[j] <- stats::coef(gamlss_fit, "sigma")
      if (lambda_id[j]) lambda[j] <- stats::coef(gamlss_fit, "nu")
      if (nu_id[j]) nu[j] <- min(exp(stats::coef(gamlss_fit, "tau")), 20)

    }

    inits <- list(beta = beta, sigma = sigma, lambda = if (any(lambda_id)) lambda[lambda_id] else NULL,
                  nu = if (any(nu_id)) nu[nu_id] else NULL, gamma = association$start(y))

  } 
  
  ### Parameter id ----
  beta_breaks <- c(0L, cumsum(unlist(lapply(inits$beta, length))))

  beta_len <- beta_breaks[length(beta_breaks)]
  lambda_len <- sum(lambda_id)
  nu_len <- sum(nu_id)
  gamma_len <- association$npar

  par_id <- list(beta = Map(function(x, y) (x + 1):y, beta_breaks[-length(beta_breaks)], beta_breaks[-1]),
                 sigma = beta_len + 1:d,
                 lambda = if (any(lambda_id)) 1:lambda_len + beta_len + d else NULL,
                 nu = if (any(nu_id)) 1:nu_len + beta_len + lambda_len + d else NULL,
                 gamma = if(gamma_len > 0) 1:gamma_len + beta_len + d + lambda_len + nu_len)

  beta_class <- factor(unlist(lapply(seq_along(par_id$beta), function(i) {
    rep(colnames(y)[i], length(par_id$beta[[i]]))})), levels = colnames(y))

  ### Log-likelihood ----
  EPS <- .Machine$double.eps^(1/1.6)

  pred <- function(X, link, beta){

    stats::make.link(link)$linkinv(X%*%beta)

  }

  inits <- c(unlist(inits$beta), unlist(inits)[-(1:beta_len)])
  inits <- inits[!is.na(inits)]

  ll <- function(theta){

    ### Parameter setting
    beta <- split(theta[1:beta_len], beta_class)

    mu <- matrix(unlist(Map(pred, X, as.list(links), beta)), ncol = d)
    sigma <- theta[par_id$sigma]
    lambda <- rep(NA, d)
    lambda[lambda_id] <- theta[par_id$lambda]
    nu <- rep(NA, d)
    nu[nu_id] <- theta[par_id$nu]
    gamma <- theta[par_id$gamma]

    ### Association matrix
    gamma_id <- FALSE
    Gamma <- association$Gamma(gamma)
    dec <- tryCatch(Rfast::cholesky(Gamma), error = function(e) e)
    if (inherits(dec, "error") | any(is.nan(dec))) dec <- NULL
    if (length(gamma) > 0L) gamma_id <- any(gamma < association$lower | gamma > association$upper)

    ### Out
    if (any(!is.finite(mu)) | any(!is.finite(sigma)) | any(!is.finite(lambda[lambda_id])) |
        any(mu < 0) | any(sigma < 0) | any(nu[nu_id] < 1) | gamma_id | is.null(dec)) {

      -Inf

    }else {

      epsilon <- log_f <- matrix(NA, n, d)
      for(j in 1:d){

        epsilon[, j] <- qPSI(pmin(pmax(get(paste0("p", margins[j]), envir = asNamespace("bcsnsm"))(q = y[, j],
                                                                    mu = mu[, j],
                                                                    sigma = sigma[j],
                                                                    lambda = lambda[j],
                                                                    nu = nu[j]), EPS), 1 - EPS))

        log_f[, j] <- get(paste0("d", margins[j]), envir = asNamespace("bcsnsm"))(x = y[, j],
                                                   mu = mu[, j],
                                                   sigma = sigma[j],
                                                   lambda = lambda[j],
                                                   nu = nu[j], log = TRUE)
      }

      tmp <- backsolve(dec, t(epsilon), transpose = TRUE)
      rss <- colSums(tmp^2)
      -n * sum(log(diag(dec))) + sum(dgf(rss, d, log = TRUE)) -
        sum(dgf(epsilon^2, 1L, log = TRUE)) + sum(log_f)

    }

  }

  ### optim out ----
  opt <- stats::optim(par = inits,
                      fn = ll,
                      method = method,
                      control = control,
                      hessian = hessian)

  ### Asymptotic covariance matrix -----
  if (hessian) {
    J <- -opt$hessian #-numDeriv::hessian(ll, opt$par)
    vcov <- try(solve(J), silent = TRUE)
    vcov <- if (unique(grepl("Error", vcov))) matrix(NA, nrow = length(opt$par), ncol = length(opt$par)) else vcov
  }

  if (opt$convergence > 0)
    warning(cat("optimization failed to converge\n"))

  opt$par_id <- par_id
  opt$inits <- inits

  ### Marginal estimates ----
  beta <- split(opt$par[1:beta_len], beta_class)

  mu <- matrix(unlist(Map(pred, X, as.list(links), beta)), ncol = d)
  colnames(mu) <- colnames(y)

  sigma <- opt$par[par_id$sigma]
  names(sigma) <- colnames(y)

  lambda <- rep(NA, d)
  names(lambda) <- colnames(y)
  lambda[lambda_id] <- opt$par[par_id$lambda]

  nu <- rep(NA, d)
  names(nu) <- colnames(y)
  nu[nu_id] <- opt$par[par_id$nu]

  ### Fitted values (fitted marginal medians)
  fv <- matrix(unlist(Map(function(mu, sigma, lambda, nu, margin){
    get(paste0("q", margin), envir = asNamespace("bcsnsm"))(p = 0.5, mu = mu, sigma = sigma, lambda = lambda, nu = nu)
  }, mu, sigma, lambda, nu, margins)), ncol = d)
  colnames(fv) <- colnames(y)

  ### Out list
  out <- list(
    coefficients = beta,
    fitted.values = fv,
    marginal.pars =  list(mu = mu, sigma = sigma,
                          lambda = if (any(lambda_id)) lambda[lambda_id] else NULL,
                          nu = if (any(nu_id)) nu[nu_id] else NULL),
    margins = margins,
    links = links,
    copula = copula,
    eta = eta,
    gamma = if (association$npar > 0L) opt$par[par_id$gamma] else NULL,
    association = association$name,
    logLik = opt$value,
    vcov = if (hessian) vcov else NULL,
    optim_pars = opt,
    y = y,
    nobs = n,
    d = d
  )

  ## Further model information ----
  out$call <- cl
  out$formula <- Formula::as.Formula(f)
  out$terms <- list(margins = mtX, full = mt)

  class(out) <- "bcsnsm"
  out

}


# Print method
#' @rdname bcsnsm
#' 
#' @param x a fitted model object of class \code{"bcsnsm"}.
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values.
#' 
#' @export
print.bcsnsm <- function(x, digits = max(3, getOption("digits") - 3), ...)
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
  
  cat(crayon::cyan("\nMultivariate BCS-NSM Regression\n"))
  
  cat(paste0(crayon::cyan("\nCall:\n"), 
             paste(deparse(x$call, width.cutoff = floor(getOption("width") * 0.7)), collapse = "\n")))
  
  
  if (x$optim_pars$convergence > 0) {
    cat("\nmodel did not converge\n")
  } else {
    
    cat(crayon::cyan("\n\nRegression coefficients:\n"))
    
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
    
    cat(crayon::cyan("\nFurther marginal parameters:\n"))
    print(t(tab), na.print = "-")
    
    if (length(gamma) > 0) {
      
      cat(crayon::cyan("\nAssociation matrix:\n", sep = ""))
      
      if (x$association == "non-associative"){
        Gamma <- diag(d)
      }else{
        Gamma <- get(tolower(x$association), envir = asNamespace("bcsnsm"))(d)$Gamma(round(gamma, digits))
      }
      
      Gamma[upper.tri(Gamma)] <- diag(Gamma) <- NA
      colnames(Gamma) <- rownames(Gamma) <- colnames(y)
      print(Gamma, na.print = ".")
    }
    
    df <- n - length(x$optim_pars$par) - as.numeric(!is.null(x$eta))
    
    cat(
      "\n---",
      crayon::cyan("\nMargins:"), margins,
      crayon::cyan("\nCopula:"), if (is.null(eta)) copula_name else paste0(copula_name, "(", round(eta, 2), ")"),
      crayon::cyan("\nAssociation:"), gamma_name,
      crayon::cyan("\nlogLik:"), round(x$logLik, digits), "on", df, "df", "\n"
    )
  }
  
  invisible(x)
}