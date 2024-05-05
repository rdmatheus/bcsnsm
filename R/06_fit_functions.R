#' @title BCS-NSM Regression for Multivariate Positive Data
#'
#' @description Fit a multivariate BCS-NSM regression via maximum likelihood for positive continuous
#'     data. Specify the marginal distributions within the class of Box-Cox symmetric (BCS) distributions
#'     (see \code{\link{bcs}}) and the dependence structure regarding the copula of a normal scale
#'     mixture distribution (see \code{\link{bcsnsm.dist}}).
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
#' @param control  a list of optimization control arguments specified via \code{\link{control_fit}}  .
#' @param ... further optimization control arguments passed to \code{\link{control_fit}}.
#'
#' @return The \code{bcsnsmreg} function returns an object of class "\code{bcsnsmreg}",
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
#'  
#'  The \code{formula} argument of the \code{formula} method uses functionalities inherited from
#'       package (Zeileis and Croissant, 2010). A basic formula for specifying a 3-variate BCS-NSM
#'       regression model is, for instance, \code{y1 + y2 + y3 ~ x1 + x2 + x3}, which specifies the
#'       same regression structure for the three marginal distributions as a function of the
#'       explanatory variables \code{x1, x2} and \code{x3}. However, using the Formula package
#'       allows different regression structures to be specified for each marginal distribution using
#'       the "\code{|}" operator. For instance, \code{formula = c(y1 + y2 + y3 ~ x1 + x2 | x1 + x3 |
#'       x2 + x3)} fits a 3-variate BCS-NSM regression that models \code{y1} as a function of
#'       \code{x1} and \code{x2}; \code{y2} as a function of \code{x1} and \code{x3}; and \code{y3}
#'       as a function of \code{x2} and \code{x3}.
#'
#'    For the \code{default} method, the \code{bcnsm} function fits a BCS-NSM distribution without
#'       regression structures. It returns a list with the same components as the \code{formula}
#'       method, except for the \code{coefficients}, \code{links}, \code{formula}, and \code{terms}
#'       components, which are set to NULL. The \code{y} argument represents a \code{data.frame} or
#'       a matrix containing non-negative observations.
#'
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' \dontrun{
#'
#' # Winscosin Breast Cancer Dataset
#' ?wdbc
#'
#' # Training set index
#' set.seed(123)
#' id <- sample(1:nrow(wdbc), 0.7 * nrow(wdbc))
#'
#' # Reference model
#' fit0 <- bcnsmreg(cbind(Texture, Area, Smoothness, Compactness, Concavity) ~ Diagnosis,
#'                  data = wdbc, subset = id)
#' fit0
#'
#' ## Marginal quantile residuals of the reference model
#' plot(fit0, "marginal", panel = c(2, 3))
#'
#' # Improve the fit on margins
#' fit_gaussian <- bcnsmreg(cbind(Texture, Area, Smoothness, Compactness, Concavity) ~ Diagnosis,
#'                          data = wdbc, subset = id, margins = c("lt", "lt", "lno", "lpe", "bct"))
#'
#' ## Marginal quantile residuals of the improved model
#' plot(fit_gaussian, "marginal", panel = c(2, 3))
#'
#'
#' ## Summary
#' summary(fit_gaussian)
#'
#' ## Overall quantile residuals of the final model
#' plot(fit_gaussian)
#'
#' ## The epsilon's transformations
#' plot(fit_gaussian, "epsilon")
#' }
#'
#' @references
#'
#'  Ferrari, S. L., and Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#'  Zeileis A., and Croissant, Y. (2010). Extended model formulas in
#'  R: multiple parts and multiple responses. \emph{Journal of Statistical Software},
#'  \bold{34}, 1–13.
#'
#'
#' @export
#'
bcsnsmreg <- function(formula, data, subset, na.action,
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
  association <- get(association, envir = parent.frame())(d)

  ## Optimization ----

  ### Specific optim parameters ----
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  trace <- control$trace

  control$method <- control$hessian <- control$inits <- NULL

  ### Initial values ----
  if (is.null(control$inits)) {

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

  } else {
    inits <- control$inits
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

        epsilon[, j] <- qPSI(pmin(pmax(get(paste0("p", margins[j]), envir = parent.frame())(q = y[, j],
                                                                    mu = mu[, j],
                                                                    sigma = sigma[j],
                                                                    lambda = lambda[j],
                                                                    nu = nu[j]), EPS), 1 - EPS))

        log_f[, j] <- get(paste0("d", margins[j]), envir = parent.frame())(x = y[, j],
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
    get(paste0("q", margin), envir = parent.frame())(p = 0.5, mu = mu, sigma = sigma, lambda = lambda, nu = nu)
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

  class(out) <- "bcsnsmreg"
  out

}


#' @title BCS-NSM Fit for Multivariate Positive Data
#'
#' @description Fit a multivariate BCS-NSM distribution via maximum likelihood for positive continuous
#'     data. Specify the marginal distributions within the class of Box-Cox symmetric (BCS) distributions
#'     (see \code{\link{bcs}}) and the dependence structure regarding the copula of a normal scale
#'     mixture distribution (see \code{\link{bcsnsm.dist}}).
#'
#' @param y a matrix with the sample observations.
#' @param association one of \code{"unstructured"} (default), \code{"uniform"}, or
#'     \code{"nonassociative"}, which specify the association matrix of the model.
#' @param copula character; informs which normal scale mixture distribution should be used to generate
#'     the NSM copula. Currently, the copulas available are: Gaussian (\code{"gaussian"}),
#'     Student's \emph{t} (\code{"t"}), slash (\code{"slash"}), and hyperbolic (\code{"hyp"}).
#' @param eta possible extra parameter induced by the copula, parameterized to take values at
#'     \code{(0, 1)}. For instance, under the \emph{t} copula, \code{eta} specifies
#'     \code{eta} / (1 - \code{eta}) degrees of freedom. To estimate the value of \code{eta}
#'     via profile log-likelihood, use the function \code{\link{choose_copula}}.
#' @param margins a character or a character vector; specifies the marginal BCS distributions.
#'     If all BCS margins are the same, it is sufficient to enter only one character. A table with
#'     the current available BCS distributions can be seen in \code{\link{bcs}}.
#' @param control a list of control arguments specified via \code{\link{control_fit}}.
#' @param ... further arguments passed to \code{\link{control_fit}}.
#'
#' @return The \code{bcsnsm} function returns an object of class \code{"bcsnsm"},
#'     which consists of a list with the following components:
#' \describe{
#'   \item{mu, sigma, lambda, nu}{vectors with the estimated values for the marginal parameters
#'       \code{mu}, \code{sigma}, \code{lambda}, and \code{nu}, respectively.}
#'   \item{gamma}{the estimated parameters of the association matrix, if any.}
#'   \item{margins}{a character vector with the specified marginal BCS distributions.}
#'   \item{association}{structure of the association matrix. It can be one of \code{"non-associative"},
#'       \code{"unstructured"}, or \code{"uniform"}.}
#'   \item{copula, eta}{\code{"copula"} is a character which informs which normal scale mixture distribution
#'       was used to generate the NSM copula and \code{"eta"} is the possible extra parameter associated with
#'       the copula.}
#'   \item{logLik}{log-likelihood of the fitted model.}
#'   \item{vcov}{asymptotic covariance matrix of the maximum likelihood estimator of the model parameters vector.
#'       Specifically, the inverse of the observed information matrix, obtained via numeric Hessian matrix.}
#'   \item{y}{the response matrix.}
#'   \item{optim_pars}{control optimization parameters used by \code{\link{control_fit}}.}
#'   \item{nobs,d}{the number of observations in the sample and the dimension of the response variable, respectively.}
#'   \item{call}{ the function call.}
#'  }
#'
#' @references Vanegas, L. H., and Paula, G. A. (2016). Log-symmetric distributions: statistical properties and
#'     parameter estimation. \emph{Brazilian Journal of Probability and Statistics}, 30, 196-220.
#'
#' Ferrari, S. L. P., and Fumes, G. (2017). Box-Cox symmetric distributions and applications to
#'     nutritional data. \emph{AStA Advances in Statistical Analysis}, 101, 321-344.
#'
#' Medeiros, R. M. R. de, and Ferrari, S. L. P. (2023). Multivariate Box-Cox symmetric distributions
#'     generated by a normal scale mixture copula.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Random sampling
#'
#' ### Sample size and dimension
#' n <- 200
#' d <- 4
#'
#' ### Association matrix
#' Gamma <- matrix(0.8, d, d)
#' diag(Gamma) <- 1
#'
#' ### Marginal specifications
#'
#' margins <- c("bcno", "lt", "bct", "lno")
#'
#' mu <- c(19, 20, 15, 20)
#' sigma <- c(0.2, 0.3, 0.4, 0.3)
#' lambda <- c(-1, NA, 1.6, NA)
#' nu <- c(NA, 4, 8, NA)
#'
#' ### Copula
#' copula <- "t"
#' eta <- 0.8
#'
#' ### Observations
#' set.seed(123) # Seed for reproducibility
#' y <- rbcsnsm(n, mu, sigma, lambda, nu, Gamma, copula, eta, margins)
#' colnames(y) <- c("y1", "y2", "y3", "y4")
#'
#' ### Visualization (based on graphics::pairs functions)
#' mvplot(y)
#'
#' ## Fit with Gaussian copula and uniform structure
#' fit <- bcsnsm(y, association = "uniform", copula = "gaussian",
#'               margins = c("bcno", "lt", "bct", "lno"))
#'
#' class(fit)
#' methods(class = "bcnsm")
#'
#' # Fit summaries
#' fit
#' summary(fit)
#'
#' # Fit visualization
#'
#' ## Bivariate fit
#' plot(fit)
#'
#' ## Marginal fit
#' plot(fit, type = "margins")
#'
#' ## Transformed vectors
#' plot(fit, "epsilon")
#'
#' # Choose the value of the extra parameter of the t copula (it can be slow)
#' fit_t <- choose_copula(fit, grid = 1:8, copula = "t")
#'
#' ## Final fit
#' final_fit <- fit_t[[4]]
#'
#' final_fit
#' plot(final_fit)
#' plot(final_fit, type = "margins")
#' }
bcsnsm <- function(y,
                   margins = "bcno",
                   association = c("unstructured", "uniform", "nonassociative"),
                   copula = c("gaussian", "t", "slash", "hyp"),
                   eta = NULL,
                   control = control_fit(...), ...)
{

  ## Function call ----
  cl <- match.call()

  if (is.vector(y)) y <- matrix(y, ncol = length(y))
  if (is.data.frame(y)) y <- as.matrix(y)

  n <- dim(y)[1L]
  d <- dim(y)[2L]

  if (any(y < 0) )
    stop("there are responses with negative values")

  if (n < 1)
    stop("empty model")


  ## Marginal information ----
  margins <- as.vector(matrix(margins, 1, d))
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

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
  association <- get(association, envir = parent.frame())(d)

  ## Optimization ----

  ### Specific optim parameters ----
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  trace <- control$trace

  control$method <- control$hessian <- control$inits <- NULL

  ### Initial values ----
  if (is.null(control$inits)) {

    mu0 <- sigma0 <- lambda0 <- nu0 <- rep(NA, d)

    for (j in 1:d) {
      marg_inits <- get(margins[j], envir = parent.frame())(y[, j])$start(y[, j])

      mu0[j] <- marg_inits[1]
      sigma0[j] <- marg_inits[2]
      lambda0[j] <- if (lambda_id[j]) marg_inits[3] else  NA
      if (nu_id[j]) nu0[j] <- utils::tail(marg_inits, 1)
    }

    inits <- c(mu0, sigma0, lambda0[lambda_id], nu0[nu_id], association$start(y))

  } else {

    inits <- control$inits

  }

  ### Parameter id ----
  par_id <- list(
    mu = 1:d,
    sigma = 1:d + d,
    lambda = if (any(lambda_id)) 1:sum(lambda_id) + 2 * d else NULL,
    nu = if (any(nu_id)) 1:sum(nu_id) + sum(lambda_id) + 2 * d else NULL,
    gamma = if (association$npar > 0) 1:association$npar + sum(nu_id) + sum(lambda_id) + 2 * d else NULL
  )

  ### Log-likelihood ----
  EPS <- .Machine$double.eps^(1/1.5)


  ll <- function(theta){

    ### Parameter setting
    mu <- theta[par_id$mu]
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
    if (inherits(dec, "error")) dec <- NULL
    if (length(gamma) > 0L) gamma_id <- any(gamma < association$lower | gamma > association$upper)

    ### Out
    if (any(!is.finite(mu)) | any(!is.finite(sigma)) | any(!is.finite(lambda[lambda_id])) |
        any(mu < 0) | any(sigma < 0) | any(nu[nu_id] < 0) | gamma_id | is.null(dec)) {

      -Inf

    }else {

      epsilon <- log_f <- matrix(NA, n, d)
      for(j in 1:d){

        epsilon[, j] <- qPSI(pmin(pmax(get(paste0("p", margins[j]), envir = parent.frame())(q = y[, j],
                                                                                            mu = mu[j],
                                                                                            sigma = sigma[j],
                                                                                            lambda = lambda[j],
                                                                                            nu = nu[j]), EPS), 1 - EPS))

        log_f[, j] <- get(paste0("d", margins[j]), envir = parent.frame())(x = y[, j],
                                                                           mu = mu[j],
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
  mu <- opt$par[par_id$mu]
  names(mu) <- colnames(y)

  sigma <- opt$par[par_id$sigma]
  names(sigma) <- colnames(y)

  lambda <- rep(NA, d)
  names(lambda) <- colnames(y)
  lambda[lambda_id] <- opt$par[par_id$lambda]

  nu <- rep(NA, d)
  names(nu) <- colnames(y)
  nu[nu_id] <- opt$par[par_id$nu]

  ### Out list
  out <- list(
    mu = mu,
    sigma = sigma,
    lambda = if (any(lambda_id)) lambda[lambda_id] else NULL,
    nu = if (any(nu_id)) nu[nu_id] else NULL,
    gamma = if (association$npar > 0L) opt$par[par_id$gamma] else NULL,
    margins = margins,
    association = association$name,
    copula = copula,
    eta = eta,
    logLik = opt$value,
    vcov = if (hessian) vcov else NULL,
    optim_pars = opt,
    y = y,
    nobs = n,
    d = d
  )

  ## Further model information ----
  out$call <- cl

  class(out) <- "bcsnsm"
  out

}
