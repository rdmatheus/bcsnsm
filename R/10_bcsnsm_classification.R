#' BCS-NSM Classification
#'
#' @param object an object of class \code{"bcsnsmreg"}.
#' @param test matrix or data frame of test set cases.
#' @param cl factor of true classifications of training set.
#'
#' @return Factor of classifications of test set.
#'
#' @export
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' \dontrun{
#' # Winscosin Breast Cancer Dataset
#' ?wdbc
#' 
#' # Training set index
#' set.seed(123)
#' id <- sample(1:nrow(wdbc), 0.7 * nrow(wdbc))
#' 
#' # BCNSM regression
#' fit <- bcsnsmreg(texture + area + smoothness + compactness + concavity ~ diagnosis,
#'                  data = wdbc[id, ], margins = c("lt", "lt", "lno", "lpe", "bct"))
#' 
#' ## Marginal quantile residuals
#' plot(fit_gaussian, "marginal", panel = c(2, 3))
#' 
#' ## Summary
#' summary(fit_gaussian)
#' 
#' ## Overall quantile residuals of the final model
#' plot(fit_gaussian)
#' 
#' ## The epsilon's transformations
#' plot(fit_gaussian, "epsilon")
#' 
#' ## Classification
#' pred <- bcsnsmclass(fit_gaussian, wdbc[-id, 1:5], wdbc$diagnosis[id])
#' table(pred = pred, true = wdbc$diagnosis[-id])
#' }
#'
bcsnsmclass <- function(object, test, cl){

  cl <- factor(cl)
  y <- as.matrix(test)
  d <- object$d
  pk <- prop.table(table(cl))

  margins <- object$margins
  lambda_id <- apply(matrix(margins, ncol = 1), 1, function(x) !grepl("Log-", as.bcs(x)$name))
  nu_id <- apply(matrix(margins, ncol = 1), 1, function(x) as.bcs(x)$extrap)

  if (object$association == "non-associative"){
    Gamma <- diag(d)
  }else{
    Gamma <- get(tolower(object$association), envir = parent.frame())(d)$Gamma(round(object$gamma, 4))
  }

  mu <- object$marginal.pars$mu
  sigma <- object$marginal.pars$sigma
  lambda <- rep(NA, d)
  lambda[lambda_id] <- object$marginal.pars$lambda
  nu <- rep(NA, d)
  nu[nu_id] <- object$marginal.pars$nu

  out <- matrix(NA, nrow = nrow(y), ncol = length(pk))
  i <- 1
  for(k in levels(cl)){

    out[, i] <- dbcsnsm(as.matrix(y), mu[which(cl == k), ][1, ], sigma, lambda, nu, Gamma,
                        copula = object$copula, eta = object$eta, margins = object$margins) * pk[i]

    i <- i + 1

  }

  aux <- function(x){
    which.max(x)
  }

  levels(cl)[apply(out, 1, which.max)]

}

