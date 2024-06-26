#' Vitamin Intake
#'
#' A subset of a study conducted by Ferreira et al. (2014) on the nutrient intake patterns of
#'     older adults in Botucatu, São Paulo, Brazil. The data consist of the vitamin intake
#'     based on the first 24h dietary recall interview for older women, where each sampling
#'     unit provides a vector with dimension 8 with the intake of vitamins D, K, C, B1, B2,
#'     B3, B6, and B12.
#'
#' @format 
#' 
#' A data frame with 228 rows and 8 columns:
#' \describe{
#'   \item{vitd}{Intake of vitamin D, in \emph{mcg}.}
#'   \item{vitk}{Intake of vitamin K, in \emph{mcg}.}
#'   \item{vitc}{Intake of vitamin C, in \emph{mg}.}
#'   \item{vitb1}{Intake of vitamin B1, in \emph{mg}.}
#'   \item{vitb2}{Intake of vitamin B2, in \emph{mg}.}
#'   \item{vitb3}{Intake of vitamin B3, in \emph{mg}.}
#'   \item{vitb6}{Intake of vitamin B6, in \emph{mg}.}
#'   \item{vitb12}{Intake of vitamin B12, in \emph{mcg}.}
#' }
#'
#' @references
#'
#' Ferreira, P. M., Papini, S. J., and Corrente, J. E. (2014). Diversity of eating patterns in older
#'     adults: A new scenario? \emph{Revista de Nutrição}, \bold{27}, 67–79.
#'
#'
#' @examples
#' mvplot(vitamins)
"vitamins"
