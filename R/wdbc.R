#' Wisconsin Diagnostic Breast Cancer
#'
#' A subset of the the data originated from the multi-disciplinary study conducted by
#'     Mangasarian et al. (1995), focusing on the diagnosis and treatment of breast cancer. The
#'     dataset comprises measurements of the cell nuclei extracted from digitized images of fine
#'     needle aspirates of breast masses. For each patient included in the study, the diagnosis
#'     indicates whether the detected tumors were malignant or benign.
#'
#' @format A data frame with 556 rows and 6 columns:
#' \describe{
#'   \item{texture}{average variation in gray levels within boundaries.}
#'   \item{area}{average nucleus area of the cells.}
#'   \item{smoothness}{average local variation of radial segments.}
#'   \item{compactness}{average cell nucleus compactness, calculated as perimeter^2 / area - 1}
#'   \item{concavity}{average severity of concave portions of the contour}
#'   \item{diagnosis}{cancer diagnosis: \code{M = malignant, B = benign}.}
#' }
#'
#' @source The Breast Cancer Wisconsin (Diagnostic) Data Set
#'     UCI Machine Learning Repository \url{https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic}.
#'
#' @references
#'
#' Mangasarian, O. L., Street, W. N., and Wolberg, W. H. (1995). Breast cancer diagnosis and
#'     prognosis via linear programming. \emph{Operations Research}, \bold{43}, 570-577.
#'
#' @examples
#' mvplot(wdbc[, 1:5])
"wdbc"
