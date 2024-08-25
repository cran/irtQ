#' Simulated 1-3-3 MST panel data
#'
#' This simulated 1-3-3 MST panel data set was used in Lim et al.' (2020) simulation study.
#'
#' @usage simMST
#'
#' @format This data set includes a list of five internal objects:
#' \describe{
#'   \item{\code{item_bank}}{The item bank metadata, containing item parameters
#'   and other information.}
#'   \item{\code{module}}{A binary matrix that maps items from the item bank
#'   to modules within the MST panel. This parameter enables precise
#'   item-to-module assignments for MST configurations, analogous to
#'   the \code{modules} argument in the \code{randomMST} function
#'   of the \pkg{mstR} package (Magis et al., 2017).}
#'   \item{\code{route_map}}{A binary square matrix that defines the MST structure,
#'   illustrating transitions between modules and stages. This concept is inspired
#'   by the \code{transMatrix} argument in the \code{randomMST}
#'   function from the \pkg{mstR} package (Magis et al., 2017).}
#'   \item{\code{cut_score}}{A list defining cut scores for routing test takers
#'   through MST stages. Each list element is a vector of cut scores for advancing
#'   test takers to subsequent stage modules.}
#'   \item{\code{theta}}{A vector of ability levels (theta) at which the MST
#'   panel's performance is assessed, allowing for the evaluation of measurement
#'   precision across a continuum of ability levels.}
#' }
#'
#' This 1-3-3 MST panel consists of 7 modules in total with 3 stages. Each module
#' contains eight items calibrated with the IRT 3 parameter logistic model.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Magis, D., Yan, D., & Von Davier, A. A. (2017). \emph{Computerized adaptive
#' and multistage testing with R: Using packages catR and mstR}. Springer.
#'
#' Lim, H., Davey, T., & Wells, C. S. (2020). A recursion-based analytical
#' approach to evaluate the performance of MST. \emph{Journal of Educational
#' Measurement, 58}(2), 154-178.
#'
"simMST"
