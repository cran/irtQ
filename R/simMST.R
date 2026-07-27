#' Simulated 1-3-3 MST Panel Data
#'
#' A simulated Multistage-Adaptive Test (MST) dataset representing a 1-3-3 panel
#' structure, assembled using module-level target test information functions
#' (TIFs) and design constraints similar to those used in the simulation study
#' by Lim et al. (2021).
#'
#' @usage simMST
#'
#' @format A list containing five internal objects:
#' \describe{
#'   \item{item_bank}{A data frame of item metadata for the 56 items assigned
#'   to the panel (one row per item; columns follow the \pkg{irtQ} item
#'   metadata format: \code{id}, \code{cats}, \code{model}, \code{par.1},
#'   \code{par.2}, \code{par.3}).}
#'
#'   \item{module}{A binary matrix that maps items in the item bank to MST
#'   modules. This structure specifies the item-to-module assignment,
#'   similar to the `modules` argument in the `randomMST()` function from
#'   the \pkg{mstR} package (Magis et al., 2017).}
#'
#'   \item{route_map}{A binary square matrix that defines the MST transition
#'   structure, showing module pathways across stages.
#'   This corresponds to the `transMatrix` argument in the `randomMST()`
#'   function in the \pkg{mstR} package.}
#'
#'   \item{cut_score}{A list of numeric vectors specifying the routing cut
#'   scores between MST stages, computed via \code{\link{find_cut}()} based
#'   on TIF-crossing points of the assembled modules.}
#'
#'   \item{theta}{A numeric vector of ability (theta) values used to evaluate
#'    the panel's measurement precision across the latent trait continuum.}
#' }
#'
#'   This 1-3-3 MST panel includes 7 modules across 3 stages, drawn from
#'   56 unique items (8 per module) with no item shared across modules. Each
#'   module contains 8 dichotomously scored items calibrated under the IRT
#'   3-parameter logistic (3PL) model.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Magis, D., Yan, D., & von Davier, A. A. (2017). *Computerized
#' adaptive and multistage testing with R: Using packages catR and mstR*.
#' Springer.
#'
#' Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
#' approach to evaluate the performance of MST. *Journal of Educational
#' Measurement, 58*(2), 154-178.
#'
"simMST"
