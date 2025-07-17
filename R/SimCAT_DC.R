#' Simulated Single-Item Format CAT Data
#'
#' A simulated dataset containing an item pool, sparse response data, and
#' examinee ability estimates, designed for single-item computerized adaptive
#' testing (CAT).
#'
#' @usage simCAT_DC
#'
#' @format A list of length three:
#' \describe{
#'   \item{item_pool}{A data frame in item metadata format containing 100 dichotomous items.
#'   - Items 1–90: Generated and calibrated under the IRT 2PL model.
#'   - Items 91–100: Generated under the IRT 3PL model but calibrated using the 2PL model.}
#'
#'   \item{response_data}{A sparse matrix of item responses from 10,000 examinees.}
#'   \item{theta_estimates}{A numeric vector of ability estimates for the 10,000 examinees.}
#' }
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
"simCAT_DC"
