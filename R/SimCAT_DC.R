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
#'   \item{item.prm}{A data frame in item metadata format containing 100 dichotomous items.
#'   - Items 1-90: Generated and calibrated under the IRT 2PL model.
#'   - Items 91-100: Generated under the IRT 3PL model but calibrated using the 2PL model.}
#'
#'   \item{res.dat}{A matrix of item responses from 10,000 examinees (rows) to 100
#'   items (columns). `NA` marks an item that was not administered to the examinee.
#'   The columns have no names; they follow the row order of `item.prm`, whose
#'   `id` values are `V1` to `V100`.}
#'   \item{score}{A numeric vector of ability estimates for the 10,000 examinees.}
#' }
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' # structure of the data
#' str(simCAT_DC, max.level = 1)
#'
#' \donttest{
#' # item fit of the first five items, using the ability estimates as scores
#' x <- simCAT_DC$item.prm[1:5, ]
#' data <- simCAT_DC$res.dat[, 1:5]
#' irtfit(
#'   x = x, score = simCAT_DC$score, data = data, group.method = "equal.freq",
#'   n.width = 10, loc.theta = "average", range.score = c(-4, 4), D = 1
#' )
#' }
#'
"simCAT_DC"
