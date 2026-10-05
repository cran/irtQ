#' Simulated Mixed-Item Format CAT Data
#'
#' A simulated dataset for computerized adaptive testing (CAT), containing an
#' item pool, sparse response data, and examinee ability estimates. The item
#' pool includes both dichotomous and polytomous items.
#'
#' @usage simCAT_MX
#'
#' @format A list of length three:
#' \describe{
#'   \item{item.prm}{A data frame in item metadata format consisting of 200
#'   dichotomous items and 30 polytomous items.
#'   - Dichotomous items: Calibrated using the IRT 2PL model.
#'   - Polytomous items: Calibrated using the Generalized Partial Credit Model
#'   (GPCM), with four score categories (0, 1, 2, 3).}
#'
#'   \item{res.dat}{A matrix of item responses from 30,000 examinees (rows) to 230
#'   items (columns). `NA` marks an item that was not administered to the examinee.
#'   The columns are named `Item.dc.1` to `Item.dc.200` and `Item.py.1` to
#'   `Item.py.30`, not by the `id` column of `item.prm` (`V1` to `V230`); the
#'   columns follow the row order of `item.prm`.}
#'   \item{score}{A numeric vector of ability estimates for the 30,000
#'   examinees.}
#' }
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' # structure of the data
#' str(simCAT_MX, max.level = 1)
#'
#' \donttest{
#' # item fit of three dichotomous items and two polytomous items
#' loc <- c(1:3, 201:202)
#' x <- simCAT_MX$item.prm[loc, ]
#' data <- simCAT_MX$res.dat[, loc]
#' irtfit(
#'   x = x, score = simCAT_MX$score, data = data, group.method = "equal.freq",
#'   n.width = 10, loc.theta = "average", range.score = c(-4, 4), D = 1
#' )
#' }
#'
"simCAT_MX"
