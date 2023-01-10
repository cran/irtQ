#' Item and Test Information Function
#'
#' @description This function computes both item and test information functions (Hambleton et al., 1991) given a set of theta values.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object of class \code{\link{est_item}}
#' obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}.
#' The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}. See below for details.
#' @param theta A vector of theta values where item and test information values are computed.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param tif A logical value. If TRUE, the test information function is computed. Default is TRUE.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details A specific form of a data frame should be used for the argument \code{x}. The first column should have item IDs,
#' the second column should contain unique score category numbers of the items, and the third column should include IRT models being fit to the items.
#' The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous item data, and "GRM" and "GPCM" for polytomous item data.
#' Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
#' response model and (generalized) partial credit model, respectively. The next columns should include the item parameters of the fitted IRT models.
#' For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
#' item guessing parameters, respectively. When "1PLM" and "2PLM" are specified in the third column, NAs should be inserted in the sixth column
#' for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be included in the
#' fourth column and the item difficulty (or threshold) parameters of category boundaries should be contained from the fifth to the last columns.
#' When the number of unique score categories differs between items, the empty cells of item parameters should be filled with NAs.
#' In the \pkg{irtQ} package, the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as
#' the item location (or overall difficulty) parameter subtracted by the threshold parameter for unique score categories of the item.
#' Note that when an GPCM item has \emph{K} unique score categories, \emph{K-1} item difficulty parameters are necessary because
#' the item difficulty parameter for the first category boundary is always 0. For example, if an GPCM item has five score categories,
#' four item difficulty parameters should be specified. An example of a data frame with a single-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
#'   ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
#'   ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
#' }
#' And an example of a data frame for a mixed-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
#'   ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
#'   ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
#'   ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
#'   ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
#'   ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
#' }
#' See \code{IRT Models} section in the page of \code{\link{irtQ-package}} for more details about the IRT models used in the \pkg{irtQ} package.
#' An easier way to create a data frame for the argument \code{x} is by using the function \code{\link{shape_df}}.
#'
#' @return This function returns an object of class \code{\link{info}}. This object contains item and test information values
#' given the specified theta values.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Hambleton, R. K., & Swaminathan, H. (1985) \emph{Item response theory: Principles and applications}.
#' Boston, MA: Kluwer.
#'
#' Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991) \emph{Fundamentals of item response theory}.
#' Newbury Park, CA: Sage.
#'
#' @seealso \code{\link{plot.info}}, \code{\link{shape_df}}, \code{\link{est_item}}
#'
#' @examples
#' ## example 1.
#' ## using the function "shape_df" to create a data frame of test metadata
#' # create a list containing the dichotomous item parameters
#' par.drm <- list(a=c(1.1, 1.2, 0.9, 1.8, 1.4),
#'                b=c(0.1, -1.6, -0.2, 1.0, 1.2),
#'                g=rep(0.2, 5))
#'
#' # create a list containing the polytomous item parameters
#' par.prm <- list(a=c(1.4, 0.6),
#'                d=list(c(0.0, -1.9, 1.2), c(0.4, -1.1, 1.5, 0.2)))
#'
#' # create a numeric vector of score categories for the items
#' cats <- c(2, 4, 2, 2, 5, 2, 2)
#'
#' # create a character vector of IRT models for the items
#' model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")
#'
#' # create an item metadata set
#' test <- shape_df(par.drm=par.drm, par.prm=par.prm,
#'                  cats=cats, model=model) # create a data frame
#'
#' # set theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # compute item and test information values given the theta values
#' info(x=test, theta=theta, D=1, tif=TRUE)
#'
#'
#' ## example 2.
#' ## using a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt",
#'                         package = "irtQ")
#'
#' # read item parameters and transform them to item metadata
#' test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # compute item and test information values given the theta values
#' info(x=test_flex, theta=theta, D=1, tif=TRUE)
#'
#' @export
#'
#' @export
info <- function(x, ...) UseMethod("info")

#' @describeIn info Default method to compute item and test information functions for
#' a data frame \code{x} containing the item metadata.
#' @export
#' @importFrom Rfast rowsums
#' @importFrom stats na.exclude
info.default <- function(x, theta, D=1, tif = TRUE, ...) {

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # extract item id
  id <- elm_item$id

  # For DRM items
  if(!is.null(idx.drm)) {

    # compute the item information
    iif_drm <- info.drm(theta=theta, a=elm_item$pars[idx.drm, 1],
                        b=elm_item$par[idx.drm, 2], g=elm_item$par[idx.drm, 3], D=D)

  } else {
    iif_drm <- NULL
  }

  # For PRM items
  if(!is.null(idx.prm)) {

    # count the number of PRM items
    n.prm <- length(idx.prm)
    iif_prm <- vector('list', n.prm)

    # compute the item information
    for(k in 1:n.prm) {
      par.tmp <- stats::na.exclude(elm_item$par[idx.prm[k], ])
      iif_prm[[k]] <- info.prm(theta=theta, a=par.tmp[1], d=par.tmp[-1],
                               D=D, pr.model=elm_item$model[idx.prm[k]])
    }
    iif_prm <- do.call(what='rbind', iif_prm)

  } else {
    iif_prm <- NULL
  }

  # create an item information matrix for all items
  iif <- rbind(iif_drm, iif_prm)

  # re-order the item information matrix along with the original order of items
  loc.item <- c(idx.drm, idx.prm)
  iif <- iif[order(loc.item), , drop = FALSE]
  rownames(iif) <- id
  colnames(iif) <- paste0("theta.", 1:length(theta))

  # compute the test information
  if(tif) {
    tif.vec <- Rfast::colsums(iif)
  } else {
    tif.vec <- NULL
  }

  # return the results
  rst <- list(iif=iif, tif=tif.vec, theta=theta)
  class(rst) <- c("info")
  rst

}

#' @describeIn info An object created by the function \code{\link{est_item}}.
#' @export
#' @importFrom Rfast rowsums
#' @importFrom stats na.exclude
info.est_item <- function(x, theta, tif = TRUE, ...) {

  # extract information from an object
  D <- x$scale.D
  x <- x$par.est

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # extract item id
  id <- elm_item$id

  # For DRM items
  if(!is.null(idx.drm)) {

    # compute the item information
    iif_drm <- info.drm(theta=theta, a=elm_item$pars[idx.drm, 1],
                        b=elm_item$par[idx.drm, 2], g=elm_item$par[idx.drm, 3], D=D)

  } else {
    iif_drm <- NULL
  }

  # For PRM items
  if(!is.null(idx.prm)) {

    # count the number of PRM items
    n.prm <- length(idx.prm)
    iif_prm <- vector('list', n.prm)

    # compute the item information
    for(k in 1:n.prm) {
      par.tmp <- stats::na.exclude(elm_item$par[idx.prm[k], ])
      iif_prm[[k]] <- info.prm(theta=theta, a=par.tmp[1], d=par.tmp[-1],
                               D=D, pr.model=elm_item$model[idx.prm[k]])
    }
    iif_prm <- do.call(what='rbind', iif_prm)

  } else {
    iif_prm <- NULL
  }

  # create an item information matrix for all items
  iif <- rbind(iif_drm, iif_prm)

  # re-order the item information matrix along with the original order of items
  loc.item <- c(idx.drm, idx.prm)
  iif <- iif[order(loc.item), , drop = FALSE]
  rownames(iif) <- id
  colnames(iif) <- paste0("theta.", 1:length(theta))

  # compute the test information
  if(tif) {
    tif.vec <- Rfast::colsums(iif)
  } else {
    tif.vec <- NULL
  }

  # return the results
  rst <- list(iif=iif, tif=tif.vec, theta=theta)
  class(rst) <- c("info")
  rst

}

#' @describeIn info An object created by the function \code{\link{est_irt}}.
#' @export
#' @importFrom Rfast rowsums
#' @importFrom stats na.exclude
info.est_irt <- function(x, theta, tif = TRUE, ...) {

  # extract information from an object
  D <- x$scale.D
  x <- x$par.est

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # extract item id
  id <- elm_item$id

  # For DRM items
  if(!is.null(idx.drm)) {

    # compute the item information
    iif_drm <- info.drm(theta=theta, a=elm_item$pars[idx.drm, 1],
                        b=elm_item$par[idx.drm, 2], g=elm_item$par[idx.drm, 3], D=D)

  } else {
    iif_drm <- NULL
  }

  # For PRM items
  if(!is.null(idx.prm)) {

    # count the number of PRM items
    n.prm <- length(idx.prm)
    iif_prm <- vector('list', n.prm)

    # compute the item information
    for(k in 1:n.prm) {
      par.tmp <- stats::na.exclude(elm_item$par[idx.prm[k], ])
      iif_prm[[k]] <- info.prm(theta=theta, a=par.tmp[1], d=par.tmp[-1],
                               D=D, pr.model=elm_item$model[idx.prm[k]])
    }
    iif_prm <- do.call(what='rbind', iif_prm)

  } else {
    iif_prm <- NULL
  }

  # create an item information matrix for all items
  iif <- rbind(iif_drm, iif_prm)

  # re-order the item information matrix along with the original order of items
  loc.item <- c(idx.drm, idx.prm)
  iif <- iif[order(loc.item), , drop = FALSE]
  rownames(iif) <- id
  colnames(iif) <- paste0("theta.", 1:length(theta))

  # compute the test information
  if(tif) {
    tif.vec <- Rfast::colsums(iif)
  } else {
    tif.vec <- NULL
  }

  # return the results
  rst <- list(iif=iif, tif=tif.vec, theta=theta)
  class(rst) <- c("info")
  rst

}


# a function to compute the expected fisher item information for DRM items
#' @importFrom Rfast Outer
info.drm <- function(theta, a, b, g, D=1, one.theta = FALSE) {

  # calculate probability of correct answers
  if(!one.theta) {
    z <- (D * a) * Rfast::Outer(x = theta, y = b, oper = "-")
  } else {
    z <- (D * a) * (theta - b)
  }
  P <- g + (1 - g) / (1 + exp(-z))

  # prevent that the probabilities are equal to 1L or g (or 0 in case of 1PLM and 2PLM)
  P[P > 9999999999e-10] <- 9999999999e-10
  lg.lessg <- P < (g + 1e-10)
  P[lg.lessg] <- P[lg.lessg] + 1e-10

  # compute the item information (Hambleton & Swaminathan, 1985, p.107)
  II <- {((D * a)^2 * (1 - P)) / P} * ((P - g)/(1 - g))^2

  # return
  II

}

# a function to compute the expected fisher item information for PRM items
#' @importFrom Rfast rowsums
info.prm <- function(theta, a, d, D=1, pr.model) {

  if(pr.model == "GRM") {

    # count the number of d parameters
    m <- length(d)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta=theta, a=a, b=d, g=0, D=D)
    allPst[allPst > 9999999999e-10] <- 9999999999e-10
    allPst[allPst < 1e-10] <- 1e-10
    allQst <- 1 - allPst[, ,drop=FALSE]

    # calculate category probabilities
    P <- cbind(1, allPst) - cbind(allPst, 0)
    P[P > 9999999999e-10] <- 9999999999e-10
    P[P < 1e-10] <- 1e-10

    # compute the component values to obtain information
    Da <- D * a
    deriv_Pstth <- (Da * allPst) * allQst
    deriv_Pth <- cbind(0, deriv_Pstth) - cbind(deriv_Pstth, 0)
    w1 <- (1 - 2 * allPst) * deriv_Pstth
    deriv2_Pth <- Da * (cbind(0, w1) - cbind(w1, 0))

    # compute the information for all score categories
    ip <- (deriv_Pth^2 / P) - deriv2_Pth

  }

  if(pr.model == "GPCM") {

    # include zero for the step parameter of the first category
    d <- c(0, d)

    # check the number of step parameters
    m <- length(d) - 1

    # calculate category probabilities
    Da <- D * a
    theta_d <- (Rfast::Outer(x =  theta, y = d, oper = "-"))
    z <- Da * theta_d
    cumsum_z <- t(Rfast::colCumSums(z))
    if(any(cumsum_z > 700)) {cumsum_z <- (cumsum_z / max(cumsum_z)) * 700}
    if(any(cumsum_z < -700)) {cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700}
    numer <- exp(cumsum_z) # numerator
    denom <- Rfast::rowsums(numer)
    P <- (numer / denom)

    # compute the component values to obtain information
    denom2 <- denom^2
    d1th_z <- matrix(Da * (1:(m+1)), nrow=length(theta), ncol=(m+1), byrow=TRUE)
    d1th_denom <- Rfast::rowsums(numer * d1th_z)
    deriv_Pth <- (numer / denom2) * (d1th_z * denom - d1th_denom)
    denom4 <- denom2^2
    d2th_denom <- Rfast::rowsums(numer * d1th_z^2)
    d1th_z_den <- d1th_z * denom
    part1 <- d1th_z_den * denom * (d1th_z_den - d1th_denom)
    part2 <- denom2 * (d1th_z * d1th_denom + d2th_denom) - 2 * denom * d1th_denom^2
    deriv2_Pth <- (numer / denom4) * (part1 - part2)

    # compute the information for all score categories
    ip <- (deriv_Pth^2 / P) - deriv2_Pth

  }

  # sum of all score category information
  II <- Rfast::rowsums(ip)

  # return the results
  II

}


# This function computes a negative expectation of second derivative of log-likelihood (fisher information)
info_score <- function(theta, elm_item, idx.drm, idx.prm,
                       method=c("ML", "MAP", "MLF"), D=1, norm.prior=c(0, 1)) {

  # For DRM items
  if(!is.null(idx.drm)) {

    # compute the fisher information
    finfo_drm <-
      sum(info.drm(theta=theta, a=elm_item$pars[idx.drm, 1],
                   b=elm_item$par[idx.drm, 2], g=elm_item$par[idx.drm, 3],
                   D=D, one.theta=TRUE))

  } else {

    # assign 0 to the fisher information
    finfo_drm <- 0

  }

  # For PRM items
  if(!is.null(idx.prm)) {

    # count the number of the PRM items
    n.prm <- length(idx.prm)

    # compute the fisher information
    finfo_prm <- c()
    for(i in 1:n.prm) {
      par.tmp <- stats::na.exclude(elm_item$par[idx.prm[i], ])
      finfo_prm[i] <-
        info.prm(theta=theta, a=par.tmp[1], d=par.tmp[-1],
                 D=D, pr.model=elm_item$model[idx.prm[i]])
    }

    # sum of the information for all PRM items
    finfo_prm <- sum(finfo_prm)

  } else {

    # assign 0 to the fisher information
    finfo_prm <- 0

  }

  # sum of the fisher information for DRM and PRM items
  finfo <- sum(finfo_drm, finfo_prm)

  # extract the fisher information when MAP method is used
  if(method == "MAP") {

    # compute a hessian of prior distribution
    rst.prior <-
      logprior_deriv(val=theta, is.aprior=FALSE, D=NULL, dist="norm",
                     par.1=norm.prior[1], par.2=norm.prior[2])

    # extract the hessian
    finfo.prior <- attributes(rst.prior)$hessian

    # add the hessian
    finfo <- sum(finfo, finfo.prior)

  }

  # return results
  finfo

}

