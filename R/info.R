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
#' par.drm <- list(
#'   a = c(1.1, 1.2, 0.9, 1.8, 1.4),
#'   b = c(0.1, -1.6, -0.2, 1.0, 1.2),
#'   g = rep(0.2, 5)
#' )
#'
#' # create a list containing the polytomous item parameters
#' par.prm <- list(
#'   a = c(1.4, 0.6),
#'   d = list(c(-1.9, 0.0, 1.2), c(0.4, -1.1, 1.5, 0.2))
#' )
#'
#' # create a numeric vector of score categories for the items
#' cats <- c(2, 4, 2, 2, 5, 2, 2)
#'
#' # create a character vector of IRT models for the items
#' model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")
#'
#' # create an item metadata set
#' test <- shape_df(
#'   par.drm = par.drm, par.prm = par.prm,
#'   cats = cats, model = model
#' ) # create a data frame
#'
#' # set theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # compute item and test information values given the theta values
#' info(x = test, theta = theta, D = 1, tif = TRUE)
#'
#'
#' ## example 2.
#' ## using a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt",
#'   package = "irtQ"
#' )
#'
#' # read item parameters and transform them to item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # compute item and test information values given the theta values
#' info(x = test_flex, theta = theta, D = 1, tif = TRUE)
#'
#' @export
info <- function(x, ...) UseMethod("info")

#' @describeIn info Default method to compute item and test information functions for
#' a data frame \code{x} containing the item metadata.
#' @export
#' @importFrom Rfast colsums
info.default <- function(x, theta, D = 1, tif = TRUE, ...) {
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
  if (!is.null(idx.drm)) {
    # compute the item information
    iif_drm <-
      info_drm(
        theta = theta, a = elm_item$pars[idx.drm, 1],
        b = elm_item$par[idx.drm, 2], g = elm_item$par[idx.drm, 3],
        D = D, one.theta = FALSE, grad = FALSE, info = TRUE
      )$II
  } else {
    iif_drm <- NULL
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # count the number of examinees
    nstd <- length(theta)

    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])
    iif_prm <- NULL
    idx.prm <- c()
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]

      # reorder the index of the PRMs
      idx.prm <- c(idx.prm, elm_item$item[[mod]])

      # compute the probabilities of endorsing each score category
      iif_all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          grad = FALSE, info = TRUE
        )$II
      iif_all <- matrix(iif_all, ncol = nstd, byrow = FALSE)
      iif_prm <- rbind(iif_prm, iif_all)
    }
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
  if (tif) {
    tif.vec <- Rfast::colsums(iif)
  } else {
    tif.vec <- NULL
  }

  # return the results
  rst <- list(iif = iif, tif = tif.vec, theta = theta)
  class(rst) <- c("info")
  rst
}


#' @describeIn info An object created by the function \code{\link{est_item}}.
#' @export
#' @importFrom Rfast colsums
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
  if (!is.null(idx.drm)) {
    # compute the item information
    iif_drm <-
      info_drm(
        theta = theta, a = elm_item$pars[idx.drm, 1],
        b = elm_item$par[idx.drm, 2], g = elm_item$par[idx.drm, 3],
        D = D, one.theta = FALSE, grad = FALSE, info = TRUE
      )$II
  } else {
    iif_drm <- NULL
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # count the number of examinees
    nstd <- length(theta)

    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])
    iif_prm <- NULL
    idx.prm <- c()
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]

      # reorder the index of the PRMs
      idx.prm <- c(idx.prm, elm_item$item[[mod]])

      # compute the probabilities of endorsing each score category
      iif_all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          grad = FALSE, info = TRUE
        )$II
      iif_all <- matrix(iif_all, ncol = nstd, byrow = FALSE)
      iif_prm <- rbind(iif_prm, iif_all)
    }
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
  if (tif) {
    tif.vec <- Rfast::colsums(iif)
  } else {
    tif.vec <- NULL
  }

  # return the results
  rst <- list(iif = iif, tif = tif.vec, theta = theta)
  class(rst) <- c("info")
  rst
}

#' @describeIn info An object created by the function \code{\link{est_irt}}.
#' @export
#' @importFrom Rfast rowsums
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
  if (!is.null(idx.drm)) {
    # compute the item information
    iif_drm <-
      info_drm(
        theta = theta, a = elm_item$pars[idx.drm, 1],
        b = elm_item$par[idx.drm, 2], g = elm_item$par[idx.drm, 3],
        D = D, one.theta = FALSE, grad = FALSE, info = TRUE
      )$II
  } else {
    iif_drm <- NULL
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # count the number of examinees
    nstd <- length(theta)

    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])
    iif_prm <- NULL
    idx.prm <- c()
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]

      # reorder the index of the PRMs
      idx.prm <- c(idx.prm, elm_item$item[[mod]])

      # compute the probabilities of endorsing each score category
      iif_all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          grad = FALSE, info = TRUE
        )$II
      iif_all <- matrix(iif_all, ncol = nstd, byrow = FALSE)
      iif_prm <- rbind(iif_prm, iif_all)
    }
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
  if (tif) {
    tif.vec <- Rfast::colsums(iif)
  } else {
    tif.vec <- NULL
  }

  # return the results
  rst <- list(iif = iif, tif = tif.vec, theta = theta)
  class(rst) <- c("info")
  rst
}

# a function to compute the expected fisher item information for DRM items
# Pi(), Ji(), and li() functions of catR (Magis & Barrada, 2017) package were referred
# to compute the first and second derivatives of P with respect to theta and
# Ji value
#' @importFrom Rfast Outer
info_drm <- function(theta, a, b, g, D = 1, one.theta = FALSE,
                     r_i, grad = FALSE, info = TRUE, ji = FALSE) {
  # calculate probability of correct answers
  Da <- D * a
  if (!one.theta) {
    z <- Da * Rfast::Outer(x = theta, y = b, oper = "-")
  } else {
    z <- Da * (theta - b)
  }
  P <- g + (1 - g) / (1 + exp(-z))

  # prevent that the probabilities are equal to 1L or g (or 0 in case of 1PLM and 2PLM)
  P[P > 9999999999e-10] <- 9999999999e-10
  lg.lessg <- P < (g + 1e-10)
  P[lg.lessg] <- P[lg.lessg] + 1e-10

  # compute the fisher information
  if (info) {
    # compute the first and second derivatives of P wts the theta
    # referred to Pi.R in the catR package (Magis and Barranda, 2017)
    expz <- exp(z)
    exp1g <- expz * (1 - g)
    dP <- Da * exp1g / (1 + expz)^2

    # compute the item information (Hambleton & Swaminathan, 1985, p.107)
    # II <- {((D * a)^2 * (1 - P)) / P} * ((P - g)/(1 - g))^2
    Q <- (1 - P)
    II <- dP^2 / (P * Q)
  } else {
    II <- NULL
  }

  # compute the gradient of the negative loglikelihood (S)
  if (grad) {
    S <- -Da * ((P - g) * (r_i - P)) / ((1 - g) * P)
  } else {
    S <- NULL
  }

  # JI should be computed only for WL method
  if (ji & info) {
    d2P <- Da^2 * expz * (1 - expz) * (1 - g) / (1 + expz)^3
    J <- dP * d2P / (P * Q)
    #   d3P <- Da^3 * expz * (1 - g) * (expz^2 - 4 * expz + 1)/(1 + expz)^4
    #   dJ <- (P * Q * (d2P^2 + dP * d3P) - dP^2 * d2P *(Q - P))/(P^2 * Q^2)
    #   dI <- 2 * dP * d2P / P - dP^3 / P^2
  } else {
    J <- NULL
  }

  # return
  list(II = II, S = S, P = P, J = J)
}


# a function to compute the expected fisher item information for PRM items
# Pi(), Ji(), and li() functions of catR (Magis & Barrada, 2017) package were referred
# to compute the first and second derivatives of P with respect to theta and
# Ji value
#' @importFrom Rfast rowsums Outer
info_prm <- function(theta, a, d, D = 1, pr.model,
                     r_i, grad = FALSE, info = TRUE, ji = FALSE) {
  # check the number of step parameters
  m <- ncol(d)

  # GRM
  if (pr.model == "GRM") {
    # convert the d matrix to a vector
    d <- c(t(d))

    # calculate the probabilities greater than equal to each threshold
    Da <- D * a
    theta_d <- Rfast::Outer(x = theta, y = d, oper = "-")
    z <- Da * matrix(theta_d, ncol = m, byrow = TRUE)
    # z <- matrix(theta_d, ncol = m, byrow = TRUE)
    allP <- 1 / (1 + exp(-z))
    allP[is.na(allP)] <- 0 # insert 0 to NAs

    # calculate the probability for endorsing each score category (P)
    P <- cbind(1, allP) - cbind(allP, 0)
    P[P > 9999999999e-10] <- 9999999999e-10
    P[P < 1e-10] <- 1e-10

    # compute the fisher information
    if (info) {
      # calculate the first and second derivative of Ps wts the theta
      allQ <- 1 - allP[, , drop = FALSE]
      deriv_Pstth <- (Da * allP) * allQ
      dP <- cbind(0, deriv_Pstth) - cbind(deriv_Pstth, 0)
      w1 <- (1 - 2 * allP) * deriv_Pstth
      d2P <- Da * (cbind(0, w1) - cbind(w1, 0))
      # w2 <- Da * w1 * (1 - 2 * allP) - 2 * deriv_Pstth^2
      # d3P <- Da * (cbind(0, w2) - cbind(w2, 0))

      # compute the information for all score categories (IP)
      IP <- (dP^2 / P) - d2P

      # weight sum of all score category information
      # which is the item information (II)
      II <- Rfast::rowsums(IP)
    } else {
      II <- NULL
    }
  }

  # GPCM
  if (pr.model == "GPCM") {
    # include zero for the step parameter of the first category
    # and convert the d matrix to a vector
    d <- c(t(cbind(0, d)))

    # calculate the probability for endorsing each score category (P)
    Da <- D * a
    theta_d <- Rfast::Outer(x = theta, y = d, oper = "-")
    z <- matrix(theta_d, nrow = m + 1, byrow = FALSE)
    cumsum_z <- Da * t(Rfast::colCumSums(z))
    # cumsum_z[is.na(cumsum_z)] <- 0
    if (any(cumsum_z > 700, na.rm = TRUE)) {
      cumsum_z <- (cumsum_z / max(cumsum_z, na.rm = TRUE)) * 700
    }
    if (any(cumsum_z < -700, na.rm = TRUE)) {
      cumsum_z <- -(cumsum_z / min(cumsum_z, na.rm = TRUE)) * 700
    }
    numer <- exp(cumsum_z) # numerator
    numer[is.na(numer)] <- 0
    denom <- Rfast::rowsums(numer, na.rm = TRUE)
    P <- (numer / denom)
    P[P < 1e-10] <- 1e-10

    # compute the fisher information
    if (info) {
      # calculate the first and second derivative of Ps wts the theta
      denom2 <- denom^2
      denom4 <- denom2^2
      d1th_z <- (1:(m + 1))
      d1th_denom <- Da * as.numeric(tcrossprod(x = d1th_z, y = numer))
      d2th_denom <- Da^2 * as.numeric(tcrossprod(x = d1th_z^2, y = numer))
      d1th_z_den <- Da * outer(X = denom, Y = d1th_z, FUN = "*")
      dP <- (numer / denom2) * (d1th_z_den - d1th_denom)
      part1 <- d1th_z_den * denom * (d1th_z_den - d1th_denom)
      part2 <-
        denom2 * (Da * outer(X = d1th_denom, Y = d1th_z, FUN = "*") + d2th_denom) -
        2 * denom * d1th_denom^2
      d2P <- (numer / denom4) * (part1 - part2)

      # compute the information for all score categories (IP)
      IP <- (dP^2 / P) - d2P

      # weight sum of all score category information
      # which is the item information (II)
      II <- Rfast::rowsums(IP)
    } else {
      dP <- dP2 <- II <- NULL
    }
  }

  # compute the gradient of the negative loglikelihood (S)
  if (grad) {
    frac_rp <- r_i / P
    S <- -Rfast::rowsums(frac_rp * dP)
  } else {
    S <- NULL
  }

  # JI should be computed only for WL method
  # and S is also adjusted accordingly using the weight function of J/2I
  if (ji & info) {
    J <- Rfast::rowsums((dP * d2P) / P)
  } else {
    J <- NULL
  }

  # return the results
  list(II = II, S = S, P = P, J = J)
}

# This function computes the gradient of the negative loglikelihood and
# a negative expectation of second derivative of log-likelihood (fisher information)
info_score <- function(theta, elm_item, freq.cat, idx.drm, idx.prm,
                       method = c("ML", "MAP", "MLF", "WL"), D = 1,
                       norm.prior = c(0, 1), grad = TRUE, ji = FALSE) {
  # For DRM items
  if (!is.null(idx.drm)) {
    # extract the response and item parameters
    r_i <- freq.cat[idx.drm, 2]
    a <- elm_item$pars[idx.drm, 1]
    b <- elm_item$pars[idx.drm, 2]
    g <- elm_item$pars[idx.drm, 3]

    # compute the gradient and fisher information
    gi_drm <-
      info_drm(
        theta = theta, a = a, b = b, g = g,
        D = D, one.theta = TRUE, r_i = r_i, grad = grad,
        info = TRUE, ji = ji
      )
    grad_drm <- gi_drm$S
    finfo_drm <- gi_drm$II
    ji_drm <- gi_drm$J
  } else {
    # assign 0 to the gradient and fisher information
    grad_drm <- 0
    finfo_drm <- 0
    ji_drm <- 0
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    pr.mod <- unique(elm_item$model[idx.prm])
    grad_prm <- c()
    finfo_prm <- c()
    ji_prm <- c()
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      r_i <- freq.cat[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]

      # compute the gradient and fisher information
      gi_prm <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          r_i = r_i, grad = grad, info = TRUE, ji = ji
        )
      grad_prm <- c(grad_prm, gi_prm$S)
      finfo_prm <- c(finfo_prm, gi_prm$II)
      ji_prm <- c(ji_prm, gi_prm$J)
    }
  } else {
    # assign 0 to the gradient and fisher information
    grad_prm <- 0
    finfo_prm <- 0
    ji_prm <- 0
  }

  # sum of the gradients and the fisher information for DRM and PRM items
  grad_sum <- sum(grad_drm, grad_prm)
  finfo_sum <- sum(finfo_drm, finfo_prm)
  if (ji) {
    ji_sum <- sum(ji_drm, ji_prm)
    grad_sum <- grad_sum - (ji_sum / (2 * finfo_sum))
  }

  # extract the fisher information when MAP method is used
  if (method == "MAP") {
    # compute a gradient and hessian of prior distribution
    rst.prior <-
      logprior_deriv(
        val = theta, is.aprior = FALSE, D = NULL, dist = "norm",
        par.1 = norm.prior[1], par.2 = norm.prior[2]
      )

    # extract the hessian and add it
    finfo.prior <- attributes(rst.prior)$hessian
    finfo_sum <- sum(finfo_sum, finfo.prior)

    # add the gradient and add it
    if (grad) {
      grad.prior <- attributes(rst.prior)$gradient
      grad_sum <- sum(grad_sum, grad.prior)
    }
  }

  # return results
  list(finfo = finfo_sum, grad = grad_sum)
}
