#' Item and Test Information Function
#'
#' This function computes item and test information functions (Hambleton et
#' al., 1991) for a given set of theta values.
#'
#' @param x A data frame containing item metadata (e.g., item parameters,
#'   number of categories, IRT model types, etc.); or an object of class
#'   `est_irt` obtained from [irtQ::est_irt()], or `est_item` from
#'   [irtQ::est_item()].
#'
#'   See [irtQ::est_irt()] or [irtQ::simdat()] for more details about the item
#'   metadata. This data frame can be easily created using the
#'   [irtQ::shape_df()] function.
#' @param theta A numeric vector of theta values at which item and test
#'   information are computed.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.7 is commonly
#'   used for this purpose. Default is 1.
#' @param tif Logical. If `TRUE`, the test information function
#'   is computed. Default is `TRUE`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' This function calculates the amount of statistical information provided by
#' each item (item information function, IIF) and the total test (test
#' information function, TIF) across a range of ability (theta) values.
#' Higher information values at a particular theta level indicate greater
#' measurement precision at that ability level.
#'
#' The input `x` must follow a specific data frame format if not already
#' an `est_irt` or `est_item` object. The structure of this data frame
#' is explained in the documentation of [irtQ::est_irt()] and [irtQ::simdat()].
#' Items of different models (e.g., 3PLM, GPCM) can be combined in a single test.
#'
#' The information is computed for each item appropriately and aggregated for
#' the TIF if `tif = TRUE`. The TIF is often used to assess where the test
#' provides the most precision, and is critical when designing adaptive tests
#' or evaluating test coverage across the ability continuum.
#'
#' The returned object is a list of class `"info"`, which contains the item
#' information matrix and the test information vector. The `plot()` method for
#' `info` objects can be used to visualize the IIFs and TIF
#' (see [irtQ::plot.info()]).
#'
#' @return This function returns an object of class `info`, which is a list
#' containing the following components:
#'
#'   \item{iif}{A matrix of item information values. Each row corresponds to an
#'   item, and each column represents the information value computed at a given
#'   theta point. The row names are the item IDs, and the column names indicate
#'   the theta points (e.g., `"theta.1"`, `"theta.2"`, ...).}
#'
#'   \item{tif}{A numeric vector containing the test information values at each
#'   theta value, computed as the sum of item information values across all
#'   items. This component is included only when `tif = TRUE`.}
#'
#'   \item{theta}{A numeric vector of theta values at which the item and test
#'   information functions are evaluated. This matches the user-supplied
#'   `theta` argument.}
#'
#' The returned object is of class `info` and can be visualized using
#' the function [irtQ::plot.info()]. This output structure is consistent across
#' input types (`data.frame`, `est_item`, `est_irt`), and facilitates
#' downstream plotting, comparison, or export of information function values.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Hambleton, R. K., & Swaminathan, H. (1985) *Item response theory:
#'   Principles and applications*. Boston, MA: Kluwer.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991) *Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#' @seealso [irtQ::plot.info()], [irtQ::shape_df()], [irtQ::est_irt()],
#' [irtQ::est_item()]
#'
#' @examples
#' ## Example 1.
#' ## Using the function `shape_df()` to create a data frame of test metadata
#'
#' # Create a list of dichotomous item parameters
#' par.drm <- list(
#'   a = c(1.1, 1.2, 0.9, 1.8, 1.4),
#'   b = c(0.1, -1.6, -0.2, 1.0, 1.2),
#'   g = rep(0.2, 5)
#' )
#'
#' # Create a list of polytomous item parameters
#' par.prm <- list(
#'   a = c(1.4, 0.6),
#'   d = list(c(-1.9, 0.0, 1.2), c(0.4, -1.1, 1.5, 0.2))
#' )
#'
#' # Create a numeric vector for the number of score categories for each item
#' cats <- c(2, 4, 2, 2, 5, 2, 2)
#'
#' # Create a character vector of IRT model types for each item
#' model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")
#'
#' # Create item metadata using `shape_df()`
#' test <- shape_df(
#'   par.drm = par.drm, par.prm = par.prm,
#'   cats = cats, model = model
#' ) # create a data frame
#'
#' # Define a sequence of theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # Compute item and test information values based on theta
#' info(x = test, theta = theta, D = 1, tif = TRUE)
#'
#'
#' ## Example 2.
#' ## Using a "-prm.txt" file exported from flexMIRT
#'
#' # Import a sample "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt",
#'   package = "irtQ"
#' )
#'
#' # Read item parameters and convert them into item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Define a sequence of theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # Compute item and test information values based on theta
#' info(x = test_flex, theta = theta, D = 1, tif = TRUE)
#'
#' @export
info <- function(x, ...) UseMethod("info")

#' @describeIn info Default method to compute item and test information functions for
#' a data frame `x` containing the item metadata.
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


#' @describeIn info An object created by the function [irtQ::est_item()].
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

#' @describeIn info An object created by the function [irtQ::est_irt()].
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
