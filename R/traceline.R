#' Compute Item/Test Characteristic Functions
#'
#' This function computes item category probabilities, item characteristic
#' curves (ICCs), and the test characteristic curve (TCC) for a given set of
#' theta values. The returned object can be used to visualize these functions
#' using [irtQ::plot.traceline()].
#'
#' @inheritParams info
#' @inheritParams est_score
#' @param theta A numeric vector of theta values at which item and test
#'   characteristic curves are computed.
#'
#' @details
#' This function computes the item and test characteristic functions commonly
#' used in IRT. For each item, the function computes the category response
#' probabilities across a specified set of theta values. These probabilities are
#' used to derive:
#'
#'  - The item characteristic curve (ICC), which represents the expected score
#'  of each item as a function of theta.
#'  - The test characteristic curve (TCC), which is the sum of expected item
#'  scores at each theta value.
#'
#' The output object can be visualized using the [irtQ::plot.traceline] to
#' inspect the relationship between ability levels (theta) and expected
#' item/test scores.
#'
#' If the input `x` is an object of class `est_item` or `est_irt`, the function
#' automatically extracts item parameter estimates and the scaling constant
#' `D` from the object. Otherwise, a properly formatted item metadata data frame
#' must be provided.
#'
#' @return This function returns an object of class `traceline`, which is a list
#' containing the following components:
#'
#'   \item{prob.cats}{A list of data frames containing the category response
#'   probabilities for each item across the specified theta values. Each data
#'   frame corresponds to an item, with rows representing theta values and
#'   columns representing response categories (e.g., `"resp.0"`, `"resp.1"`, ...).}
#'
#'   \item{icc}{A numeric matrix representing ICCs. Each column corresponds to
#'   an item, and each row represents the expected item score at a given theta value.
#'   The column names are the item IDs.}
#'
#'   \item{tcc}{A numeric vector representing the TCC, computed as the sum of
#'   expected item scores across all items at each theta value.}
#'
#'   \item{theta}{A numeric vector of theta values at which the item and test
#'   information functions are evaluated. This matches the user-supplied
#'   `theta` argument.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::plot.traceline()], [irtQ::est_irt()], [irtQ::est_item()]
#'
#' @examples
#' ## Example using a "-prm.txt" file exported from flexMIRT
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read the item parameters and convert them into item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Define a sequence of theta values
#' theta <- seq(-3, 3, 0.5)
#'
#' # Compute item category probabilities, ICCs,
#' # and the TCC for the given theta values
#' traceline(x = test_flex, theta, D = 1)
#'
#' @export
traceline <- function(x, ...) UseMethod("traceline")

#' @describeIn traceline Default method to compute the item category probabilities,
#' item characteristic function, and test characteristic function for a data frame
#' `x` containing the item metadata.
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
traceline.default <- function(x, theta, D = 1, ...) {
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several components
  elm_item <- breakdown(x)

  # set indices of the DRM and PRM items
  idx.all <- idxfinder(elm_item)
  idx.drm <- idx.all$idx.drm
  idx.prm <- idx.all$idx.prm

  # count the number of thetas
  n.theta <- length(theta)

  # count the total number of items
  n.item <- sum(length(idx.drm), length(idx.prm))

  # make the empty list and data frame to contain icc and tcc
  prob.cats <- vector("list", n.item)
  icc_df <- array(NA, c(n.theta, n.item))
  names(prob.cats) <- x$id
  colnames(icc_df) <- x$id

  # For DRM items
  if (!is.null(idx.drm)) {
    # compute the probabilities of correct answers
    p.mat <- drm(
      theta = theta, a = elm_item$pars[idx.drm, 1],
      b = elm_item$pars[idx.drm, 2], g = elm_item$pars[idx.drm, 3],
      D = D
    )
    p.vec <- c(p.mat)
    q.vec <- 1 - p.vec

    # split the probabilities into each item group
    prob.drm <-
      split.data.frame(
        cbind(resp.0 = q.vec, resp.1 = p.vec),
        rep(idx.drm, each = n.theta)
      )

    # fill the empty list
    prob.cats[idx.drm] <- prob.drm

    # insert icc
    icc_df[, idx.drm] <- p.mat
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])

    # check the maximum score category
    max.cats <- max(elm_item$cats) - 1
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]
      cat.tmp <- elm_item$cats[lg.prm]

      # reorder the index of the PRMs
      idx.tmp <- elm_item$item[[mod]]

      # compute the probabilities of endorsing each score category
      P.all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          grad = FALSE, info = FALSE
        )$P
      colnames(P.all) <- paste0("resp.", 0:max.cats)
      prob.cats[idx.tmp] <-
        split.data.frame(
          P.all,
          rep(1:length(a), n.theta)
        ) %>%
        purrr::map2(
          .y = cat.tmp,
          .f = ~ {
            .x[, 1:(.y)]
          }
        )

      # insert icc
      icc_df[, idx.tmp] <-
        matrix(P.all %*% c(0:max.cats), nrow = n.theta, byrow = TRUE)
    }
  }

  # compute tcc
  tcc.vec <- Rfast::rowsums(icc_df)

  # return results
  rst <- list(prob.cats = prob.cats, icc = icc_df, tcc = tcc.vec, theta = theta)
  class(rst) <- "traceline"
  rst
}



#' @describeIn traceline An object created by the function [irtQ::est_item()].
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
traceline.est_item <- function(x, theta, ...) {
  # extract information from an object
  D <- x$scale.D
  x <- x$par.est

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several components
  elm_item <- breakdown(x)

  # set indices of the DRM and PRM items
  idx.all <- idxfinder(elm_item)
  idx.drm <- idx.all$idx.drm
  idx.prm <- idx.all$idx.prm

  # count the number of thetas
  n.theta <- length(theta)

  # count the total number of items
  n.item <- sum(length(idx.drm), length(idx.prm))

  # make the empty list and data frame to contain icc and tcc
  prob.cats <- vector("list", n.item)
  icc_df <- array(NA, c(n.theta, n.item))
  names(prob.cats) <- x$id
  colnames(icc_df) <- x$id

  # For DRM items
  if (!is.null(idx.drm)) {
    # compute the probabilities of correct answers
    p.mat <- drm(
      theta = theta, a = elm_item$pars[idx.drm, 1],
      b = elm_item$pars[idx.drm, 2], g = elm_item$pars[idx.drm, 3],
      D = D
    )
    p.vec <- c(p.mat)
    q.vec <- 1 - p.vec

    # split the probabilities into each item group
    prob.drm <-
      split.data.frame(
        cbind(resp.0 = q.vec, resp.1 = p.vec),
        rep(idx.drm, each = n.theta)
      )

    # fill the empty list
    prob.cats[idx.drm] <- prob.drm

    # insert icc
    icc_df[, idx.drm] <- p.mat
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])

    # check the maximum score category
    max.cats <- max(elm_item$cats) - 1
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]
      cat.tmp <- elm_item$cats[lg.prm]

      # reorder the index of the PRMs
      idx.tmp <- elm_item$item[[mod]]

      # compute the probabilities of endorsing each score category
      P.all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          grad = FALSE, info = FALSE
        )$P
      colnames(P.all) <- paste0("resp.", 0:max.cats)
      prob.cats[idx.tmp] <-
        split.data.frame(
          P.all,
          rep(1:length(a), n.theta)
        ) %>%
        purrr::map2(
          .y = cat.tmp,
          .f = ~ {
            .x[, 1:(.y)]
          }
        )

      # insert icc
      icc_df[, idx.tmp] <-
        matrix(P.all %*% c(0:max.cats), nrow = n.theta, byrow = TRUE)
    }
  }

  # compute tcc
  tcc.vec <- Rfast::rowsums(icc_df)

  # return results
  rst <- list(prob.cats = prob.cats, icc = icc_df, tcc = tcc.vec, theta = theta)
  class(rst) <- "traceline"
  rst
}


#' @describeIn traceline An object created by the function [irtQ::est_irt()].
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
traceline.est_irt <- function(x, theta, ...) {
  # extract information from an object
  D <- x$scale.D
  x <- x$par.est

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several components
  elm_item <- breakdown(x)

  # set indices of the DRM and PRM items
  idx.all <- idxfinder(elm_item)
  idx.drm <- idx.all$idx.drm
  idx.prm <- idx.all$idx.prm

  # count the number of thetas
  n.theta <- length(theta)

  # count the total number of items
  n.item <- sum(length(idx.drm), length(idx.prm))

  # make the empty list and data frame to contain icc and tcc
  prob.cats <- vector("list", n.item)
  icc_df <- array(NA, c(n.theta, n.item))
  names(prob.cats) <- x$id
  colnames(icc_df) <- x$id

  # For DRM items
  if (!is.null(idx.drm)) {
    # compute the probabilities of correct answers
    p.mat <- drm(
      theta = theta, a = elm_item$pars[idx.drm, 1],
      b = elm_item$pars[idx.drm, 2], g = elm_item$pars[idx.drm, 3],
      D = D
    )
    p.vec <- c(p.mat)
    q.vec <- 1 - p.vec

    # split the probabilities into each item group
    prob.drm <-
      split.data.frame(
        cbind(resp.0 = q.vec, resp.1 = p.vec),
        rep(idx.drm, each = n.theta)
      )

    # fill the empty list
    prob.cats[idx.drm] <- prob.drm

    # insert icc
    icc_df[, idx.drm] <- p.mat
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])

    # check the maximum score category
    max.cats <- max(elm_item$cats) - 1
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]
      cat.tmp <- elm_item$cats[lg.prm]

      # reorder the index of the PRMs
      idx.tmp <- elm_item$item[[mod]]

      # compute the probabilities of endorsing each score category
      P.all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          grad = FALSE, info = FALSE
        )$P
      colnames(P.all) <- paste0("resp.", 0:max.cats)
      prob.cats[idx.tmp] <-
        split.data.frame(
          P.all,
          rep(1:length(a), n.theta)
        ) %>%
        purrr::map2(
          .y = cat.tmp,
          .f = ~ {
            .x[, 1:(.y)]
          }
        )

      # insert icc
      icc_df[, idx.tmp] <-
        matrix(P.all %*% c(0:max.cats), nrow = n.theta, byrow = TRUE)
    }
  }

  # compute tcc
  tcc.vec <- Rfast::rowsums(icc_df)

  # return results
  rst <- list(prob.cats = prob.cats, icc = icc_df, tcc = tcc.vec, theta = theta)
  class(rst) <- "traceline"
  rst
}


# This function computes item characteristic functions
# This function is used in item parameter estimation
#' @importFrom purrr map2
#' @importFrom Rfast rowsums
#' @import dplyr
trace <- function(elm_item, theta, D = 1, tcc = TRUE) {
  # set indices of the DRM and PRM items
  idx.all <- idxfinder(elm_item)
  idx.drm <- idx.all$idx.drm
  idx.prm <- idx.all$idx.prm

  # count the number of thetas
  n.theta <- length(theta)

  # count the total number of items
  n.item <- sum(length(idx.drm), length(idx.prm))

  # make the empty list and data frame to contain icc and tcc
  prob.cats <- vector("list", n.item)
  if (tcc) {
    # icc_df <- matrix(NA, nrow=n.theta, ncol=n.item)
    icc_df <- array(NA, c(n.theta, n.item))
  } else {
    icc_df <- NULL
  }

  # For DRM items
  if (!is.null(idx.drm)) {
    # compute the probabilities of correct answers
    p.mat <- drm(
      theta = theta, a = elm_item$pars[idx.drm, 1],
      b = elm_item$pars[idx.drm, 2], g = elm_item$pars[idx.drm, 3],
      D = D
    )
    p.vec <- c(p.mat)
    q.vec <- 1 - p.vec

    # split the probabilities into each item group
    # prob.drm <-
    #   split.data.frame(cbind(resp.0 = q.vec, resp.1 = p.vec),
    #                    rep(idx.drm, each=n.theta))
    prob.drm <-
      split.data.frame(
        cbind(q.vec, p.vec),
        rep(idx.drm, each = n.theta)
      )

    # fill the empty list
    prob.cats[idx.drm] <- prob.drm

    # if tcc = TRUE
    if (tcc) {
      icc_df[, idx.drm] <- p.mat
    }
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])

    # check the maximum score category
    max.cats <- max(elm_item$cats) - 1
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]
      cat.tmp <- elm_item$cats[lg.prm]

      # reorder the index of the PRMs
      idx.tmp <- elm_item$item[[mod]]

      # compute the probabilities of endorsing each score category
      P.all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          grad = FALSE, info = FALSE
        )$P
      prob.cats[idx.tmp] <-
        split.data.frame(
          P.all,
          rep(1:length(a), n.theta)
        ) %>%
        purrr::map2(
          .y = cat.tmp,
          .f = ~ {
            .x[, 1:(.y)]
          }
        )

      # if tcc = TRUE
      if (tcc) {
        icc_df[, idx.tmp] <-
          matrix(P.all %*% c(0:max.cats), nrow = n.theta, byrow = TRUE)
      }
    }
  }

  # if tcc = TRUE
  if (tcc) {
    tcc.vec <- Rfast::rowsums(icc_df)
  } else {
    tcc.vec <- NULL
  }

  # return results
  rst <- list(prob.cats = prob.cats, icc = icc_df, tcc = tcc.vec, theta = theta)
  rst
}
