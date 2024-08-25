#' Compute Item/Test Characteristic Functions
#'
#' @description This function computes the item category probabilities, item characteristic function, and
#' test characteristic function given a set of theta values. The returned object of this function can be used
#' to draw the item or test characteristic curve using the function \code{\link{plot.traceline}}.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object
#' of class \code{\link{est_item}} obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}}
#' obtained from the function \code{\link{est_irt}}. See \code{\link{irtfit}}, \code{\link{info}}, or \code{\link{simdat}}
#' for more details about the item metadata. The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}.
#' @param theta A vector of theta values.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#'          Default is 1.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return This function returns an object of class \code{\link{traceline}}. This object contains a list containing
#' the item category probabilities, item characteristic function, and test characteristic function.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{plot.traceline}}, \code{\link{est_item}}
#'
#' @examples
#' ## example
#' ## using a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # read item parameters and transform them to item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-3, 3, 0.5)
#'
#' # compute the item category probabilities and item/test
#' # characteristic functions given the theta values
#' traceline(x = test_flex, theta, D = 1)
#'
#' @export
traceline <- function(x, ...) UseMethod("traceline")

#' @describeIn traceline Default method to compute the item category probabilities, item characteristic function, and
#' test characteristic function for a data frame \code{x} containing the item metadata.
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



#' @describeIn traceline An object created by the function \code{\link{est_item}}.
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


#' @describeIn traceline An object created by the function \code{\link{est_irt}}.
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
