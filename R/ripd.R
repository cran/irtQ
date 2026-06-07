#' Residual-based Item Parameter Drift (RIPD) Detection Framework
#'
#' This function computes three RIPD statistics—\eqn{RIPD_{R}}, \eqn{RIPD_{S}},
#' and \eqn{RIPD_{RS}}—for each item. \eqn{RIPD_{R}} captures differences in
#' mean raw residuals between groups, which is typically indicative of uniform
#' item parameter drift (IPD). \eqn{RIPD_{S}} captures differences in mean
#' squared residuals between groups, reflecting nonuniform IPD. \eqn{RIPD_{RS}},
#' a combined chi-square-based index, is sensitive to both uniform and
#' nonuniform IPD.
#'
#' @inheritParams rdif
#' @param score A numeric vector containing examinees' ability estimates (theta
#'   values). If not provided, [irtQ::ripd()] will estimate ability parameters
#'   internally before computing the RIPD statistics. See [irtQ::est_score()]
#'   for more information on scoring methods. Default is `NULL`.
#' @param group A numeric or character vector indicating group membership of
#'   examinees. The length of vector should the same with the number of rows in
#'   the response data matrix.
#' @param item.skip A numeric vector of item indices to exclude from IPD
#'   analysis. If NULL, all items are included. Useful for omitting specific
#'   items based on prior insights.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for hypothesis testing using the RIPD statistics. Default is `0.05`.
#' @param missing A value indicating missing values in the response data set.
#'   Default is NA.
#' @param purify.by A character string specifying which RIPD statistic is used
#'   to perform the purification. Available options are "ripdrs" for
#'   \eqn{RIPD_{RS}}, "ripdr" for \eqn{RIPD_{R}}, and "ripds" for
#'   \eqn{RIPD_{S}}.
#' @param min.resp A positive integer specifying the minimum number of valid
#'   item responses required from an examinee in order to compute an ability
#'   estimate. Default is `NULL`.
#'
#' @return This function returns a list containing four main components:
#'
#' \item{no_purify}{A list of sub-objects containing the results of IPD analysis
#' without applying a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{ipd_stat}{A data frame summarizing the RIPD analysis results for all
#'     items. The columns include: item ID, \eqn{RIPD_{R}} statistic, standardized
#'     \eqn{RIPD_{R}}, \eqn{RIPD_{S}} statistic, standardized \eqn{RIPD_{S}},
#'     \eqn{RIPD_{RS}} statistic, p-values for \eqn{RIPD_{R}}, \eqn{RIPD_{S}}, and
#'     \eqn{RIPD_{RS}}, sample sizes for the reference and focal groups, and total
#'     sample size. Note that \eqn{RIPD_{RS}} does not have a standardized value
#'     because it is a \eqn{\chi^{2}}-based statistic.}
#'     \item{moments}{A data frame reporting the first and second moments of the
#'     RIPD statistics. The columns include: item ID, mean and standard deviation
#'     of \eqn{RIPD_{R}}, mean and standard deviation of \eqn{RIPD_{S}}, and the
#'     covariance between \eqn{RIPD_{R}} and \eqn{RIPD_{S}}.}
#'     \item{ipd_item}{A list of three numeric vectors identifying items flagged
#'     as drifting by each RIPD statistic: \eqn{RIPD_{R}}, \eqn{RIPD_{S}}, and
#'     \eqn{RIPD_{RS}}.}
#'     \item{score}{A numeric vector of ability estimates used to compute the RIPD
#'     statistics.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether the purification procedure
#' was applied.}
#'
#' \item{with_purify}{A list of sub-objects containing the results of IPD analysis
#' with a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{purify.by}{A character string indicating the RIPD statistic used for
#'     purification. Possible values are "ripdr", "ripds", and "ripdrs",
#'     corresponding to \eqn{RIPD_{R}}, \eqn{RIPD_{S}}, and \eqn{RIPD_{RS}},
#'     respectively.}
#'     \item{ipd_stat}{A data frame reporting the RIPD analysis results for
#'     all items from the final iteration. Same structure as in \code{no_purify},
#'     with one additional column indicating the iteration number in which each
#'     result was obtained.}
#'     \item{moments}{A data frame reporting the moments of RIPD statistics
#'     from the final iteration. Includes the same columns as in
#'     \code{no_purify}, with an additional column for the iteration number.}
#'     \item{ipd_item}{A list of numeric item indices identified as IPD items
#'     across iterations.}
#'     \item{n.iter}{An integer indicating the total number of iterations
#'     performed during the purification process.}
#'     \item{score}{A numeric vector of purified ability estimates used to
#'     compute the final RIPD statistics.}
#'     \item{complete}{A logical value indicating whether the purification
#'     process converged. If \code{FALSE}, the maximum number of iterations
#'     was reached without convergence.}
#'   }
#' }
#'
#' \item{alpha}{A numeric value indicating the significance level (\eqn{\alpha})
#' used in hypothesis testing for RIPD statistics.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::rdif()], [irtQ::est_irt()], [irtQ::est_item()],
#'   [irtQ::simdat()], [irtQ::shape_df()], [irtQ::est_score()]
#'
#' @references Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based
#'   differential item functioning detection framework in item response theory.
#'   *Journal of Educational Measurement, 59*(1), 80-104.
#'   \doi{doi:10.1111/jedm.12313}.
#'
#'   Lim, H., & H, K. T. (2025, April). *IRT residual-based approach to
#'   detecting item parameter drift in CAT*. Paper presented at the annual
#'   conference of the National Council on Measurement in Education (NCME),
#'   Denver, CO
#'
#'@export
ripd <- function(x, ...) UseMethod("ripd")

#' @describeIn ripd Default method for computing the three RIPD statistics using
#' a data frame `x` that contains item metadata
#'
#'@export
ripd.default <- function(x,
                         data,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         D = 1,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("ripdrs", "ripdr", "ripds"),
                         max.iter = 10,
                         min.resp = NULL,
                         method = "ML",
                         range = c(-5, 5),
                         norm.prior = c(0, 1),
                         nquad = 41,
                         weights = NULL,
                         ncore = 1,
                         verbose = TRUE,
                         ...) {

  # match.call
  cl <- match.call()

  ## ----------------------------------
  ## (1) prepare IPD analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # stop when the model includes any polytomous response data
  if(any(data > 1, na.rm=TRUE)) {
    stop("The current version only supports dichotomous response data.", call.=FALSE)
  }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
  } else {
    # if min.resp is not NULL, find the examinees of the studied group who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- Rfast::rowsums(!is.na(data[group == focal.name, ]))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[group == focal.name, ][loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...)$est.theta
  }

  # a) when no purification is set
  # do only one iteration of IPD analysis
  ipd_rst <- ripd_one(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha
  )

  # create two empty lists to contain the results
  no_purify <- list(ipd_stat = NULL, moments = NULL, ipd_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, ipd_stat = NULL, moments = NULL,
    ipd_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first IPD detection results into the no purification list
  no_purify$ipd_stat <- ipd_rst$ipd_stat
  no_purify$ipd_item <- ipd_rst$ipd_item
  no_purify$moments <- data.frame(
    id = x$id,
    ipd_rst$moments$ripdr[, c(1, 3)],
    ipd_rst$moments$ripds[, c(1, 3)],
    ipd_rst$covariance, stringsAsFactors = FALSE
  )
  names(no_purify$moments) <- c("id", "mu.ripdr", "sigma.ripdr", "mu.ripds",
                                "sigma.ripds", "covariance")
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected IPD items, statistics, and moments
    ipd_item <- NULL
    ipd_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)), ripdr = NA, z.ripdr = NA,
        ripds = NA, z.ripds = NA, ripdrs = NA, p.ripdr = NA, p.ripds = NA, p.ripdrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA, stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.ripdr = NA, sigma.ripdr = NA,
        mu.ripds = NA, sigma.ripds = NA, covariance = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )

    # extract the first IPD analysis results
    # and check if at least one IPD item is detected
    ipd_item_tmp <- ipd_rst$ipd_item[[purify.by]]
    ipd_stat_tmp <- ipd_rst$ipd_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # start the iteration if any item is detected as an IPD item
    if (!is.null(ipd_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # in case when at least one IPD item is detected from the no purification IPD analysis
      # in this case, the maximum number of iteration must be greater than 0.
      # if not, stop and return an error message
      if (max.iter < 1) stop("The maximum iteration (i.e., max.iter) must be greater than 0 when purify = TRUE.", call. = FALSE)

      # print a message
      if (verbose) {
        cat("Purification started...", "\n")
      }

      for (i in 1:max.iter) {
        # print a message
        if (verbose) {
          cat("\r", paste0("Iteration: ", i))
        }

        # a flagged item which has the largest significant IPD statistic
        flag_max <-
          switch(purify.by,
            ripdr = which.max(abs(ipd_stat_tmp$z.ripdr)),
            ripds = which.max(abs(ipd_stat_tmp$z.ripds)),
            ripdrs = which.max(ipd_stat_tmp$ripdrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the IPD item
        ipd_item <- c(ipd_item, del_item)

        # add the IPD statistics and moments for the detected IPD item
        ipd_stat[del_item, 1:12] <- ipd_stat_tmp[flag_max, ]
        ipd_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # remove the detected IPD item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected IPD item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # if min.resp is not NULL, find the examinees of the studied group who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri[group == focal.name, ]))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[group == focal.name, ][loc_less, ] <- NA
        }

        # compute the updated ability estimates for the studied groups
        # after deleting the detected IPD item data
        score_puri <-
          est_score(
            x = x_puri, data = data_puri, D = D, method = method,
            range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
            ncore = ncore, ...)$est.theta

        # replace the scores of the studied group with the purified scores
        score <- score_puri

        # do IPD analysis using the updated ability estimates
        ipd_rst_tmp <- ripd_one(
          x = x_puri, data = data_puri, score = score, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
        )

        # extract the IPD analysis results
        # and check if at least one IPD item is detected
        ipd_item_tmp <- ipd_rst_tmp$ipd_item[[purify.by]]
        ipd_stat_tmp <- ipd_rst_tmp$ipd_stat
        mmt_df_tmp <- data.frame(
          id = ipd_rst_tmp$ipd_stat$id,
          ipd_rst_tmp$moments$ripdr[, c(1, 3)],
          ipd_rst_tmp$moments$ripds[, c(1, 3)],
          ipd_rst_tmp$covariance, stringsAsFactors = FALSE
        )
        names(mmt_df_tmp) <- c("id", "mu.ripdr", "sigma.ripdr", "mu.ripds",
                               "sigma.ripds", "covariance")

        # check if a further IPD item is flagged
        if (is.null(ipd_item_tmp)) {
          # add no additional IPD item
          ipd_item <- ipd_item

          # add the IPD statistics for rest of items
          ipd_stat[item_num, 1:12] <- ipd_stat_tmp
          ipd_stat[item_num, 13] <- i
          mmt_df[item_num, 1:6] <- mmt_df_tmp
          mmt_df[item_num, 7] <- i

          break
        }
      }

      # print a message
      if (verbose) {
        cat("", "\n")
      }

      # record the actual number of iteration
      n_iter <- i

      # if the iteration reached out the maximum number of iteration but the purification is incomplete,
      # then, return a warning message
      if (max.iter == n_iter & !is.null(ipd_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged IPD item at the last iteration
        ipd_item <- c(ipd_item, item_num[ipd_item_tmp])

        # add the IPD statistics for rest of items
        ipd_stat[item_num, 1:12] <- ipd_stat_tmp
        ipd_stat[item_num, 13] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # record the final IPD detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- ipd_stat
      with_purify$moments <- mmt_df
      with_purify$ipd_item <- sort(ipd_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score
      with_purify$complete <- complete
    } else {
      # in case when no IPD item is detected from the first IPD analysis results
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- cbind(no_purify$ipd_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the IPD detection results
  class(rst) <- "ripd"
  rst$call <- cl
  rst
}

#' @describeIn ripd An object created by the function [irtQ::est_irt()].
#'
#'@export
ripd.est_irt <- function(x,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("ripdrs", "ripdr", "ripds"),
                         max.iter = 10,
                         min.resp = NULL,
                         method = "ML",
                         range = c(-5, 5),
                         norm.prior = c(0, 1),
                         nquad = 41,
                         weights = NULL,
                         ncore = 1,
                         verbose = TRUE,
                         ...) {

  # match.call
  cl <- match.call()

  # extract information from an object
  data <- x$data
  D <- x$scale.D
  x <- x$par.est

  ## ----------------------------------
  ## (1) prepare IPD analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # stop when the model includes any polytomous response data
  # if(any(data > 1, na.rm=TRUE)) {
  #   stop("The current version only supports dichotomous response data.", call.=FALSE)
  # }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
  } else {
    # if min.resp is not NULL, find the examinees of the studied group who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- Rfast::rowsums(!is.na(data[group == focal.name, ]))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[group == focal.name, ][loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of IPD analysis
  ipd_rst <- ripd_one(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha
  )

  # create two empty lists to contain the results
  no_purify <- list(ipd_stat = NULL, moments = NULL, ipd_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, ipd_stat = NULL, moments = NULL,
    ipd_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first IPD detection results into the no purification list
  no_purify$ipd_stat <- ipd_rst$ipd_stat
  no_purify$ipd_item <- ipd_rst$ipd_item
  no_purify$moments <- data.frame(
    id = x$id,
    ipd_rst$moments$ripdr[, c(1, 3)],
    ipd_rst$moments$ripds[, c(1, 3)],
    ipd_rst$covariance, stringsAsFactors = FALSE
  )
  names(no_purify$moments) <- c("id", "mu.ripdr", "sigma.ripdr", "mu.ripds", "sigma.ripds", "covariance")
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected IPD items, statistics, and moments
    ipd_item <- NULL
    ipd_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)), ripdr = NA, z.ripdr = NA,
        ripds = NA, z.ripds = NA, ripdrs = NA, p.ripdr = NA, p.ripds = NA, p.ripdrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA, stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.ripdr = NA, sigma.ripdr = NA,
        mu.ripds = NA, sigma.ripds = NA, covariance = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )

    # extract the first IPD analysis results
    # and check if at least one IPD item is detected
    ipd_item_tmp <- ipd_rst$ipd_item[[purify.by]]
    ipd_stat_tmp <- ipd_rst$ipd_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # start the iteration if any item is detected as an IPD item
    if (!is.null(ipd_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # in case when at least one IPD item is detected from the no purification IPD analysis
      # in this case, the maximum number of iteration must be greater than 0.
      # if not, stop and return an error message
      if (max.iter < 1) stop("The maximum iteration (i.e., max.iter) must be greater than 0 when purify = TRUE.", call. = FALSE)

      # print a message
      if (verbose) {
        cat("Purification started...", "\n")
      }

      for (i in 1:max.iter) {
        # print a message
        if (verbose) {
          cat("\r", paste0("Iteration: ", i))
        }

        # a flagged item which has the largest significant IPD statistic
        flag_max <-
          switch(purify.by,
                 ripdr = which.max(abs(ipd_stat_tmp$z.ripdr)),
                 ripds = which.max(abs(ipd_stat_tmp$z.ripds)),
                 ripdrs = which.max(ipd_stat_tmp$ripdrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the IPD item
        ipd_item <- c(ipd_item, del_item)

        # add the IPD statistics and moments for the detected IPD item
        ipd_stat[del_item, 1:12] <- ipd_stat_tmp[flag_max, ]
        ipd_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # remove the detected IPD item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected IPD item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # if min.resp is not NULL, find the examinees of the studied group who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri[group == focal.name, ]))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[group == focal.name, ][loc_less, ] <- NA
        }

        # compute the updated ability estimates for the studied groups
        # after deleting the detected IPD item data
        # score_puri <- est_score(x=x_puri, data=data_puri[group == focal.name, ], D=D, method=method,
        #                         range=range, norm.prior=norm.prior, nquad=nquad, weights=weights,
        #                         ncore=ncore, ...)$est.theta
        score_puri <-
          est_score(
            x = x_puri, data = data_puri, D = D, method = method,
            range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
            ncore = ncore, ...)$est.theta

        # replace the scores of the studied group with the purified scores
        # score[group == focal.name] <- score_puri
        score <- score_puri

        # do IPD analysis using the updated ability estimates
        ipd_rst_tmp <- ripd_one(
          x = x_puri, data = data_puri, score = score, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
        )

        # extract the first IPD analysis results
        # and check if at least one IPD item is detected
        ipd_item_tmp <- ipd_rst_tmp$ipd_item[[purify.by]]
        ipd_stat_tmp <- ipd_rst_tmp$ipd_stat
        mmt_df_tmp <- data.frame(
          id = ipd_rst_tmp$ipd_stat$id,
          ipd_rst_tmp$moments$ripdr[, c(1, 3)],
          ipd_rst_tmp$moments$ripds[, c(1, 3)],
          ipd_rst_tmp$covariance, stringsAsFactors = FALSE
        )
        names(mmt_df_tmp) <- c("id", "mu.ripdr", "sigma.ripdr", "mu.ripds", "sigma.ripds", "covariance")

        # check if a further IPD item is flagged
        if (is.null(ipd_item_tmp)) {
          # add no additional IPD item
          ipd_item <- ipd_item

          # add the IPD statistics for rest of items
          ipd_stat[item_num, 1:12] <- ipd_stat_tmp
          ipd_stat[item_num, 13] <- i
          mmt_df[item_num, 1:6] <- mmt_df_tmp
          mmt_df[item_num, 7] <- i

          break
        }
      }

      # print a message
      if (verbose) {
        cat("", "\n")
      }

      # record the actual number of iteration
      n_iter <- i

      # if the iteration reached out the maximum number of iteration but the purification is incomplete,
      # then, return a warning message
      if (max.iter == n_iter & !is.null(ipd_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged IPD item at the last iteration
        ipd_item <- c(ipd_item, item_num[ipd_item_tmp])

        # add the IPD statistics for rest of items
        ipd_stat[item_num, 1:12] <- ipd_stat_tmp
        ipd_stat[item_num, 13] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # record the final IPD detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- ipd_stat
      with_purify$moments <- mmt_df
      with_purify$ipd_item <- sort(ipd_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score
      with_purify$complete <- complete
    } else {
      # in case when no IPD item is detected from the first IPD analysis results
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- cbind(no_purify$ipd_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the IPD detection results
  class(rst) <- "ripd"
  rst$call <- cl
  rst
}

#' @describeIn ripd An object created by the function [irtQ::est_item()].
#'
#'@export
ripd.est_item <- function(x,
                          group,
                          focal.name,
                          item.skip = NULL,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("ripdrs", "ripdr", "ripds"),
                          max.iter = 10,
                          min.resp = NULL,
                          method = "ML",
                          range = c(-5, 5),
                          norm.prior = c(0, 1),
                          nquad = 41,
                          weights = NULL,
                          ncore = 1,
                          verbose = TRUE,
                          ...) {

  # match.call
  cl <- match.call()

  # extract information from an object
  data <- x$data
  score <- x$score
  D <- x$scale.D
  x <- x$par.est

  ## ----------------------------------
  ## (1) prepare IPD analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # stop when the model includes any polytomous response data
  # if(any(data > 1, na.rm=TRUE)) {
  #   stop("The current version only supports dichotomous response data.", call.=FALSE)
  # }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
  } else {
    # if min.resp is not NULL, find the examinees of the studied group who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- Rfast::rowsums(!is.na(data[group == focal.name, ]))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[group == focal.name, ][loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of IPD analysis
  ipd_rst <- ripd_one(
    x = x, data = data, score = score, group = group, focal.name = focal.name,
    item.skip = item.skip, D = D, alpha = alpha
  )

  # create two empty lists to contain the results
  no_purify <- list(ipd_stat = NULL, moments = NULL, ipd_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, ipd_stat = NULL, moments = NULL,
    ipd_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first IPD detection results into the no purification list
  no_purify$ipd_stat <- ipd_rst$ipd_stat
  no_purify$ipd_item <- ipd_rst$ipd_item
  no_purify$moments <- data.frame(
    id = x$id,
    ipd_rst$moments$ripdr[, c(1, 3)],
    ipd_rst$moments$ripds[, c(1, 3)],
    ipd_rst$covariance, stringsAsFactors = FALSE
  )
  names(no_purify$moments) <- c("id", "mu.ripdr", "sigma.ripdr", "mu.ripds", "sigma.ripds", "covariance")
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected IPD items, statistics, and moments
    ipd_item <- NULL
    ipd_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)), ripdr = NA, z.ripdr = NA,
        ripds = NA, z.ripds = NA, ripdrs = NA, p.ripdr = NA, p.ripds = NA, p.ripdrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA, stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.ripdr = NA, sigma.ripdr = NA,
        mu.ripds = NA, sigma.ripds = NA, covariance = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )

    # extract the first IPD analysis results
    # and check if at least one IPD item is detected
    ipd_item_tmp <- ipd_rst$ipd_item[[purify.by]]
    ipd_stat_tmp <- ipd_rst$ipd_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # start the iteration if any item is detected as an IPD item
    if (!is.null(ipd_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # in case when at least one IPD item is detected from the no purification IPD analysis
      # in this case, the maximum number of iteration must be greater than 0.
      # if not, stop and return an error message
      if (max.iter < 1) stop("The maximum iteration (i.e., max.iter) must be greater than 0 when purify = TRUE.", call. = FALSE)

      # print a message
      if (verbose) {
        cat("Purification started...", "\n")
      }

      for (i in 1:max.iter) {
        # print a message
        if (verbose) {
          cat("\r", paste0("Iteration: ", i))
        }

        # a flagged item which has the largest significant IPD statistic
        flag_max <-
          switch(purify.by,
                 ripdr = which.max(abs(ipd_stat_tmp$z.ripdr)),
                 ripds = which.max(abs(ipd_stat_tmp$z.ripds)),
                 ripdrs = which.max(ipd_stat_tmp$ripdrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the IPD item
        ipd_item <- c(ipd_item, del_item)

        # add the IPD statistics and moments for the detected IPD item
        ipd_stat[del_item, 1:12] <- ipd_stat_tmp[flag_max, ]
        ipd_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # remove the detected IPD item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected IPD item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # if min.resp is not NULL, find the examinees of the studied group who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri[group == focal.name, ]))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[group == focal.name, ][loc_less, ] <- NA
        }

        # compute the updated ability estimates for the studied groups
        # after deleting the detected IPD item data
        # score_puri <- est_score(x=x_puri, data=data_puri[group == focal.name, ], D=D, method=method,
        #                         range=range, norm.prior=norm.prior, nquad=nquad, weights=weights,
        #                         ncore=ncore, ...)$est.theta
        score_puri <-
          est_score(
            x = x_puri, data = data_puri, D = D, method = method,
            range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
            ncore = ncore, ...)$est.theta

        # replace the scores of the studied group with the purified scores
        # score[group == focal.name] <- score_puri
        score <- score_puri

        # do IPD analysis using the updated ability estimates
        ipd_rst_tmp <- ripd_one(
          x = x_puri, data = data_puri, score = score, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
        )

        # extract the first IPD analysis results
        # and check if at least one IPD item is detected
        ipd_item_tmp <- ipd_rst_tmp$ipd_item[[purify.by]]
        ipd_stat_tmp <- ipd_rst_tmp$ipd_stat
        mmt_df_tmp <- data.frame(
          id = ipd_rst_tmp$ipd_stat$id,
          ipd_rst_tmp$moments$ripdr[, c(1, 3)],
          ipd_rst_tmp$moments$ripds[, c(1, 3)],
          ipd_rst_tmp$covariance, stringsAsFactors = FALSE
        )
        names(mmt_df_tmp) <- c("id", "mu.ripdr", "sigma.ripdr", "mu.ripds", "sigma.ripds", "covariance")

        # check if a further IPD item is flagged
        if (is.null(ipd_item_tmp)) {
          # add no additional IPD item
          ipd_item <- ipd_item

          # add the IPD statistics for rest of items
          ipd_stat[item_num, 1:12] <- ipd_stat_tmp
          ipd_stat[item_num, 13] <- i
          mmt_df[item_num, 1:6] <- mmt_df_tmp
          mmt_df[item_num, 7] <- i

          break
        }
      }

      # print a message
      if (verbose) {
        cat("", "\n")
      }

      # record the actual number of iteration
      n_iter <- i

      # if the iteration reached out the maximum number of iteration but the purification is incomplete,
      # then, return a warning message
      if (max.iter == n_iter & !is.null(ipd_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged IPD item at the last iteration
        ipd_item <- c(ipd_item, item_num[ipd_item_tmp])

        # add the IPD statistics for rest of items
        ipd_stat[item_num, 1:12] <- ipd_stat_tmp
        ipd_stat[item_num, 13] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # record the final IPD detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- ipd_stat
      with_purify$moments <- mmt_df
      with_purify$ipd_item <- sort(ipd_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score
      with_purify$complete <- complete
    } else {
      # in case when no IPD item is detected from the first IPD analysis results
      with_purify$purify.by <- purify.by
      with_purify$ipd_stat <- cbind(no_purify$ipd_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the IPD detection results
  class(rst) <- "ripd"
  rst$call <- cl
  rst
}


# This function conducts one iteration of IPD analysis using the IRT residual based statistics
ripd_one <- function(x, data, score, group, focal.name, item.skip = NULL, D = 1, alpha = 0.05) {
  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # check the unique score categories
  cats <- elm_item$cats

  # check the number of items
  nitem <- length(cats)

  ## ---------------------------------
  # compute the two statistics
  ## ---------------------------------
  # find the location of examinees for the reference and the focal groups
  loc_ref <- which(group != focal.name)
  loc_foc <- which(group == focal.name)

  # divide the response data into the two group data
  resp_ref <- data[loc_ref, , drop = FALSE]
  resp_foc <- data[loc_foc, , drop = FALSE]

  # check sample size
  n_ref <- Rfast::colsums(!is.na(resp_ref))
  n_foc <- Rfast::colsums(!is.na(resp_foc))

  # check if an item has all missing data for either of two groups
  all_miss <- sort(unique(c(which(n_ref == 0), which(n_foc == 0))))

  # divide the thetas into the two group data
  score_ref <- score[loc_ref]
  score_foc <- score[loc_foc]

  # compute the model-predicted probability of answering correctly (a.k.a. model-expected item score)
  extscore_ref <- trace(elm_item = elm_item, theta = score_ref, D = D, tcc = TRUE)$icc
  extscore_foc <- trace(elm_item = elm_item, theta = score_foc, D = D, tcc = TRUE)$icc

  # compute the model probability of score categories
  prob_ref <- trace(elm_item = elm_item, theta = score_ref, D = D, tcc = FALSE)$prob.cats
  prob_foc <- trace(elm_item = elm_item, theta = score_foc, D = D, tcc = FALSE)$prob.cats

  # replace NA values into the missing data location
  extscore_ref[is.na(resp_ref)] <- NA
  extscore_foc[is.na(resp_foc)] <- NA

  # compute the raw residuals
  resid_ref <- resp_ref - extscore_ref
  resid_foc <- resp_foc - extscore_foc

  # compute the residual-based IPD statistic
  ripdr <- colMeans(resid_foc, na.rm = TRUE) - colMeans(resid_ref, na.rm = TRUE) # the difference of mean raw residuals (DMRR)
  ripds <- colMeans(resid_foc^2, na.rm = TRUE) - colMeans(resid_ref^2, na.rm = TRUE) # the difference of mean squared residuals (DMSR)

  # compute the means and variances of the two statistics for the hypothesis testing
  moments <- resid_moments(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats = cats
  )
  moments_ripdr <- moments$rdifr
  moments_ripds <- moments$rdifs
  covar <- moments$covariance
  # poolsd <- moments$poolsd

  # compute the chi-square statistics
  chisq <- c()
  for (i in 1:nitem) {
    if (i %in% all_miss) {
      chisq[i] <- NaN
    } else {
      # create a var-covariance matrix between ripdr and ripds
      cov_mat <- array(NA, c(2, 2))

      # replace NAs with the analytically computed covariance
      cov_mat[col(cov_mat) != row(cov_mat)] <- covar[i]

      # replace NAs with the analytically computed variances
      diag(cov_mat) <- c(moments_ripdr[i, 2], moments_ripds[i, 2])

      # create a vector of mean ripdr and mean ripds
      mu_vec <- cbind(moments_ripdr[i, 1], moments_ripds[i, 1])

      # create a vector of ripdr and ripds
      est_mu_vec <- cbind(ripdr[i], ripds[i])

      # compute the chi-square statistic
      inv_cov <- suppressWarnings(tryCatch(
        {
          solve(cov_mat, tol = 1e-200)
        },
        error = function(e) {
          NULL
        }
      ))
      if (is.null(inv_cov)) {
        inv_cov <- suppressWarnings(tryCatch(
          {
            solve(cov_mat + 1e-15, tol = 1e-200)
          },
          error = function(e) {
            NULL
          }
        ))
        if (is.null(inv_cov)) {
          inv_cov <- suppressWarnings(tryCatch(
            {
              solve(cov_mat + 1e-10, tol = 1e-200)
            },
            error = function(e) {
              NULL
            }
          ))
        }
      }
      chisq[i] <- as.numeric((est_mu_vec - mu_vec) %*% inv_cov %*% t(est_mu_vec - mu_vec))
    }
  }

  # standardize the two statistics of ripdr and ripds
  z_stat_ripdr <- (ripdr - moments_ripdr$mu) / moments_ripdr$sigma
  z_stat_ripds <- (ripds - moments_ripds$mu) / moments_ripds$sigma

  # calculate p-values for all three statistics
  p_ripdr <- round(2 * stats::pnorm(q = abs(z_stat_ripdr), mean = 0, sd = 1, lower.tail = FALSE), 4)
  p_ripds <- round(2 * stats::pnorm(q = abs(z_stat_ripds), mean = 0, sd = 1, lower.tail = FALSE), 4)
  p_ripdrs <- round(stats::pchisq(chisq, df = 2, lower.tail = FALSE), 4)

  # compute total sample size
  n_total <- n_foc + n_ref

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      ripdr = round(ripdr, 4), z.ripdr = round(z_stat_ripdr, 4),
      ripds = round(ripds, 4), z.ripds = round(z_stat_ripds, 4),
      ripdrs = round(chisq, 4),
      p.ripdr = p_ripdr, p.ripds = p_ripds, p.ripdrs = p_ripdrs,
      n.ref = n_ref, n.foc = n_foc, n.total = n_total, stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # when there are items that should be skipped for the IPD analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2:9] <- NA
    p_ripdr[item.skip] <- NA
    p_ripds[item.skip] <- NA
    p_ripdrs[item.skip] <- NA
    moments_ripdr[item.skip, ] <- NA
    moments_ripds[item.skip, ] <- NA
    covar[item.skip] <- NA
  }

  # find the flagged items
  ipd_item_ripdr <- as.numeric(which(p_ripdr <= alpha))
  ipd_item_ripds <- as.numeric(which(p_ripds <= alpha))
  ipd_item_ripdrs <- which(p_ripdrs <= alpha)
  if (length(ipd_item_ripdr) == 0) ipd_item_ripdr <- NULL
  if (length(ipd_item_ripds) == 0) ipd_item_ripds <- NULL
  if (length(ipd_item_ripdrs) == 0) ipd_item_ripdrs <- NULL

  # summarize the results
  rst <- list(
    ipd_stat = stat_df,
    ipd_item = list(ripdr = ipd_item_ripdr, ripds = ipd_item_ripds, ripdrs = ipd_item_ripdrs),
    moments = list(ripdr = moments_ripdr, ripds = moments_ripds), covariance = covar, alpha = alpha
  )

  # return the results
  rst
}
