#' Residual-Based DIF Detection Framework Using Categorical Residuals (RDIF-CR)
#'
#' This function computes three statistics of the residual-based DIF detection
#' framework using categorical residuals (RDIF-CR)—\eqn{RDIF_{R}-CR},
#' \eqn{RDIF_{S}-CR}, and \eqn{RDIF_{RS}-CR}—for detecting global differential
#' item functioning (DIF), particularly in polytomously scored items. The RDIF-CR
#' framework evaluates DIF by comparing categorical residual vectors, which are
#' calculated as the difference between a one-hot encoded response vector
#' (with 1 for the selected category and 0 for all others) and the IRT
#' model–predicted probability vector across all score categories.
#' This approach enables fine-grained detection of global DIF patterns at the
#' category level.
#'
#' @inheritParams rdif
#' @param score A numeric vector containing examinees' ability estimates (theta
#'   values). If not provided, [irtQ::crdif()] will estimate ability parameters
#'   internally before computing the RDIF statistics. See [irtQ::est_score()]
#'   for more information on scoring methods. Default is `NULL`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for hypothesis testing using the CRDIF statistics. Default is `0.05`.
#' @param purify.by A character string specifying which RDIF statistic is used
#'   to perform the purification. Available options are `"crdifrs"` for
#'   \eqn{RDIF_{RS}-CR}, `"crdifr"` for \eqn{RDIF_{R}-CR}, and
#'   `"crdifs"` for \eqn{RDIF_{S}-CR}.
#' @param min.resp A positive integer specifying the minimum number of valid
#'   item responses required from an examinee in order to compute an ability
#'   estimate. Default is `NULL`.
#'
#' @details
#' According to Penfield (2010), differential item functioning (DIF) in
#' polytomously scored items can be conceptualized in two forms: global DIF and
#' net DIF. Global DIF refers to differences between groups in the conditional
#' probabilities of responding in specific score categories, thus offering a
#' fine-grained view of DIF at the category level. In contrast, net DIF
#' summarizes these differences into a single value representing the overall
#' impact of DIF on the item’s expected score.
#'
#' The RDIF framework using categorical residuals (RDIF-CR), implemented in
#' [irtQ::crdif()], extends the original residual-based DIF framework proposed
#' by Lim et al. (2022) to detect global DIF in polytomous items. This framework
#' includes three statistics: \eqn{RDIF_{R}-CR}, \eqn{RDIF_{S}-CR}, and
#' \eqn{RDIF_{RS}-CR}, each designed to capture different aspects of group-level
#' differences in categorical response patterns.
#'
#' To illustrate how the RDIF-CR framework operates, consider an item with five
#' ordered score categories (\eqn{k \in \{0,1,2,3,4\}}). Suppose an examinee with
#' latent ability \eqn{\theta} responds with category 2. The one-hot encoded
#' response vector for this response is \eqn{(0,0,1,0,0)^T}. Assume that the IRT
#' model estimates the examinee’s expected score as 2.5 and predicts the
#' category probabilities as \eqn{(0.1, 0.2, 0.4, 0.25, 0.05)^T}.
#' In the RDIF-CR framework, the categorical residual vector is calculated by
#' subtracting the predicted probability vector from the one-hot response vector,
#' resulting in \eqn{(-0.1, -0.2, 0.6, -0.25, -0.05)^T}.
#'
#' In contrast to the RDIF-CR framework, net DIF is assessed using a
#' unidimensional item score residual. In this example, the residual would be
#' \eqn{2 - 2.5 = -0.5}. For detecting net DIF, the [irtQ::rdif()] function
#' should be used instead.
#'
#' Note that for dichotomous items, [irtQ::crdif()] and [irtQ::rdif()] yield
#' identical results. This is because the categorical probability vector for a
#' binary item reduces to a scalar difference, making the global and net DIF
#' evaluations mathematically equivalent.
#'
#' @return This function returns a list containing four main components:
#'
#' \item{no_purify}{A list of sub-objects containing the results of DIF analysis
#' without applying a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{dif_stat}{A data frame summarizing the RDIF-CR analysis results for all
#'     items. Columns include item ID, \eqn{RDIF_{R}-CR}, degrees of freedom,
#'     \eqn{RDIF_{S}-CR}, degrees of freedom, \eqn{RDIF_{RS}-CR}, degrees of freedom,
#'     associated p-values, and sample sizes for the reference and focal groups.}
#'     \item{moments}{A list containing the first and second moments (means and
#'     covariance matrices) of the RDIF-CR statistics. The elements include:
#'       \code{mu.crdifr}, \code{mu.crdifs}, \code{mu.crdifrs} (means), and
#'       \code{cov.crdifr}, \code{cov.crdifs}, \code{cov.crdifrs} (covariances),
#'       each indexed by item ID.}
#'     \item{dif_item}{A list of three numeric vectors identifying items flagged
#'     as DIF based on each statistic: \code{crdifr}, \code{crdifs}, and \code{crdifrs}.}
#'     \item{score}{A numeric vector of ability estimates used to compute the RDIF-CR
#'     statistics. These may be user-supplied or internally estimated.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether a purification procedure
#' was applied.}
#'
#' \item{with_purify}{A list of sub-objects containing the results of DIF analysis
#' after applying the purification procedure. The sub-objects include:
#'   \describe{
#'     \item{purify.by}{A character string indicating the RDIF-CR statistic used for
#'     purification. Possible values are "crdifr", "crdifs", or "crdifrs".}
#'     \item{dif_stat}{A data frame summarizing the final RDIF-CR statistics after
#'     purification. Same structure as in \code{no_purify}, with an additional
#'     column indicating the iteration in which the result was obtained.}
#'     \item{moments}{A list of moments (means and covariance matrices) of the
#'     RDIF-CR statistics for all items, updated based on the final iteration.}
#'     \item{dif_item}{A list of three numeric vectors identifying items flagged
#'     as DIF at any iteration, by each statistic.}
#'     \item{n.iter}{An integer indicating the number of iterations performed during
#'     the purification procedure.}
#'     \item{score}{A numeric vector of updated ability estimates used in the final
#'     iteration.}
#'     \item{complete}{A logical value indicating whether the purification process
#'     converged. If \code{FALSE}, the maximum number of iterations was reached
#'     before convergence.}
#'   }
#' }
#'
#' \item{alpha}{A numeric value indicating the significance level (\eqn{\alpha})
#' used for hypothesis testing with RDIF-CR statistics.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::rdif()], [irtQ::est_irt()], [irtQ::est_item()],
#'  [irtQ::est_score()]
#'
#' @references
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{doi:10.1111/jedm.12313}.
#'
#'   Penfield, R. D. (2010). Distinguishing between net and global DIF in
#'   polytomous items. *Journal of Educational Measurement, 47*(2), 129–149.
#'
#' @examples
#' \donttest{
#'
#' ############################################################################
#' # This example demonstrates how to detect global DIF in polytomous items
#' # using the RDIF-CR framework implemented in `irtQ::crdif()`.
#' # Simulated response data are generated from 5 GRM items with 4 score
#' # categories. DIF is introduced in the 1st and 5th items.
#' ############################################################################
#'
#' ###############################################
#' # (1) Simulate response data with DIF
#' ###############################################
#'
#' set.seed(1)
#'
#' # Generate ability parameters for 1000 examinees in each group
#' # Reference and focal groups follow N(0, 1.5^2)
#' theta_ref <- rnorm(1000, 0, 1.5)
#' theta_foc <- rnorm(1000, 0, 1.5)
#'
#' # Combine abilities from both groups
#' theta_all <- c(theta_ref, theta_foc)
#'
#' # Define item parameters using `irtQ::shape_df()`
#' # Items 1 and 5 are intentionally modified to exhibit DIF
#' par_ref <- irtQ::shape_df(
#'   par.prm = list(
#'     a = c(1, 1, 1, 2, 2),
#'     d = list(c(-2, 0, 1),
#'              c(-2, 0, 2),
#'              c(-2, 0, 1),
#'              c(-1, 0, 2),
#'              c(-2, 0, 0.5))
#'   ),
#'   cats = 4, model = "GRM"
#' )
#'
#' par_foc <- irtQ::shape_df(
#'   par.prm = list(
#'     a = c(2, 1, 1, 2, 0.5),
#'     d = list(c(-0.5, 0, 0.5),
#'              c(-2, 0, 2),
#'              c(-2, 0, 1),
#'              c(-1, 0, 2),
#'              c(-1.5, -1, 0))
#'   ),
#'   cats = 4, model = "GRM"
#' )
#'
#' # Generate response data
#' resp_ref <- irtQ::simdat(x = par_ref, theta = theta_ref, D = 1)
#' resp_foc <- irtQ::simdat(x = par_foc, theta = theta_foc, D = 1)
#'
#' # Combine response data across groups
#' data <- rbind(resp_ref, resp_foc)
#'
#' ###############################################
#' # (2) Estimate item and ability parameters
#' ###############################################
#'
#' # Estimate GRM item parameters using `irtQ::est_irt()`
#' fit_mod <- irtQ::est_irt(data = data, D = 1, model = "GRM", cats = 4)
#'
#' # Extract estimated item parameters
#' x <- fit_mod$par.est
#'
#' # Estimate ability scores using ML method
#' score <- est_score(x = x, data = data, method = "ML")$est.theta
#'
#' ###############################################
#' # (3) Perform RDIF-CR DIF analysis
#' ###############################################
#'
#' # Define group membership: 1 = focal group
#' group <- c(rep(0, 1000), rep(1, 1000))
#'
#' # (a) DIF detection without purification
#' dif_nopuri <- crdif(
#'   x = x, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05
#' )
#' print(dif_nopuri)
#'
#' # (b) DIF detection with purification using RDIF_{R}-CR
#' dif_puri_1 <- crdif(
#'   x = x, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "crdifr"
#' )
#' print(dif_puri_1)
#'
#' # (c) DIF detection with purification using RDIF_{S}-CR
#' dif_puri_2 <- crdif(
#'   x = x, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "crdifs"
#' )
#' print(dif_puri_2)
#'
#' # (d) DIF detection with purification using RDIF_{RS}-CR
#' dif_puri_3 <- crdif(
#'   x = x, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "crdifrs"
#' )
#' print(dif_puri_3)
#'
#' }
#' @export
crdif <- function(x, ...) UseMethod("crdif")

#' @describeIn crdif Default method for computing the three RDIF-CR statistics using
#' a data frame `x` that contains item metadata
#'
#' @export
crdif.default <- function(x,
                          data,
                          score = NULL,
                          group,
                          focal.name,
                          item.skip = NULL,
                          D = 1,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("crdifrs", "crdifr", "crdifs"),
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

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
  } else {
    # if min.resp is not NULL, find the examinees who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- Rfast::rowsums(!is.na(data))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <-
    crdif_one(x = x, data = data, score = score, group = group,
              focal.name = focal.name, item.skip = item.skip, D = D,
              alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- dif_rst$moments
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        crdifr = NA, df.crdifr = NA,
        crdifs = NA, df.crdifs = NA,
        crdifrs = NA, df.crdifrs = NA,
        p.crdifr = NA, p.crdifs = NA, p.crdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE)
    mmt_list <-
      purrr::map(
        .x = dif_rst$moments,
        .f = ~ {
          purrr::map(.x = .x, .f = function(k) {
            NA
          })
        }
      )

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_list_tmp <- no_purify$moments
    dif_pval_tmp <- dif_rst$p_val

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # start the iteration if any item is detected as an DIF item
    if (!is.null(dif_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # in case when at least one DIF item is detected from the no purification DIF analysis
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

        # a flagged item which has the smallest significant p-value
        flag_min <-
          switch(purify.by,
                 crdifr = which.min(dif_pval_tmp$crdifr),
                 crdifs = which.min(dif_pval_tmp$crdifs),
                 crdifrs = which.min(dif_pval_tmp$crdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_min]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:13] <- dif_stat_tmp[flag_min, ]
        dif_stat[del_item, 14] <- i - 1
        mmt_list$mu.crdifr[[del_item]] <- mmt_list_tmp$mu.crdifr[[flag_min]]
        mmt_list$mu.crdifs[[del_item]] <- mmt_list_tmp$mu.crdifs[[flag_min]]
        mmt_list$mu.crdifrs[[del_item]] <- mmt_list_tmp$mu.crdifrs[[flag_min]]
        mmt_list$cov.crdifr[[del_item]] <- mmt_list_tmp$cov.crdifr[[flag_min]]
        mmt_list$cov.crdifs[[del_item]] <- mmt_list_tmp$cov.crdifs[[flag_min]]
        mmt_list$cov.crdifrs[[del_item]] <- mmt_list_tmp$cov.crdifrs[[flag_min]]

        # refine the leftover items
        item_num <- item_num[-flag_min]

        # remove the detected DIF item data which has the smallest p-value from the item metadata
        x_puri <- x_puri[-flag_min, ]

        # remove the detected DIF item data which has the smallest p-value from the response data
        data_puri <- data_puri[, -flag_min]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range,
          norm.prior = norm.prior, nquad = nquad, weights = weights,
          ncore = ncore, ...)$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <-
          crdif_one(
            x = x_puri, data = data_puri, score = score_puri, group = group,
            focal.name = focal.name, item.skip = item.skip.puri, D = D,
            alpha = alpha)

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_list_tmp <- dif_rst_tmp$moments
        dif_pval_tmp <- dif_rst_tmp$p_val

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:13] <- dif_stat_tmp
          dif_stat[item_num, 14] <- i
          mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
          mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
          mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
          mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
          mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
          mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs

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
      if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:13] <- dif_stat_tmp
        dif_stat[item_num, 14] <- i
        mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
        mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
        mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
        mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
        mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
        mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # record the final DIF detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- dif_stat
      with_purify$moments <- mmt_list
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- no_purify$moments
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "crdif"
  rst$call <- cl
  rst
}

#' @describeIn crdif An object created by the function [irtQ::est_irt()].
#'
#' @export
#'
crdif.est_irt <- function(x,
                          score = NULL,
                          group,
                          focal.name,
                          item.skip = NULL,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("crdifrs", "crdifr", "crdifs"),
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

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
  } else {
    # if min.resp is not NULL, find the examinees who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- Rfast::rowsums(!is.na(data))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <-
    crdif_one(x = x, data = data, score = score, group = group,
              focal.name = focal.name, item.skip = item.skip, D = D,
              alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- dif_rst$moments
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        crdifr = NA, df.crdifr = NA,
        crdifs = NA, df.crdifs = NA,
        crdifrs = NA, df.crdifrs = NA,
        p.crdifr = NA, p.crdifs = NA, p.crdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE)
    mmt_list <-
      purrr::map(
        .x = dif_rst$moments,
        .f = ~ {
          purrr::map(.x = .x, .f = function(k) {
            NA
          })
        }
      )

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_list_tmp <- no_purify$moments
    dif_pval_tmp <- dif_rst$p_val

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # start the iteration if any item is detected as an DIF item
    if (!is.null(dif_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # in case when at least one DIF item is detected from the no purification DIF analysis
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

        # a flagged item which has the smallest significant p-value
        flag_min <-
          switch(purify.by,
                 crdifr = which.min(dif_pval_tmp$crdifr),
                 crdifs = which.min(dif_pval_tmp$crdifs),
                 crdifrs = which.min(dif_pval_tmp$crdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_min]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:13] <- dif_stat_tmp[flag_min, ]
        dif_stat[del_item, 14] <- i - 1
        mmt_list$mu.crdifr[[del_item]] <- mmt_list_tmp$mu.crdifr[[flag_min]]
        mmt_list$mu.crdifs[[del_item]] <- mmt_list_tmp$mu.crdifs[[flag_min]]
        mmt_list$mu.crdifrs[[del_item]] <- mmt_list_tmp$mu.crdifrs[[flag_min]]
        mmt_list$cov.crdifr[[del_item]] <- mmt_list_tmp$cov.crdifr[[flag_min]]
        mmt_list$cov.crdifs[[del_item]] <- mmt_list_tmp$cov.crdifs[[flag_min]]
        mmt_list$cov.crdifrs[[del_item]] <- mmt_list_tmp$cov.crdifrs[[flag_min]]

        # refine the leftover items
        item_num <- item_num[-flag_min]

        # remove the detected DIF item data which has the smallest p-value from the item metadata
        x_puri <- x_puri[-flag_min, ]

        # remove the detected DIF item data which has the smallest p-value from the response data
        data_puri <- data_puri[, -flag_min]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range,
          norm.prior = norm.prior, nquad = nquad, weights = weights,
          ncore = ncore, ...)$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <-
          crdif_one(
            x = x_puri, data = data_puri, score = score_puri, group = group,
            focal.name = focal.name, item.skip = item.skip.puri, D = D,
            alpha = alpha)

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_list_tmp <- dif_rst_tmp$moments
        dif_pval_tmp <- dif_rst_tmp$p_val

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:13] <- dif_stat_tmp
          dif_stat[item_num, 14] <- i
          mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
          mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
          mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
          mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
          mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
          mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs

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
      if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:13] <- dif_stat_tmp
        dif_stat[item_num, 14] <- i
        mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
        mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
        mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
        mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
        mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
        mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # record the final DIF detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- dif_stat
      with_purify$moments <- mmt_list
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- no_purify$moments
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "crdif"
  rst$call <- cl
  rst
}

#' @describeIn crdif An object created by the function [irtQ::est_item()].
#'
#' @export
#'
crdif.est_item <- function(x,
                           group,
                           focal.name,
                           item.skip = NULL,
                           alpha = 0.05,
                           missing = NA,
                           purify = FALSE,
                           purify.by = c("crdifrs", "crdifr", "crdifs"),
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

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
  } else {
    # if min.resp is not NULL, find the examinees who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- Rfast::rowsums(!is.na(data))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <-
    crdif_one(x = x, data = data, score = score, group = group,
              focal.name = focal.name, item.skip = item.skip, D = D,
              alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- dif_rst$moments
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        crdifr = NA, df.crdifr = NA,
        crdifs = NA, df.crdifs = NA,
        crdifrs = NA, df.crdifrs = NA,
        p.crdifr = NA, p.crdifs = NA, p.crdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE)
    mmt_list <-
      purrr::map(
        .x = dif_rst$moments,
        .f = ~ {
          purrr::map(.x = .x, .f = function(k) {
            NA
          })
        }
      )

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_list_tmp <- no_purify$moments
    dif_pval_tmp <- dif_rst$p_val

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # start the iteration if any item is detected as an DIF item
    if (!is.null(dif_item_tmp)) {
      # record unique item numbers
      item_num <- 1:nrow(x)

      # in case when at least one DIF item is detected from the no purification DIF analysis
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

        # a flagged item which has the smallest significant p-value
        flag_min <-
          switch(purify.by,
                 crdifr = which.min(dif_pval_tmp$crdifr),
                 crdifs = which.min(dif_pval_tmp$crdifs),
                 crdifrs = which.min(dif_pval_tmp$crdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_min]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:13] <- dif_stat_tmp[flag_min, ]
        dif_stat[del_item, 14] <- i - 1
        mmt_list$mu.crdifr[[del_item]] <- mmt_list_tmp$mu.crdifr[[flag_min]]
        mmt_list$mu.crdifs[[del_item]] <- mmt_list_tmp$mu.crdifs[[flag_min]]
        mmt_list$mu.crdifrs[[del_item]] <- mmt_list_tmp$mu.crdifrs[[flag_min]]
        mmt_list$cov.crdifr[[del_item]] <- mmt_list_tmp$cov.crdifr[[flag_min]]
        mmt_list$cov.crdifs[[del_item]] <- mmt_list_tmp$cov.crdifs[[flag_min]]
        mmt_list$cov.crdifrs[[del_item]] <- mmt_list_tmp$cov.crdifrs[[flag_min]]

        # refine the leftover items
        item_num <- item_num[-flag_min]

        # remove the detected DIF item data which has the smallest p-value from the item metadata
        x_puri <- x_puri[-flag_min, ]

        # remove the detected DIF item data which has the smallest p-value from the response data
        data_puri <- data_puri[, -flag_min]

        # update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range,
          norm.prior = norm.prior, nquad = nquad, weights = weights,
          ncore = ncore, ...)$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <-
          crdif_one(
            x = x_puri, data = data_puri, score = score_puri, group = group,
            focal.name = focal.name, item.skip = item.skip.puri, D = D,
            alpha = alpha)

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_list_tmp <- dif_rst_tmp$moments
        dif_pval_tmp <- dif_rst_tmp$p_val

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:13] <- dif_stat_tmp
          dif_stat[item_num, 14] <- i
          mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
          mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
          mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
          mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
          mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
          mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs

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
      if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:13] <- dif_stat_tmp
        dif_stat[item_num, 14] <- i
        mmt_list$mu.crdifr[item_num] <- mmt_list_tmp$mu.crdifr
        mmt_list$mu.crdifs[item_num] <- mmt_list_tmp$mu.crdifs
        mmt_list$mu.crdifrs[item_num] <- mmt_list_tmp$mu.crdifrs
        mmt_list$cov.crdifr[item_num] <- mmt_list_tmp$cov.crdifr
        mmt_list$cov.crdifs[item_num] <- mmt_list_tmp$cov.crdifs
        mmt_list$cov.crdifrs[item_num] <- mmt_list_tmp$cov.crdifrs
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # record the final DIF detection results with the purification procedure
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- dif_stat
      with_purify$moments <- mmt_list
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- no_purify$moments
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "crdif"
  rst$call <- cl
  rst
}

# This function conducts one iteration of DIF analysis using the IRT residual based statistics
crdif_one <- function(x,
                      data,
                      score,
                      group,
                      focal.name,
                      item.skip = NULL,
                      D = 1,
                      alpha = 0.05) {

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # check the unique score categories
  cats <- cats.mod <- elm_item$cats

  # modify the score category for the binary items
  cats.mod[cats.mod == 2] <- 1

  # check the number of items
  nitem <- length(cats)

  ## ---------------------------------
  # compute the CRDIF statistics
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

  # compute the model probability of score categories
  # and keep the probabilities of the correct answers only for the binary items
  prob_ref <-
    trace(elm_item = elm_item, theta = score_ref, D = D, tcc = FALSE)$prob.cats %>%
    purrr::map2(,
                .y = cats.mod,
                .f = ~ {
                  if (.y == 1) {
                    .x[, 2, drop = FALSE]
                  } else {
                    .x
                  }
                }
    )
  prob_foc <-
    trace(elm_item = elm_item, theta = score_foc, D = D, tcc = FALSE)$prob.cats %>%
    purrr::map2(,
                .y = cats.mod,
                .f = ~ {
                  if (.y == 1) {
                    .x[, 2, drop = FALSE]
                  } else {
                    .x
                  }
                }
    )

  # create a list of the factorized response values across all items
  resp_ref <-
    purrr::map2(
      .x = data.frame(resp_ref, stringsAsFactors = FALSE), .y = cats,
      .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
    )
  resp_foc <-
    purrr::map2(
      .x = data.frame(resp_foc, stringsAsFactors = FALSE), .y = cats,
      .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
    )

  # omit the rows including NA values
  prob_ref <-
    purrr::map2(
      .x = prob_ref, .y = resp_ref,
      .f = ~ {
        .x[!is.na(.y), , drop = FALSE]
      }
    )
  prob_foc <-
    purrr::map2(
      .x = prob_foc, .y = resp_foc,
      .f = ~ {
        .x[!is.na(.y), , drop = FALSE]
      }
    )

  # create a contingency table of score categories for each item
  # and then, transform the table to a matrix format
  # also, omit the rows including NAs
  ref_id <- 1:length(score_ref)
  foc_id <- 1:length(score_foc)
  freqtab_ref <-
    purrr::map2(
      .x = resp_ref, .y = cats,
      .f = function(k, y) {
        tmp.tb <-
          stats::xtabs(~ ref_id + k,
                       na.action = stats::na.omit, addNA = FALSE
          ) %>%
          matrix(ncol = y)
        if (y == 2) {
          tmp.tb[, 2, drop = FALSE]
        } else {
          tmp.tb
        }
      }
    )
  freqtab_foc <-
    purrr::map2(
      .x = resp_foc, .y = cats,
      .f = function(k, y) {
        tmp.tb <-
          stats::xtabs(~ foc_id + k,
                       na.action = stats::na.omit, addNA = FALSE
          ) %>%
          matrix(ncol = y)
        if (y == 2) {
          tmp.tb[, 2, drop = FALSE]
        } else {
          tmp.tb
        }
      }
    )

  # compute the raw residuals
  resid_ref <- purrr::map2(.x = freqtab_ref, .y = prob_ref, .f = ~ {
    .x - .y
  })
  resid_foc <- purrr::map2(.x = freqtab_foc, .y = prob_foc, .f = ~ {
    .x - .y
  })

  # compute the residual-based DIF statistic
  rdifr <-
    purrr::map2(
      .x = resid_foc, .y = resid_ref,
      .f = ~ {
        colMeans(.x, na.rm = TRUE) - colMeans(.y, na.rm = TRUE)
      }
    )
  rdifs <-
    purrr::map2(
      .x = resid_foc, .y = resid_ref,
      .f = ~ {
        colMeans(.x^2, na.rm = TRUE) - colMeans(.y^2, na.rm = TRUE)
      }
    )

  # compute the means and variances of the two statistics for the hypothesis testing
  moments <- resid_moments_cr(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats.mod = cats.mod
  )
  mu_rdifr <- moments$mu.rdifr
  mu_rdifs <- moments$mu.rdifs
  cov_rdifr <- moments$cov.rdifr
  cov_rdifs <- moments$cov.rdifs
  cov_rdifrs <- moments$cov.rdifrs

  # compute the chi-square statistics
  chisq_r <- c()
  chisq_s <- c()
  chisq_rs <- c()
  for (i in 1:nitem) {
    if (i %in% all_miss) {
      chisq_r[i] <- NaN
      chisq_s[i] <- NaN
      chisq_rs[i] <- NaN
    } else {
      # covariance matrices for rdifr and rdifs
      cov_r <- cov_rdifr[[i]]
      cov_s <- cov_rdifs[[i]]
      cov_rs <- cov_rdifrs[[i]]

      # create a vector of mean rdifr and mean rdifs
      mu_r <- cbind(mu_rdifr[[i]])
      mu_s <- cbind(mu_rdifs[[i]])
      mu_rs <- cbind(c(mu_rdifr[[i]], mu_rdifs[[i]]))

      # create a vector of rdifr and rdifs
      est_r <- cbind(rdifr[[i]])
      est_s <- cbind(rdifs[[i]])
      est_rs <- cbind(c(rdifr[[i]], rdifs[[i]]))

      # compute the inverse of covariance matrix
      inv_cov_r <- suppressWarnings(tryCatch(
        {
          solve(cov_r, tol = 1e-15)
        },
        error = function(e) {
          NULL
        }
      ))
      inv_cov_s <- suppressWarnings(tryCatch(
        {
          solve(cov_s, tol = 1e-15)
        },
        error = function(e) {
          NULL
        }
      ))
      inv_cov_rs <- suppressWarnings(tryCatch(
        {
          solve(cov_rs, tol = 1e-15)
        },
        error = function(e) {
          NULL
        }
      ))
      if (is.null(inv_cov_r)) {
        inv_cov_r <- suppressWarnings(tryCatch(
          {
            solve(cov_r + 1e-15, tol = 1e-25)
          },
          error = function(e) {
            NULL
          }
        ))
        if (is.null(inv_cov_r)) {
          inv_cov_r <- suppressWarnings(tryCatch(
            {
              solve(cov_r + 1e-10, tol = 1e-25)
            },
            error = function(e) {
              NULL
            }
          ))
        }
      }
      if (is.null(inv_cov_s)) {
        inv_cov_s <- suppressWarnings(tryCatch(
          {
            solve(cov_s + 1e-15, tol = 1e-25)
          },
          error = function(e) {
            NULL
          }
        ))
        if (is.null(inv_cov_s)) {
          inv_cov_s <- suppressWarnings(tryCatch(
            {
              solve(cov_s + 1e-10, tol = 1e-25)
            },
            error = function(e) {
              NULL
            }
          ))
        }
      }
      if (is.null(inv_cov_rs)) {
        inv_cov_rs <- suppressWarnings(tryCatch(
          {
            solve(cov_rs + 1e-15, tol = 1e-25)
          },
          error = function(e) {
            NULL
          }
        ))
        if (is.null(inv_cov_rs)) {
          inv_cov_rs <- suppressWarnings(tryCatch(
            {
              solve(cov_rs + 1e-10, tol = 1e-25)
            },
            error = function(e) {
              NULL
            }
          ))
        }
      }

      # compute the chi-square statistic
      chisq_r[i] <-
        as.numeric(t(est_r - mu_r) %*% inv_cov_r %*% (est_r - mu_r))
      chisq_s[i] <-
        as.numeric(t(est_s - mu_s) %*% inv_cov_s %*% (est_s - mu_s))
      chisq_rs[i] <-
        as.numeric(t(est_rs - mu_rs) %*% inv_cov_rs %*% (est_rs - mu_rs))
    }
  }

  # degree of freedom for chi-square statistics
  df.1 <- cats
  df.1[df.1 == 2] <- 1
  df.2 <- df.1  * 2

  # calculate p-values for all three statistics
  # p_crdifr <- round(stats::pchisq(chisq_r, df=df.1, lower.tail=FALSE), 4)
  # p_crdifs <- round(stats::pchisq(chisq_s, df=df.1, lower.tail=FALSE), 4)
  # p_crdifrs <- round(stats::pchisq(chisq_rs, df=df.2, lower.tail=FALSE), 4)
  p_crdifr <- stats::pchisq(chisq_r, df = df.1, lower.tail = FALSE)
  p_crdifs <- stats::pchisq(chisq_s, df = df.1, lower.tail = FALSE)
  p_crdifrs <- stats::pchisq(chisq_rs, df = df.2, lower.tail = FALSE)
  p_val <- list(crdifr = p_crdifr, crdifs = p_crdifs, crdifrs = p_crdifrs)

  # compute total sample size
  n_total <- n_foc + n_ref

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      crdifr = round(chisq_r, 4),
      df.crdifr = df.1,
      crdifs = round(chisq_s, 4),
      df.crdifs = df.1,
      crdifrs = round(chisq_rs, 4),
      df.crdifrs = df.2,
      p.crdifr = round(p_crdifr, 4),
      p.crdifs = round(p_crdifs, 4),
      p.crdifrs = round(p_crdifrs, 4),
      n.ref = n_ref,
      n.foc = n_foc,
      n.total = n_total,
      stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # data frame of the moments
  names(mu_rdifr) <-
    names(mu_rdifs) <-
    names(cov_rdifr) <-
    names(cov_rdifs) <-
    names(cov_rdifrs) <- x$id
  mmt_list <- list(
    mu.crdifr = mu_rdifr,
    mu.crdifs = mu_rdifs,
    mu.crdifrs = purrr::map2(
      .x = mu_rdifr, .y = mu_rdifs,
      .f = ~ {
        c(.x, .y)
      }),
    cov.crdifr = cov_rdifr,
    cov.crdifs = cov_rdifs,
    cov.crdifrs = cov_rdifrs
  )

  # when there are items that should be skipped for the DIF analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2:10] <- NA
    p_crdifr[item.skip] <- NA
    p_crdifs[item.skip] <- NA
    p_crdifrs[item.skip] <- NA
    p_val <- p_val %>%
      purrr::map(.f = ~{
        .x[item.skip] <- NA
        .x})
    mmt_list$mu.crdifr[item.skip] <- NA
    mmt_list$mu.crdifs[item.skip] <- NA
    mmt_list$mu.crdifrs[item.skip] <- NA
    mmt_list$cov.crdifr[item.skip] <- NA
    mmt_list$cov.crdifs[item.skip] <- NA
    mmt_list$cov.crdifrs[item.skip] <- NA
  }

  # find the flagged items
  dif_item_crdifr <- which(p_crdifr <= alpha)
  dif_item_crdifs <- which(p_crdifs <= alpha)
  dif_item_crdifrs <- which(p_crdifrs <= alpha)
  if (length(dif_item_crdifr) == 0) dif_item_crdifr <- NULL
  if (length(dif_item_crdifs) == 0) dif_item_crdifs <- NULL
  if (length(dif_item_crdifrs) == 0) dif_item_crdifrs <- NULL

  # summarize the results
  rst <- list(
    dif_stat = stat_df,
    dif_item = list(
      crdifr = dif_item_crdifr,
      crdifs = dif_item_crdifs,
      crdifrs = dif_item_crdifrs
    ),
    moments = mmt_list,
    alpha = alpha,
    p_val = p_val
  )

  # return the results
  rst

}

# This function computes the mean and covariance of the RDIF statistics across all category response data
resid_moments_cr <- function(p_ref, p_foc, n_ref, n_foc, resp_ref, resp_foc, cats.mod) {
  # check the number of items
  nitem <- length(cats.mod)

  # check the number of rows and columns in the two-group response data
  nrow_ref <- length(resp_ref[[1]])
  nrow_foc <- length(resp_foc[[1]])

  # compute the mean, variance, and covariance of rdif statistics
  mu_rdifr <-
    mu_rdifs <-
    sigma2_rdifr <-
    sigma2_rdifs <-
    cov_rdifr <-
    cov_rdifs <-
    cov_rs_diag <-
    cov_rs <-
    cov_all <- vector("list", nitem)
  for (i in 1:nitem) {
    # extract the requirements to compute the moments
    p.f <- p_foc[[i]]
    p.r <- p_ref[[i]]
    q.f <- 1 - p.f
    q.r <- 1 - p.r
    p.f2 <- p.f^2
    p.r2 <- p.r^2
    n.f <- n_foc[i]
    n.r <- n_ref[i]
    n.f2 <- n.f^2
    n.r2 <- n.r^2

    # compute the moments of statistics
    # (a) means of the rdifr and rdifs
    mu_rdifr[[i]] <- rep(0, cats.mod[i])
    mu_rdifs[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j], na.rm = TRUE) / n.f) -
            (sum(p.r[, j] * q.r[, j], na.rm = TRUE) / n.r)
        }
      )

    # (b) variances of the rdifr and rdifs
    sigma2_rdifr[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j], na.rm = TRUE) / n.f2) +
            (sum(p.r[, j] * q.r[, j], na.rm = TRUE) / n.r2)
        }
      )
    sigma2_rdifs[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j] * (1 - 2 * p.f[, j])^2, na.rm = TRUE) / n.f2) +
            (sum(p.r[, j] * q.r[, j] * (1 - 2 * p.r[, j])^2, na.rm = TRUE) / n.r2)
        }
      )

    # (c) covariances between the rdifr and rdifs for the same score categories
    cov_rs_diag[[i]] <-
      purrr::map_dbl(
        .x = 1:cats.mod[i],
        .f = function(j) {
          (sum(p.f[, j] * q.f[, j] * (1 - 2 * p.f[, j]), na.rm = TRUE) / n.f2) +
            (sum(p.r[, j] * q.r[, j] * (1 - 2 * p.r[, j]), na.rm = TRUE) / n.r2)
        }
      )

    # (d) covariances of rdifr and rdifs for the different score categories
    if (cats.mod[i] == 1) {
      cov_rdifr[[i]] <- matrix(sigma2_rdifr[[i]])
      cov_rdifs[[i]] <- matrix(sigma2_rdifs[[i]])
      cov_rs[[i]] <- matrix(cov_rs_diag[[i]])
    } else {
      cprod.p.f <- crossprod(p.f)
      cprod.p.r <- crossprod(p.r)
      cprod.p.f2 <- crossprod(p.f2)
      cprod.p.r2 <- crossprod(p.r2)
      cprod.pp2.f <- crossprod(x = p.f, y = p.f2)
      cprod.pp2.r <- crossprod(x = p.r, y = p.r2)
      # note that
      # t(p.f2) %*% p.f == crossprod(x = p.f2, y = p.f)
      # t(p.f) %*% p.f2 == crossprod(x = p.f, y = p.f2)
      # t(cprod.pp2.f) == t(p.f2) %*% p.f

      cov.r.tmp <- -{
        cprod.p.f / n.f2 + cprod.p.r / n.r2
      }
      diag(cov.r.tmp) <- sigma2_rdifr[[i]]
      cov_rdifr[[i]] <- cov.r.tmp

      cov.s.tmp <-
        {
          (-cprod.p.f - 4 * cprod.p.f2 + 2 * (t(cprod.pp2.f) + cprod.pp2.f)) / n.f2
        } + {
          (-cprod.p.r - 4 * cprod.p.r2 + 2 * (t(cprod.pp2.r) + cprod.pp2.r)) / n.r2
        }
      diag(cov.s.tmp) <- sigma2_rdifs[[i]]
      cov_rdifs[[i]] <- cov.s.tmp

      cov.rs.tmp <-
        {
          (-cprod.p.f + 2 * cprod.pp2.f) / n.f2
        } + {
          (-cprod.p.r + 2 * cprod.pp2.r) / n.r2
        }
      diag(cov.rs.tmp) <- cov_rs_diag[[i]]
      cov_rs[[i]] <- cov.rs.tmp
    }

    # (e) covariances of both rdifr and rdifs for the different score categories
    cov_all[[i]] <-
      rbind(
        cbind(cov_rdifr[[i]], cov_rs[[i]]),
        cbind(t(cov_rs[[i]]), cov_rdifs[[i]])
      )
  }

  # return the results
  rst <-
    list(
      mu.rdifr = mu_rdifr, mu.rdifs = mu_rdifs,
      cov.rdifr = cov_rdifr, cov.rdifs = cov_rdifs,
      cov.rdifrs = cov_all
    )
  rst
}
