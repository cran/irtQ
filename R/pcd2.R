#' Pseudo-count D2 method
#'
#' This function calculates the Pseudo-count \eqn{D^{2}} statistic
#' to evaluate item parameter drift, as described by Cappaert et al. (2018) and
#' Stone (2000). The Pseudo-count \eqn{D^{2}} statistic is designed to detect
#' item parameter drift efficiently without requiring item recalibration, making
#' it especially valuable in computerized adaptive testing (CAT) environments.
#' This method compares observed and expected response frequencies across
#' quadrature points, which represent latent ability levels. The expected
#' frequencies are computed using the posterior distribution of each examinee's
#' ability (Stone, 2000), providing a robust and sensitive measure of item
#' parameter drift, ensuring the stability and accuracy of the test over time.
#'
#' @inheritParams rdif
#' @inheritParams est_irt
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function.
#' @param item.skip A numeric vector of item indices to exclude from IPD analysis.
#'  If `NULL`, all items are included. Useful for omitting specific items based on
#'  prior insights.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and their corresponding weights (in the second
#'   column) for the latent variable prior distribution. If not `NULL`, the
#'   scale of the latent ability distribution is fixed to match the scale of the
#'   provided quadrature points and weights. The weights and points can be
#'   conveniently generated using the function [irtQ::gen.weight()].
#'
#'   If `NULL`, a normal prior density is used instead, based on the
#'   information provided in the `Quadrature`, `group.mean`, and `group.var`
#'   arguments. Default is `NULL`.
#' @param group.mean A numeric value specifying the mean of the latent variable
#'   prior distribution when `weights = NULL`. Default is 0. This value is fixed
#'   to resolve the indeterminacy of the item parameter scale during
#'   calibration.
#' @param group.var A positive numeric value specifying the variance of the
#'   latent variable prior distribution when `weights = NULL`. Default is 1.
#'   This value is fixed to resolve the indeterminacy of the item parameter
#'   scale during calibration.
#' @param crit.val A critical value applied in hypothesis testing using
#'   the Pseudo-count \eqn{D^{2}} statistic. Default is `NULL`.
#' @param min.resp A positive integer specifying the minimum required number of
#'   responses for each evaluated item. Defaults to `NULL`.
#'
#' @details
#' The Pseudo-count \eqn{D^{2}} statistic quantifies item parameter drift (IPD) by
#' computing the weighted squared differences between the observed and expected
#' response frequencies for each score category across ability levels. The expected
#' frequencies are determined using the posterior distribution of each examinee's
#' ability (Stone, 2000).
#'
#' The Pseudo-count \eqn{D^{2}} statistic is calculated as:
#' \deqn{
#' Pseudo-count D^{2} = \sum_{k=1}^{Q} \left( \frac{r_{0k} + r_{1k}}{N}\right)
#' \left( \frac{r_{1k}}{r_{0k} + r_{1k}} - E_{1k} \right)^2
#' }
#'
#' where \eqn{r_{0k}} and \eqn{r_{1k}} are the pseudo-counts for the incorrect
#' and correct responses at each ability level \eqn{k}, \eqn{E_{1k}} is the
#' expected proportion of correct responses at each ability level \eqn{k},
#' calculated using item parameters from the item bank, and \eqn{N} is the total
#' count of examinees who received each item.
#'
#' \strong{Critical Value (`crit.val`)}:
#' The `crit.val` argument specifies the threshold used to flag an item as
#' exhibiting potential parameter drift. If an item's Pseudo-count \eqn{D^{2}}
#' value exceeds this threshold, it is identified as a drifted item. If
#' `crit.val = NULL`, the function reports the raw statistic without flagging.
#'
#' \strong{Minimum Response Count (`min.resp`)}:
#' The `min.resp` argument sets a lower bound on the number of responses required
#' for an item to be included in the analysis. Items with fewer responses than
#' `min.resp` are automatically excluded by replacing all their responses
#' with `NA`. This avoids unreliable estimates based on small sample sizes.
#'
#' \strong{Purification Procedure}:
#' Although Cappaert et al. (2018) did not incorporate purification into their method,
#' [irtQ::pcd2()] implements an optional iterative purification process similar
#' to Lim et al. (2022). When `purify = TRUE` and a `crit.val` is provided:
#' \itemize{
#'   \item The procedure begins by identifying items flagged for drift using the initial
#'   Pseudo-count \eqn{D^{2}} statistics.
#'   \item In each subsequent iteration, the item with the highest flagged Pseudo-count
#'   \eqn{D^{2}} value is removed from the item set, and the statistics are recalculated
#'   using only the remaining items.
#'   \item The process continues until no additional items are flagged or the
#'   number of iterations reaches `max.iter`.
#'   \item All flagged items and statistics are saved, and convergence status is reported.
#' }
#' This process ensures that drift detection is not distorted by already-flagged items,
#' improving the robustness of the results.
#'
#' @return This function returns a list containing four main components:
#'
#' \item{no_purify}{A list containing the results of Pseudo-count \eqn{D^{2}} analysis
#' without applying the purification procedure. It includes:
#'   \describe{
#'     \item{ipd_stat}{A data frame summarizing the Pseudo-count \eqn{D^{2}} statistics
#'     for all items. The columns include:
#'     `id` (item ID),
#'     `pcd2` (the computed \eqn{D^{2}} value), and
#'     `N` (the number of valid examinee responses per item).}
#'     \item{ipd_item}{A numeric vector of item indices that were flagged as exhibiting
#'     item parameter drift (IPD), based on the specified critical value `crit.val`.
#'     If no items are flagged or `crit.val = NULL`, this is `NULL`.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether the iterative purification
#' procedure was applied (`TRUE`) or not (`FALSE`).}
#'
#' \item{with_purify}{A list containing the results of Pseudo-count \eqn{D^{2}} analysis
#' after applying the purification procedure. This list is populated only when both
#' `purify = TRUE` and `crit.val` is not `NULL`. It includes:
#'   \describe{
#'     \item{ipd_stat}{A data frame reporting the final Pseudo-count \eqn{D^{2}} statistics
#'     after purification. Columns include:
#'     `id` (item ID),
#'     `pcd2` (the computed \eqn{D^{2}} value),
#'     `N` (the number of valid responses), and
#'     `n.iter` (the iteration number in which each item was evaluated).}
#'     \item{ipd_item}{A numeric vector of item indices flagged as IPD items during
#'     purification. Items are ordered by the iteration in which they were flagged.}
#'     \item{n.iter}{An integer indicating the number of purification iterations
#'     completed.}
#'     \item{complete}{A logical value indicating whether the purification procedure
#'     converged before reaching the maximum number of iterations (`max.iter`).
#'     If `FALSE`, the iteration limit was reached before convergence.}
#'   }
#' }
#'
#' \item{crit.val}{A numeric value indicating the critical threshold used to flag
#' items for parameter drift. If not specified by the user, this will be `NULL`.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Cappaert, K. J., Wen, Y., & Chang, Y. F. (2018). Evaluating CAT-adjusted
#' approaches for suspected item parameter drift detection. *Measurement:
#' Interdisciplinary Research and Perspectives, 16*(4), 226-238.
#'
#' Stone, C. A. (2000). Monte Carlo based null distribution for an alternative
#' goodness-of-fit test statistic in IRT models. *Journal of educational
#' measurement, 37*(1), 58-75.
#'
#' @examples
#' ## Example 1: No critical value specified
#' ## Compute the Pseudo-count D² statistics for dichotomous items
#' ## Import the "-prm.txt" output file generated by flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Extract metadata for the first 30 3PLM items
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[1:30, 1:6]
#'
#' # Generate abilities for 500 examinees from N(0, 1)
#' set.seed(25)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # Simulate response data using the item metadata and ability values
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' # Compute the Pseudo-count D² statistics (no purification applied)
#' ps_d2 <- pcd2(x = x, data = data)
#' print(ps_d2)
#'
#' ## Example 2: Applying a critical value with purification
#' # Compute the Pseudo-count D² statistics with purification enabled
#' ps_d2_puri <- pcd2(x = x, data = data, crit.val = 0.002, purify = TRUE)
#' print(ps_d2_puri)
#'
#' @import dplyr
#' @export
pcd2 <- function(x,
                 data,
                 D = 1,
                 item.skip = NULL,
                 missing = NA,
                 Quadrature = c(49, 6.0),
                 weights = NULL,
                 group.mean = 0.0,
                 group.var = 1.0,
                 crit.val = NULL,
                 min.resp = NULL,
                 purify = FALSE,
                 max.iter = 10,
                 verbose = TRUE
) {

  # match.call
  cl <- match.call()

  # Transform a data set to matrix
  data <- data.matrix(data)

  # Confirm and correct all item metadata information
  x <- confirm_df(x)

  # Stop when the model includes any polytomous model
  if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

  # Recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # count the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    memo1 <- paste0(paste0("Item(s) ", loc_allmiss, collapse = ", "), " has(have) no item response data. \n")
    warning(memo1, call. = FALSE)
  }

  # If min.resp is not NULL, find the items which have the number of responses
  # less than the specified value (e.g., 5), and then, replace all responses of
  # those items with NA
  if (!is.null(min.resp)) {
    loc_less <- which(n.resp < min.resp & n.resp > 0)
    data[, loc_less] <- NA
    memo2 <- paste0(paste0("Item(s) ", loc_less, collapse = ", "), " is(are) not be analyzed ",
                    "because it(they) did not meet the minimum response count criterion. \n")
    warning(memo2, call. = FALSE)
  }

  # (1) when no purification is implemented
  # Do only one iteration of IPD analysis
  ipd_rst <-
    pcd2_one(x = x, data =data, D = D, item.skip = item.skip, missing = missing,
             Quadrature = Quadrature, weights = weights, group.mean = group.mean,
             group.var = group.var, crit.val = crit.val)

  # Create two empty lists to contain the results
  no_purify <- list(ipd_stat = NULL, ipd_item = NULL)
  with_purify <-
    list(
      ipd_stat = NULL, ipd_item = NULL, n.iter = NULL, complete = NULL
    )

  # Record the first IPD detection results into the no purification list
  no_purify$ipd_stat <- ipd_rst$ipd_stat
  if(!is.null(ipd_rst$ipd_item)) {no_purify$ipd_item <- ipd_rst$ipd_item}

  # (2) When purification is implemented
  if (!is.null(crit.val) && purify) {

    # Create an empty vector and empty data frames
    # to contain the detected IPD items and statistics
    ipd_item <- NULL
    ipd_stat <-
      data.frame(id = rep(NA_character_, nrow(x)), pcd2 = NA, N = NA,
                 n.iter = NA, stringsAsFactors = FALSE)

    # Extract the first IPD analysis results
    # and check if at least one IPD item is detected
    ipd_item_tmp <- ipd_rst$ipd_item
    ipd_stat_tmp <- ipd_rst$ipd_stat

    # Copy the response data and item meta data
    x_puri <- x
    data_puri <- data

    # Start the iteration if any item is detected as an IPD item
    if (!is.null(ipd_item_tmp)) {

      # Record unique item numbers
      item_num <- 1:nrow(x)

      # In case when at least one IPD item is detected from the no purification IPD analysis:
      # In this case, the maximum number of iteration must be greater than 0.
      # if not, stop and return an error message
      if (max.iter < 1) stop("The maximum iteration (i.e., max.iter) must be greater than 0 when purify = TRUE.", call. = FALSE)

      # Print a message
      if (verbose) {
        cat("Purification started...", "\n")
      }

      for (i in 1:max.iter) {

        # Print a message
        if (verbose) {
          cat("\r", paste0("Iteration: ", i))
        }

        # A flagged item which has the largest significant IPD statistic
        flag_max <- which.max(ipd_stat_tmp$pcd2)

        # Check an item that is deleted
        del_item <- item_num[flag_max]

        # Add the deleted item as the IPD item
        ipd_item <- c(ipd_item, del_item)

        # Add the IPD statistics for the detected IPD item
        ipd_stat[del_item, 1:3] <- ipd_stat_tmp[flag_max, ]
        ipd_stat[del_item, 4] <- i - 1

        # Refine the leftover items
        item_num <- item_num[-flag_max]

        # Remove the detected IPD item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # Remove the detected IPD item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # Update the locations of the items that should be skipped in the purified data
        if (!is.null(item.skip)) {
          item.skip.puri <- c(1:length(item_num))[item_num %in% item.skip]
        } else {
          item.skip.puri <- NULL
        }

        # Conduct IPD analysis using the purified data
        ipd_rst_tmp <-
          pcd2_one(x = x_puri, data = data_puri, D = D, item.skip = item.skip.puri,
                   missing = missing, Quadrature = Quadrature, weights = weights,
                   group.mean = group.mean, group.var = group.var,
                   crit.val = crit.val)

        # Extract the IPD analysis results
        # and check if at least one IPD item is detected
        ipd_item_tmp <- ipd_rst_tmp$ipd_item
        ipd_stat_tmp <- ipd_rst_tmp$ipd_stat

        # Check if a further IPD item is flagged
        if (is.null(ipd_item_tmp)) {

          # Add no additional IPD item
          ipd_item <- ipd_item

          # Add the IPD statistics for rest of items
          ipd_stat[item_num, 1:3] <- ipd_stat_tmp
          ipd_stat[item_num, 4] <- i

          break
        }
      }

      # Print a message
      if (verbose) {
        cat("", "\n")
      }

      # Record the actual number of iteration
      n_iter <- i

      # If the iteration reached out the maximum number of iteration but the purification is incomplete,
      # then, return a warning message
      if(max.iter == n_iter & !is.null(ipd_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged IPD item at the last iteration
        ipd_item <- c(ipd_item, item_num[ipd_item_tmp])

        # add the IPD statistics for rest of items
        ipd_stat[item_num, 1:3] <- ipd_stat_tmp
        ipd_stat[item_num, 4] <- i
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # Record the final IPD detection results with the purification procedure
      with_purify$ipd_stat <- ipd_stat
      with_purify$ipd_item <- sort(ipd_item)
      with_purify$n.iter <- n_iter
      with_purify$complete <- complete

    } else {

      # In case when no IPD item is detected from the first IPD analysis results
      with_purify$ipd_stat <- cbind(no_purify$ipd_stat, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE

    }
  }

  # Summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, crit.val = crit.val)

  # return the IPD detection results
  # class(rst) <- "pcd2"
  rst$call <- cl
  rst

}


# Pseudo-count D2 method for a single iteration
pcd2_one <- function(x, data, D = 1, item.skip = NULL,
                     missing = NA, Quadrature = c(49, 6.0), weights = NULL,
                     group.mean = 0.0, group.var = 1.0, crit.val = NULL) {

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # extract item ids for all items
  id <- x$id

  # extract score categories for all items
  cats <- x$cats

  # extract model names for all items
  model <- x$model

  # count the total number of item in the response data set
  nitem <- ncol(data)

  # count the total number of examinees
  nstd <- nrow(data)

  # check whether the data have the same number of items as indicated in the x argument
  if (nrow(x) != nitem) {
    stop("The number of items included in 'x' and 'data' must be the same.", call. = FALSE)
  }

  # count the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # create the initial weights of prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create a vector of quad-points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])

    # create a two column data frame to contain the quad-points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
    n.quad <- length(quadpt)
  } else {
    quadpt <- weights[, 1]
    n.quad <- length(quadpt)
    moments.tmp <- cal_moment(node = quadpt, weight = weights[, 2])
    group.mean <- moments.tmp[1]
    group.var <- moments.tmp[2]
  }

  # factorize the response values
  resp <-
    purrr::map2(
      .x = data.frame(data, stringsAsFactors = FALSE), .y = cats,
      .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
    )

  # create a contingency table of score categories for each item
  # and then, transform the table to a matrix format
  std.id <- 1:nstd
  freq.cat <-
    purrr::map(
      .x = resp,
      .f = function(k) {
        stats::xtabs(~ std.id + k,
                     na.action = stats::na.pass, addNA = FALSE
        ) %>%
          # as.numeric() %>%
          matrix(nrow = length(k))
      }
    )

  # delete 'resp' object
  rm(resp, envir = environment(), inherits = FALSE)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # divide the data set for the mixed-item format
  datlist <- divide_data(data = data, idx.item = idx.item, freq.cat = freq.cat)
  data_drm <- cbind(datlist$data_drm_q, datlist$data_drm_p)
  data_prm <- datlist$data_prm
  data_all <- datlist$data_all

  # delete 'datlist' object
  rm(datlist, envir = environment(), inherits = FALSE)

  # find the columns of the frequency matrix corresponding to all items
  cols.item <- cols4item(nitem, cats, NULL)

  # compute the likelihood matrix
  lmat <-
    likelihood(elm_item,
               idx.drm = idx.drm, idx.prm = idx.prm,
               data_drm = data_drm, data_prm = data_prm, theta = weights[, 1], D = D)$L

  # compute posterior distribution for all examinees
  post_dist <- posterior(likehd = lmat, weights = weights, idx.std = NULL)

  # compute the expected frequency of score categories across all items
  # : this is the conditional expectation of item responses with respect
  #   to posterior likelihood distribution
  freq_exp_all <-
    base::as.matrix(Matrix::crossprod(post_dist, data_all))

  # list including the expected frequency of score categories for each item
  freq_exp_item <-
    purrr::map(.x = cols.item$cols.all,
               .f = ~ {freq_exp_all[, .x]})

  # sum of the expected frequencies across all score categories for each item
  freq_exp_sum <-
    purrr::map(.x = freq_exp_item,
               .f = ~ {Rfast::rowsums(.x)})

  # compute the expected probabilities of endorsing each score categories for each item
  exp_prob_item <-
    trace(elm_item = elm_item, theta = weights[, 1], D = D,
          tcc = FALSE)$prob.cats

  # compute the pseudo-count D2 values across all items
  pc_d2 <-
    purrr::pmap_dbl(.l = list(x = freq_exp_item, y = freq_exp_sum,
                              z = n.resp, h = exp_prob_item),
                    .f = function(x, y, z, h) {
                      sum((y / z) * ((x[, 2] / y) - h[, 2])^2)
                    })

  # create a data frame to contain the computed pcd2 statistics
  stat_df <- data.frame(id = id, pcd2 = pc_d2, N = n.resp)

  # when there are items that should be skipped for the IPD analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2] <- NA
    pc_d2[item.skip] <- NA
  }

  # find the flagged items based on the provided critical value
  if(is.null(crit.val)) {
    ipd_item <- NULL
  } else {
    ipd_item <- as.numeric(which(pc_d2 >= crit.val))
    if(length(ipd_item) == 0) {ipd_item <- NULL}
  }

  # return the results
  rst <- list(ipd_stat = stat_df, ipd_item = ipd_item, crit.val = crit.val)
  rst

}
