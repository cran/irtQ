#' CATSIB DIF Detection Procedure
#'
#' This function performs DIF analysis on items using the CATSIB procedure
#' (Nandakumar & Roussos, 2004), a modified version of SIBTEST (Shealy & Stout,
#' 1993). The CATSIB procedure is suitable for computerized adaptive testing
#' (CAT) environments. In CATSIB, examinees are matched on IRT-based ability
#' estimates that have been adjusted using a regression correction method
#' (Shealy & Stout, 1993) to reduce statistical bias in the CATSIB statistic
#' caused by impact.
#'
#' @inheritParams rdif
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function.
#' @param score A numeric vector containing examinees' ability estimates (theta
#'   values). If not provided, [irtQ::catsib()] will estimate ability parameters
#'   internally before computing the CATSIB statistics. See [irtQ::est_score()]
#'   for more information on scoring methods. Default is `NULL`.
#' @param se A vector of standard errors corresponding to the ability estimates.
#'   The order of the standard errors must match the order of the ability
#'   estimates provided in the `score` argument. Default is `NULL`.
#' @param n.bin A numeric vector of two positive integers specifying the maximum
#'   and minimum numbers of bins (or intervals) on the ability scale. The first
#'   and second values represent the maximum and minimum numbers of bins,
#'   respectively. Default is `c(80, 10)`. See the **Details** section below for
#'   more information.
#' @param min.binsize A positive integer specifying the minimum number of
#'   examinees required in each bin. To ensure stable statistical estimation,
#'   each bin must contain at least the specified number of examinees from both
#'   the reference and focal groups in order to be included in the calculation
#'   of \eqn{\hat{\beta}}. Bins that do not meet this minimum are excluded from
#'   the computation. Default is 3. See the **Details** section for further
#'   explanation.
#' @param max.del A numeric value specifying the maximum allowable proportion of
#'   examinees that may be excluded from either the reference or focal group
#'   during the binning process. This threshold is used when determining the
#'   number of bins on the ability scale automatically. Default is 0.075. See
#'   the **Details** section for more information.
#' @param weight.group A character string specifying the target ability
#'   distribution used to compute the expected DIF measure \eqn{\hat{\beta}} and
#'   its corresponding standard error. Available options are: `"comb"` for the
#'   combined distribution of both the reference and focal groups, `"foc"` for
#'   the focal group's distribution, and `"ref"` for the reference group's
#'   distribution. Default is `"comb"`. See the **Details** section below for
#'   more information.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for the hypothesis test associated with the CATSIB (*beta*) statistic.
#'   Default is 0.05.
#' @param ... Additional arguments passed to the [irtQ::est_score()] function.
#'
#' @details In the CATSIB procedure (Nandakumar & Roussos, 2004),
#' \eqn{\hat{\theta}^{\ast}}— the expected value of \eqn{\theta} regressed on
#' \eqn{\hat{\theta}}—is a continuous variable. The range of
#' \eqn{\hat{\theta}^{\ast}} is divided into *K* equal-width intervals, and
#' examinees are classified into one of these *K* intervals based on their
#' \eqn{\hat{\theta}^{\ast}} values. Any interval containing fewer than three
#' examinees from either the reference or focal group is excluded from the
#' computation of \eqn{\hat{\beta}}, the DIF effect size, to ensure statistical
#' stability. According to Nandakumar and Roussos (2004), the default minimum
#' bin size is 3, which can be controlled via the `min.binsize` argument.
#'
#' To determine an appropriate number of intervals (*K*), [irtQ::catsib()]
#' automatically decreases *K* from a large starting value (e.g., 80) based on
#' the rule proposed by Nandakumar and Roussos (2004). Specifically, if more
#' than 7.5\% of examinees in either the reference or focal group would be
#' excluded due to small bin sizes, the number of bins is reduced by one and the
#' process is repeated. This continues until the retained examinees in each
#' group comprise at least 92.5\% of the total. However, to prevent having too
#' few bins, they recommended a minimum of *K* = 10. Therefore, the default
#' maximum and minimum number of bins are set to 80 and 10, respectively, via
#' `n.bin`. Likewise, the maximum allowable proportion of excluded examinees is
#' set to 0.075 by default through the `max.del` argument.
#'
#' When it comes to the target ability distribution used to compute
#' \eqn{\hat{\beta}}, Li and Stout (1996) and Nandakumar and Roussos (2004)
#' employed the combined-group target ability distribution, which is the default
#' option in `weight.group`. See Nandakumar and Roussos (2004) for further
#' details about the CATSIB method.
#'
#' Although Nandakumar and Roussos (2004) did not propose a purification
#' procedure for DIF analysis using CATSIB, [irtQ::catsib()] can implement an
#' iterative purification process in a manner similar to that of Lim et al.
#' (2022). Specifically, at each iteration, examinees' latent abilities are
#' recalculated using the purified set of items and the scoring method specified
#' in the `method` argument. The iterative purification process terminates
#' either when no additional DIF items are detected or when the number of
#' iterations reaches the limit set by `max.iter`. See Lim et al. (2022) for
#' more details on the purification procedure.
#'
#' Scoring based on a limited number of items may result in large standard
#' errors, which can negatively affect the effectiveness of DIF detection using
#' the CATSIB procedure. The `min.resp` argument can be used to prevent the use
#' of scores with large standard errors, particularly during the purification
#' process. For example, if `min.resp` is not NULL (e.g., `min.resp = 5`), item
#' responses from examinees whose total number of valid responses is below the
#' specified threshold are treated as missing (i.e., NA). As a result, their
#' ability estimates are also treated as missing and are excluded from the
#' CATSIB statistic computation. If `min.resp = NULL`, a score will be computed
#' for any examinee with at least one valid item response.
#'
#' @return This function returns a list consisting of four elements:
#'
#' \item{no_purify}{A list containing the results of the DIF analysis without
#' applying a purification procedure. This list includes:
#'   \describe{
#'     \item{dif_stat}{A data frame containing the results of the CATSIB
#'     statistics for all evaluated items. The columns include the item ID,
#'     CATSIB (*beta*) statistic, standard error of *beta*, standardized *beta*,
#'     p-value for *beta*, sample size of the reference group, sample size of
#'     the focal group, and total sample size.}
#'     \item{dif_item}{A numeric vector identifying items flagged as potential
#'     DIF items based on the CATSIB statistic.}
#'     \item{contingency}{A list of contingency tables used for computing the
#'     CATSIB statistics for each item.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether a purification procedure was
#' applied.}
#'
#' \item{with_purify}{A list containing the results of the DIF analysis with
#' a purification procedure. This list includes:
#'   \describe{
#'     \item{dif_stat}{A data frame containing the results of the CATSIB
#'     statistics for all evaluated items. The columns include the item ID,
#'     CATSIB (*beta*) statistic, standard error of *beta*, standardized *beta*,
#'     p-value for *beta*, sample size of the reference group, sample size of
#'     the focal group, total sample size, and the iteration number (*n*) in
#'     which the CATSIB statistics were computed.}
#'     \item{dif_item}{A numeric vector identifying items flagged as potential
#'     DIF items based on the CATSIB statistic.}
#'     \item{n.iter}{An integer indicating the total number of iterations
#'     performed during the purification process.}
#'     \item{complete}{A logical value indicating whether the purification
#'     process was completed. If FALSE, the process reached the maximum
#'     number of iterations without full convergence.}
#'     \item{contingency}{A list of contingency tables used for computing the
#'     CATSIB statistics for each item during the purification process.}
#'   }
#' }
#'
#' \item{alpha}{The significance level \eqn{\alpha} used to compute the p-values
#' of the CATSIB statistics.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::rdif()], [irtQ::est_irt], [irtQ::est_item()],
#'   [irtQ::simdat()], [irtQ::shape_df()], [irtQ::est_score()]
#'
#' @references Li, H. H., & Stout, W. (1996). A new procedure for detection of
#'   crossing DIF. *Psychometrika, 61*(4), 647-677.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement*.
#'
#'   Nandakumar, R., & Roussos, L. (2004). Evaluation of the CATSIB DIF
#'   procedure in a pretest setting. *Journal of Educational and Behavioral
#'   Statistics, 29*(2), 177-199.
#'
#'   Shealy, R. T., & Stout, W. F. (1993). A model-based standardization
#'   approach that separates true bias/DIF from group ability differences and
#'   detects test bias/DIF as well as item bias/DIF. *Psychometrika, 58*,
#'   159–194.
#'
#'
#' @examples
#' \donttest{
#' # Load required package
#' library("dplyr")
#'
#' ## Uniform DIF Detection
#' ###############################################
#' # (1) Simulate data with true uniform DIF
#' ###############################################
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select 36 3PLM items that are non-DIF
#' par_nstd <-
#'   bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
#'   dplyr::filter(.data$model == "3PLM") %>%
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>%
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#'
#' # Generate four new items to contain uniform DIF
#' difpar_ref <-
#'   shape_df(
#'     par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = 0.15),
#'     item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
#'   )
#'
#' # Introduce uniform DIF in the focal group by shifting b-parameters
#' difpar_foc <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))
#'
#' # Combine the 4 DIF and 36 non-DIF items for both reference and focal groups
#' # Threfore, the first four items now exhibit uniform DIF
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc <- rbind(difpar_foc, par_nstd)
#'
#' # Generate true theta values
#' set.seed(123)
#' theta_ref <- rnorm(500, 0.0, 1.0)
#' theta_foc <- rnorm(500, 0.0, 1.0)
#'
#' # Simulate response data
#' resp_ref <- simdat(par_ref, theta = theta_ref, D = 1)
#' resp_foc <- simdat(par_foc, theta = theta_foc, D = 1)
#' data <- rbind(resp_ref, resp_foc)
#'
#' ###############################################
#' # (2) Estimate item and ability parameters
#' #     using the aggregated data
#' ###############################################
#'
#' # Estimate item parameters
#' est_mod <- est_irt(data = data, D = 1, model = "3PLM")
#' est_par <- est_mod$par.est
#'
#' # Estimate ability parameters using ML
#' theta_est <- est_score(x = est_par, data = data, method = "ML")
#' score <- theta_est$est.theta
#' se <- theta_est$se.theta
#'
#' ###############################################
#' # (3) Conduct DIF analysis
#' ###############################################
#' # Create a vector of group membership indicators
#' # where '1' indicates the focal group
#' group <- c(rep(0, 500), rep(1, 500))
#'
#' # (a)-1 Compute the CATSIB statistic using provided scores,
#' #       without purification
#' dif_1 <- catsib(
#'   x = NULL, data = data, D = 1, score = score, se = se, group = group, focal.name = 1,
#'   weight.group = "comb", alpha = 0.05, missing = NA, purify = FALSE
#' )
#' print(dif_1)
#'
#' # (a)-2 Compute the CATSIB statistic using provided scores,
#' #       with purification
#' dif_2 <- catsib(
#'   x = est_par, data = data, D = 1, score = score, se = se, group = group, focal.name = 1,
#'   weight.group = "comb", alpha = 0.05, missing = NA, purify = TRUE
#' )
#' print(dif_2)
#' }
#'
#' @import dplyr
#' @importFrom janitor adorn_totals
#' @export
catsib <- function(x = NULL,
                   data,
                   score = NULL,
                   se = NULL,
                   group,
                   focal.name,
                   item.skip = NULL,
                   D = 1,
                   n.bin = c(80, 10),
                   min.binsize = 3,
                   max.del = 0.075,
                   weight.group = c("comb", "foc", "ref"),
                   alpha = 0.05,
                   missing = NA,
                   purify = FALSE,
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
  ## (1) prepare DIF analysis
  ## ----------------------------------
  # when purify = TRUE, the item metadata should be provided.
  # if not, stop.
  if (purify & is.null(x)) {
    stop("To implement a purification process, the item metadata must be provided in the argument 'x'.", call. = FALSE)
  }

  # clean the data frame of the item metadata
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # stop when the model includes any polytomous model
    if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
      stop("The current version only supports dichotomous response data.", call. = FALSE)
    }
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # stop when the model includes any polytomous response data
  if (any(data > 1, na.rm = TRUE)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

  # create an item id
  if (is.null(x)) {
    item.id <- paste0("item.", 1:ncol(data))
  } else {
    item.id <- x$id
  }

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }
    # transform scores to a vector form
    if (is.matrix(se) | is.data.frame(se)) {
      score <- as.numeric(data.matrix(se))
    }
  } else {
    # if min.resp is not NULL, find the examinees who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- Rfast::rowsums(!is.na(data))
      loc_less <- which(n_resp < min.resp & n_resp > 0)
      data[loc_less, ] <- NA
    }
    score_rst <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )
    score <- score_rst$est.theta
    se <- score_rst$se.theta
  }

  # confirm the target ability density function
  weight.group <- match.arg(weight.group)
  weight.group <- switch(weight.group,
                         foc = "foc",
                         ref = "ref",
                         comb = "comb"
  )

  # set the maximum & minimum number of bins (intervals)
  max.bin <- n.bin[1]
  min.bin <- n.bin[2]

  # a) when no purification is set
  # do only one iteration of DIF analysis:
  # compute the beta statistic and its SE for all items using
  # corrected theta scores
  dif_rst <-
    catsib_one(
      data = data, group = group, focal.name = focal.name, item.skip = item.skip,
      score = score, se = se, range = range, item.id = item.id, max.bin = max.bin,
      min.bin = min.bin, min.binsize = min.binsize, max.del = max.del,
      weight.group = weight.group, alpha = alpha
    )

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, dif_item = NULL, contingency = NULL)
  with_purify <- list(
    dif_stat = NULL, dif_item = NULL, n.iter = NULL,
    complete = NULL, contingency = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$contingency <- dif_rst$contingency

  # when purification is used
  if (purify) {
    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and contingency tables
    dif_item <- NULL
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)), beta = NA, se = NA,
        z.beta = NA, p = NA, n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )
    contingency <- vector("list", nrow(x))
    names(contingency) <- item.id

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item
    dif_stat_tmp <- dif_rst$dif_stat
    contingency_tmp <- dif_rst$contingency

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

        # a flagged item which has the largest significant DIF statistic
        flag_max <- which.max(abs(dif_stat_tmp$z.beta))

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:8] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, 9] <- i - 1
        contingency[del_item] <- contingency_tmp[flag_max]

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

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
        score_rst_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )
        score_puri <- score_rst_puri$est.theta
        se_puri <- score_rst_puri$se.theta

        # do DIF analysis using the updated ability estimates
        item.id <- x_puri$id
        dif_rst_tmp <-
          catsib_one(
            data = data_puri, group = group, focal.name = focal.name,
            item.skip = item.skip.puri, score = score_puri, se = se_puri,
            range = range, item.id = item.id, max.bin = max.bin, min.bin = min.bin,
            min.binsize = min.binsize, max.del = max.del, weight.group = weight.group,
            alpha = alpha
          )

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        contingency_tmp <- dif_rst_tmp$contingency

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:8] <- dif_stat_tmp
          dif_stat[item_num, 9] <- i
          contingency[item_num] <- contingency_tmp

          break
        }
      }

      # print a message
      if (verbose) {
        cat("", "\n")
      }

      # record the actual number of iteration
      n_iter <- i

      # if the iteration reached out the maximum number of iteration but the purification was incomplete,
      # then, return a warning message
      if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:8] <- dif_stat_tmp
        dif_stat[item_num, 9] <- i
        contingency[item_num] <- contingency_tmp
      } else {
        complete <- TRUE

        # print a message
        if (verbose) {
          cat("Purification is finished.", "\n")
        }
      }

      # record the final DIF detection results with the purification procedure
      with_purify$dif_stat <- dif_stat
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$complete <- complete
      with_purify$contingency <- contingency
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
      with_purify$contingency <- no_purify$contingency
    }
  }


  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify,
              with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "catsib"
  rst$call <- cl
  rst
}


# This function performs a regression correction for ability estimates and, then
# computes the beta statistic and its SE for all items
catsib_one <- function(data,
                       group,
                       focal.name,
                       item.skip = NULL,
                       score,
                       se,
                       range,
                       item.id,
                       max.bin,
                       min.bin,
                       min.binsize = 3,
                       max.del = 0.075,
                       weight.group,
                       alpha = 0.05) {

  ## ------------------------------------------
  ## prepare data sets
  ## ------------------------------------------
  # count the number of items
  nitem <- length(item.id)

  # find the location of examinees for the reference and the focal groups
  loc_ref <- which(group != focal.name)
  loc_foc <- which(group == focal.name)

  # divide the response data into the two group data
  resp_ref <- data[loc_ref, , drop = FALSE]
  resp_foc <- data[loc_foc, , drop = FALSE]

  # divide the thetas into the two group data
  score_ref <- score[loc_ref]
  score_foc <- score[loc_foc]

  # divide the SEs into the two group data
  se_ref <- se[loc_ref]
  se_foc <- se[loc_foc]

  ## ------------------------------------------
  ## a regression correction for theta scores
  ## ------------------------------------------
  # compute the mean and variance of scores for each group
  # note that, the bound scores are excluded.
  lb_score <- range[1]
  up_score <- range[2]
  mu_ref <- mean(score_ref[score_ref > lb_score & score_ref < up_score], na.rm = TRUE)
  mu_foc <- mean(score_foc[score_foc > lb_score & score_foc < up_score], na.rm = TRUE)
  sigma2_ref <- stats::var(score_ref[score_ref > lb_score & score_ref < up_score], na.rm = TRUE)
  sigma2_foc <- stats::var(score_foc[score_foc > lb_score & score_foc < up_score], na.rm = TRUE)

  # compute the error variance of scores for each group
  errvar_ref <- mean((se_ref^2)[score_ref > lb_score & score_ref < up_score], na.rm = TRUE)
  errvar_foc <- mean((se_foc^2)[score_foc > lb_score & score_foc < up_score], na.rm = TRUE)

  # compute the squared correlation (a.k.a. reliability) between theta estimate and true theta
  rho_ref2 <- suppressWarnings(1 - errvar_ref / sigma2_ref)
  rho_foc2 <- suppressWarnings(1 - errvar_foc / sigma2_foc)

  # apply a regression correction to the ability estimates
  crscore_ref <- mu_ref + rho_ref2 * (score_ref - mu_ref)
  crscore_foc <- mu_foc + rho_foc2 * (score_foc - mu_foc)

  ## ------------------------------------------
  ## compute CATSIB statistic
  ## ------------------------------------------
  # conduct a DIF analysis
  catsib_dif <- purrr::map(
    .x = 1:nitem,
    .f = function(i) {
      catsib_item(
        crscore_ref = crscore_ref, crscore_foc = crscore_foc,
        resp.ref = resp_ref[, i], resp.foc = resp_foc[, i],
        max.bin = max.bin, min.bin = min.bin, min.binsize = min.binsize,
        max.del = max.del, weight.group = weight.group
      )
    }
  )

  # extract the DIF analysis results
  dif_stat <-
    purrr::map(.x = catsib_dif, "dif_stat") %>%
    do.call(what = "rbind") %>%
    dplyr::mutate_if("is.numeric", "round", digits = 4) %>%
    dplyr::mutate(id = item.id) %>%
    dplyr::relocate("id", .before = "beta")
  contingency <- purrr::map(.x = catsib_dif, "contingency")
  names(contingency) <- item.id

  # when there are items that should be skipped for the DIF analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    dif_stat[item.skip, 2:5] <- NA
    contingency[item.skip] <- NA
  }

  # find the flagged items
  dif_item <- as.numeric(which(dif_stat$p <= alpha))
  if (length(dif_item) == 0) dif_item <- NULL

  # return the results
  rst <- list(dif_stat = dif_stat, dif_item = dif_item, contingency = contingency)
  rst
}


# This function computes the beta statistic and its SE for an item
catsib_item <- function(crscore_ref, crscore_foc, resp.ref, resp.foc,
                        max.bin, min.bin, min.binsize = 3, max.del = 0.075,
                        weight.group) {
  # combine all corrected theta scores
  crscore <- c(crscore_ref, crscore_foc)

  # set the range of the ability scale
  min.crscore <- min(crscore, na.rm = TRUE)
  max.crscore <- max(crscore, na.rm = TRUE)

  # decide the number of bins and create an initial frequency table
  for (num.bin in max.bin:min.bin) {
    # compute the cut-scores to divide the theta scale into the bins
    cutscore <- seq(from = min.crscore, to = max.crscore, length.out = num.bin + 1)

    # assign a group variable to each score
    bin_ref <- cut(crscore_ref, breaks = cutscore, include.lowest = TRUE, dig.lab = 7)
    bin_foc <- cut(crscore_foc, breaks = cutscore, include.lowest = TRUE, dig.lab = 7)

    # exclude bins where candidates' responses are NAs
    non.na.ref <- !is.na(resp.ref)
    non.na.foc <- !is.na(resp.foc)
    bin_ref <- bin_ref[non.na.ref]
    bin_foc <- bin_foc[non.na.foc]

    # create a temporary data frame of bin frequency for both groups
    bin.n.ref <- stats::xtabs(~bin_ref, drop.unused.levels = FALSE)
    bin.n.foc <- stats::xtabs(~bin_foc, drop.unused.levels = FALSE)

    # check if the counts of remaining sample is greater than equal to minimum a criterion
    isok_ref <- (sum(bin.n.ref[bin.n.ref >= min.binsize & bin.n.foc >= min.binsize]) /
                   sum(bin.n.ref)) >= 1 - max.del
    isok_foc <- (sum(bin.n.foc[bin.n.ref >= min.binsize & bin.n.foc >= min.binsize]) /
                   sum(bin.n.foc)) >= 1 - max.del

    # if the criterion is met, then break out the loop
    if (all(isok_ref, isok_foc)) {
      break
    }
  }

  # final data frame containing all components to compute the beta statistic
  prop.ref <- as.numeric(table(bin_ref, resp.ref[non.na.ref])[, 2] / bin.n.ref)
  prop.foc <- as.numeric(table(bin_foc, resp.foc[non.na.foc])[, 2] / bin.n.foc)
  var.ref <- as.numeric(by(
    data = resp.ref[non.na.ref],
    INDICES = bin_ref,
    FUN = stats::var, na.rm = TRUE
  ))
  var.foc <- as.numeric(by(
    data = resp.foc[non.na.foc],
    INDICES = bin_foc,
    FUN = stats::var, na.rm = TRUE
  ))
  ref.df <-
    as.data.frame(bin.n.ref, stringsAsFactors = FALSE) %>%
    stats::setNames(nm = c("bin", "n.ref")) %>%
    data.frame(prop.ref = prop.ref, var.ref = var.ref)
  foc.df <-
    as.data.frame(bin.n.foc, stringsAsFactors = FALSE) %>%
    stats::setNames(nm = c("bin", "n.foc")) %>%
    data.frame(prop.foc = prop.foc, var.foc = var.foc)
  n.ref <- n.foc <- weight <- NULL
  item_df <-
    merge(x = ref.df, y = foc.df, by = "bin", all = TRUE, sort = FALSE) %>%
    subset(n.ref >= 3 & n.foc >= 3) %>%
    transform(n.total = n.ref + n.foc) %>%
    dplyr::mutate(
      weight = dplyr::case_when(
        weight.group == "comb" ~ .data$n.total / sum(.data$n.total),
        weight.group == "foc" ~ .data$n.foc / sum(.data$n.foc),
        weight.group == "ref" ~ .data$n.ref / sum(.data$n.ref)
      )
    ) %>%
    transform(
      beta = ((prop.ref - prop.foc) * weight),
      var.beta = ((var.ref / n.ref + var.foc / n.foc) * weight^2)
    )

  # compute the beta statistic and its SE
  beta <- sum(item_df$beta)
  se_beta <- sqrt(sum(item_df$var.beta))
  z_beta <- beta / se_beta
  pval_beta <- 2 * stats::pnorm(q = abs(z_beta), mean = 0, sd = 1, lower.tail = FALSE)
  n.ref <- sum(item_df$n.ref)
  n.foc <- sum(item_df$n.foc)
  stat_df <- data.frame(
    beta = beta, se = se_beta, z.beta = z_beta, p = pval_beta,
    n.ref = n.ref, n.foc = n.foc, n.total = n.ref + n.foc
  )

  # round the numbers of the data frame
  item_df2 <-
    janitor::adorn_totals(
      dat = item_df, where = "row", fill = NA, , , "n.ref", "n.foc", "n.total",
      "beta", "var.beta"
    ) %>%
    dplyr::mutate_at(.vars = c(3, 4, 6, 7), "round", digits = 4)

  # return the results
  list(dif_stat = stat_df, contingency = item_df2)
}
