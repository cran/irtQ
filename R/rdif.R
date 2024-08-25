#' IRT residual-based differential item functioning (RDIF) detection framework
#'
#' @description This function computes three RDIF statistics (Lim & Choe, In press; Lim, Choe, & Han, 2022),
#' which are \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, for each item. \eqn{RDIF_{R}} primarily
#' captures the typical contrast in raw residual pattern between two groups caused by uniform DIF whereas
#' \eqn{RDIF_{S}} primarily captures the typical contrast in squared residual pattern between two groups caused
#' by nonuniform DIF. \eqn{RDIF_{RS}} can reasonably capture both types of DIF.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object of class \code{\link{est_item}}
#' obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}.
#' The item metadata can be easily created using the function \code{\link{shape_df}}. See \code{\link{est_irt}}, \code{\link{irtfit}},
#' \code{\link{info}} or \code{\link{simdat}} for more details about the item metadata.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param score A vector of examinees' ability estimates. If the abilities are not provided, \code{\link{rdif}} estimates the abilities before
#' computing the RDIF statistics. See \code{\link{est_score}} for more details about scoring methods. Default is NULL.
#' @param group A numeric or character vector indicating group membership of examinees. The length of the vector should be the same
#' with the number of rows in the response data matrix.
#' @param focal.name A single numeric or character scalar representing the level associated with the focal group. For instance,
#' given \code{group = c(0, 1, 0, 1, 1)} and '1' indicating the focal group, set \code{focal.name = 1}.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param alpha A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test using the RDIF statistics.
#' Default is .05.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param purify A logical value indicating whether a purification process will be implemented or not. Default is FALSE.
#' @param purify.by A character string specifying a RDIF statistic with which the purification is implemented. Available statistics
#' are "rdifrs" for \eqn{RDIF_{RS}}, "rdifr" for \eqn{RDIF_{R}}, and "rdifs" for \eqn{RDIF_{S}}.
#' @param max.iter A positive integer value specifying the maximum number of iterations for
#' the purification process. Default is 10.
#' @param min.resp A positive integer value specifying the minimum number of item responses for an examinee
#' required to compute the ability estimate. Default is NULL. See details below for more information.
#' @param method A character string indicating a scoring method. Available methods are "ML" for the maximum likelihood estimation,
#' "WL" for the weighted likelihood estimation, "MAP" for the maximum a posteriori estimation, and "EAP" for the expected
#' a posteriori estimation. Default method is "ML".
#' @param range A numeric vector of two components to restrict the range of ability scale for the ML, WL, EAP, and MAP scoring methods.
#' Default is c(-5, 5).
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
#' c(0,1). Ignored if \code{method} is "ML" or "WL".
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
#' Ignored if \code{method} is "ML", "WL", or "MAP".
#' @param weights A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
#' using the function \code{\link{gen.weight}}. If NULL and \code{method} is "EAP", default values are used (see the arguments
#' of \code{norm.prior} and \code{nquad}). Ignored if \code{method} is "ML", "WL" or "MAP".
#' @param ncore The number of logical CPU cores to use. Default is 1. See \code{\link{est_score}} for details.
#' @param verbose A logical value. If TRUE, the progress messages of purification procedure are suppressed. Default is TRUE.
#' @param ... Additional arguments that will be forwarded to the \code{\link{est_score}} function.
#'
#' @details The RDIF framework (Lim et al., 2022) consists of three IRT residual-based statistics: \eqn{RDIF_{R}}, \eqn{RDIF_{S}},
#' and \eqn{RDIF_{RS}}. Under the null hypothesis that a test contains no DIF items, \eqn{RDIF_{R}} and \eqn{RDIF_{S}} follow
#' normal distributions asymptotically. \eqn{RDIF_{RS}} is a based on a bivariate normal distribution of \eqn{RDIF_{R}} and
#' \eqn{RDIF_{S}} statistics. Under the null hypothesis of no DIF items, it follows a \eqn{\chi^{2}} distribution asymptotically
#' with 2 degrees of freedom. See Lim et al. (2022) for more details about RDIF framework.
#'
#' The \code{\link{rdif}} function computes all three RDIF statistics of \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}. The current
#' version of \code{\link{rdif}} function supports both dichotomous and polytomous item response data. To compute the three statistics, the \code{\link{rdif}} function
#' requires (1) item parameter estimates obtained from aggregate data regardless of group membership, (2) examinees' ability estimates
#' (e.g., ML), and (3) examinees' item response data. Note that the ability estimates need to be computed using the aggregate data-based
#' item parameter estimates. The item parameter estimates should be provided in the \code{x} argument, the ability estimates should
#' be provided in the \code{score} argument, and the response data should be provided in the \code{data} argument. When the abilities
#' are not given in the \code{score} argument (i.e., \code{score = NULL}), the \code{\link{rdif}} function estimates examinees' abilities
#' automatically using the scoring method specified in the \code{method} argument (e.g., \code{method = "ML"}).
#'
#' The \code{group} argument accepts a vector of either two distinct numeric or character variables. Between two distinct variable, one is to
#' represent the reference group and another one is to represent the focal group. The length of the vector should be the same with the number
#' of rows in the response data and each value in the vector should indicate each examinee of the response data. Once the \code{gruop} is
#' specified, a single numeric or character value needs to be provided in the \code{focal.name} argument to define which group variable in
#' the \code{group} argument represents the focal group.
#'
#' As other DIF detection approaches, an iterative purification process can be implemented for the RDIF framework.
#' When \code{purify = TRUE}, the purification process is implemented based on one of RDIF statistics specified in the \code{purify.by}
#' argument (e.g, \code{purify.by="rdifrs"}). At each iterative purification, examinees' latent abilities are computed using purified items and
#' scoring method specified in the \code{method} argument. The iterative purification process stops when no further DIF items are found or
#' the process reaches a predetermined limit of iteration, which can be specified in the \code{max.iter} argument. See Lim et al. (2022)
#' for more details about the purification procedure.
#'
#' Scoring with a limited number of items can result in large standard errors, which may impact the effectiveness of DIF detection within
#' the RDIF framework. The \code{min.resp} argument can be employed to avoid using scores with significant standard errors when calculating
#' the RDIF statistics, particularly during the purification process. For instance, if \code{min.resp} is not NULL (e.g., \code{min.resp=5}),
#' item responses from examinees whose total item responses fall below the specified minimum number are treated as missing values (i.e., NA).
#' Consequently, their ability estimates become missing values and are not utilized in computing the RDIF statistics. If \code{min.resp=NULL},
#' an examinee's score will be computed as long as there is at least one item response for the examinee.
#'
#'
#' @return This function returns a list of four internal objects. The four objects are:
#' \item{no_purify}{A list of several sub-objects containing the results of DIF analysis without a purification procedure. The sub-objects are:
#'     \describe{
#'       \item{dif_stat}{A data frame containing the results of three RDIF statistics across all evaluated items. From the first column, each column
#'        indicates item's ID, \eqn{RDIF_{R}} statistic, standardized \eqn{RDIF_{R}}, \eqn{RDIF_{S}} statistic, standardized \eqn{RDIF_{S}},
#'        \eqn{RDIF_{RS}} statistic, p-value of the \eqn{RDIF_{R}}, p-value of the \eqn{RDIF_{S}}, p-value of the \eqn{RDIF_{RS}}, sample size of
#'        the reference group, sample size of the focal group, and total sample size, respectively. Note that \eqn{RDIF_{RS}} does not have its standardized
#'        value because it is a \eqn{\chi^{2}} statistic.}
#'       \item{moments}{A data frame containing the moments of three RDIF statistics. From the first column, each column indicates item's ID,
#'        mean of \eqn{RDIF_{R}}, standard deviation of \eqn{RDIF_{R}}, mean of \eqn{RDIF_{S}}, standard deviation of \eqn{RDIF_{S}}, and
#'        covariance of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}, respectively.}
#'       \item{dif_item}{A list of three numeric vectors showing potential DIF items flagged by each of the RDIF statistics. Each of the numeric vector
#'        means the items flagged by \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
#'       \item{score}{A vector of ability estimates used to compute the RDIF statistics.}
#'    }
#' }
#' \item{purify}{A logical value indicating whether the purification process was used.}
#' \item{with_purify}{A list of several sub-objects containing the results of DIF analysis with a purification procedure. The sub-objects are:
#'     \describe{
#'       \item{purify.by}{A character string indicating which RDIF statistic is used for the purification. "rdifr", "rdifs", and "rdifrs" refers to
#'        \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
#'       \item{dif_stat}{A data frame containing the results of three RDIF statistics across all evaluated items. From the first column, each column
#'        indicates item's ID, \eqn{RDIF_{R}} statistic, standardized \eqn{RDIF_{R}}, \eqn{RDIF_{S}} statistic, standardized \eqn{RDIF_{S}},
#'        \eqn{RDIF_{RS}} statistic, p-value of the \eqn{RDIF_{R}}, p-value of the \eqn{RDIF_{S}}, p-value of the \eqn{RDIF_{RS}}, sample size of
#'        the reference group, sample size of the focal group, total sample size, and \emph{n}th iteration where the RDIF statistics were computed,
#'        respectively.}
#'       \item{moments}{A data frame containing the moments of three RDIF statistics. From the first column, each column indicates item's ID,
#'        mean of \eqn{RDIF_{R}}, standard deviation of \eqn{RDIF_{R}}, mean of \eqn{RDIF_{S}}, standard deviation of \eqn{RDIF_{S}}, covariance
#'        of \eqn{RDIF_{R}} and \eqn{RDIF_{S}}, and \emph{n}th iteration where the RDIF statistics were computed, respectively.}
#'       \item{dif_item}{A list of three numeric vectors showing potential DIF items flagged by each of the RDIF statistics. Each of the numeric vector
#'        means the items flagged by \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}, respectively.}
#'       \item{n.iter}{A total number of iterations implemented for the purification.}
#'       \item{score}{A vector of final purified ability estimates used to compute the RDIF statistics.}
#'       \item{complete}{A logical value indicating whether the purification process was completed. If FALSE, it means that the purification process
#'        reached the maximum iteration number but it was not complete.}
#'     }
#' }
#' \item{alpha}{A significance \eqn{\alpha}-level used to compute the p-values of RDIF statistics.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{est_item}}, \code{\link{info}}, \code{\link{simdat}}, \code{\link{shape_df}},
#' \code{\link{gen.weight}}, \code{\link{est_score}}
#'
#' @references
#' Lim, H., & Choe, E. M. (2023). Detecting differential item functioning in CAT using IRT residual DIF approach.
#' \emph{Journal of Educational Measurement}. \doi{doi:10.1111/jedm.12366}.
#'
#' Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential item functioning detection framework in
#' item response theory. \emph{Journal of Educational Measurement, 59}(1), 80-104. \doi{doi:10.1111/jedm.12313}.
#'
#' @examples
#' \donttest{
#' # call library
#' library("dplyr")
#'
#' ## Uniform DIF detection
#' ###############################################
#' # (1) manipulate true uniform DIF data
#' ###############################################
#' # import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # select 36 of 3PLM items which are non-DIF items
#' par_nstd <-
#'   bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
#'   dplyr::filter(.data$model == "3PLM") %>%
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>%
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#'
#' # generate four new items to inject uniform DIF
#' difpar_ref <-
#'   shape_df(
#'     par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = 0.15),
#'     item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
#'   )
#'
#' # manipulate uniform DIF on the four new items by adding constants to b-parameters
#' # for the focal group
#' difpar_foc <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))
#'
#' # combine the 4 DIF and 36 non-DIF items for both reference and focal groups
#' # thus, the first four items have uniform DIF
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc <- rbind(difpar_foc, par_nstd)
#'
#' # generate the true thetas
#' set.seed(123)
#' theta_ref <- rnorm(500, 0.0, 1.0)
#' theta_foc <- rnorm(500, 0.0, 1.0)
#'
#' # generate the response data
#' resp_ref <- simdat(par_ref, theta = theta_ref, D = 1)
#' resp_foc <- simdat(par_foc, theta = theta_foc, D = 1)
#' data <- rbind(resp_ref, resp_foc)
#'
#' ###############################################
#' # (2) estimate the item and ability parameters
#' #     using the aggregate data
#' ###############################################
#' # estimate the item parameters
#' est_mod <- est_irt(data = data, D = 1, model = "3PLM")
#' est_par <- est_mod$par.est
#'
#' # estimate the ability parameters using ML
#' score <- est_score(x = est_par, data = data, method = "ML")$est.theta
#'
#' ###############################################
#' # (3) conduct DIF analysis
#' ###############################################
#' # create a vector of group membership indicators
#' # where '1' indicates the focal group
#' group <- c(rep(0, 500), rep(1, 500))
#'
#' # (a)-1 compute RDIF statistics by providing scores,
#' #       and without a purification
#' dif_nopuri_1 <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05
#' )
#' print(dif_nopuri_1)
#'
#' # (a)-2 compute RDIF statistics by not providing scores
#' #       and without a purification
#' dif_nopuri_2 <- rdif(
#'   x = est_par, data = data, score = NULL,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   method = "ML"
#' )
#' print(dif_nopuri_2)
#'
#' # (b)-1 compute RDIF statistics with a purification
#' #       based on RDIF(R)
#' dif_puri_r <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "rdifr"
#' )
#' print(dif_puri_r)
#'
#' # (b)-2 compute RDIF statistics with a purification
#' #       based on RDIF(S)
#' dif_puri_s <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "rdifs"
#' )
#' print(dif_puri_s)
#'
#' # (b)-3 compute RDIF statistics with a purification
#' #       based on RDIF(RS)
#' dif_puri_rs <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "rdifrs"
#' )
#' print(dif_puri_rs)
#' }
#'
#' @export
rdif <- function(x, ...) UseMethod("rdif")

#' @describeIn rdif Default method to computes three RDIF statistics using a data frame \code{x} containing the item metadata.
#'
#' @export
rdif.default <- function(x, data, score = NULL, group, focal.name, D = 1, alpha = 0.05, missing = NA, purify = FALSE,
                         purify.by = c("rdifrs", "rdifr", "rdifs"), max.iter = 10, min.resp = NULL, method = "ML",
                         range = c(-5, 5), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1, verbose = TRUE, ...) {
  # match.call
  cl <- match.call()

  ## ----------------------------------
  ## (1) prepare DIF analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  # if(any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
  #   stop("The current version only supports dichotomous response data.", call.=FALSE)
  # }

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
  dif_rst <- rdif_one(x = x, data = data, score = score, group = group, focal.name = focal.name, D = D, alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- data.frame(
    id = x$id,
    dif_rst$moments$rdifr[, c(1, 3)],
    dif_rst$moments$rdifs[, c(1, 3)],
    dif_rst$covariance, stringsAsFactors = FALSE
  )
  names(no_purify$moments) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs", "sigma.rdifs", "covariance")
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
        id = rep(NA_character_, nrow(x)), rdifr = NA, z.rdifr = NA,
        rdifs = NA, z.rdifs = NA, rdifrs = NA, p.rdifr = NA, p.rdifs = NA, p.rdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA, stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.rdifr = NA, sigma.rdifr = NA,
        mu.rdifs = NA, sigma.rdifs = NA, covariance = NA, n.iter = NA, stringsAsFactors = FALSE
      )

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_df_tmp <- no_purify$moments

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
        flag_max <-
          switch(purify.by,
            rdifr = which.max(abs(dif_stat_tmp$z.rdifr)),
            rdifs = which.max(abs(dif_stat_tmp$z.rdifs)),
            rdifrs = which.max(dif_stat_tmp$rdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:12] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, D = D, alpha = alpha
        )

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <- data.frame(
          id = dif_rst_tmp$dif_stat$id,
          dif_rst_tmp$moments$rdifr[, c(1, 3)],
          dif_rst_tmp$moments$rdifs[, c(1, 3)],
          dif_rst_tmp$covariance, stringsAsFactors = FALSE
        )
        names(mmt_df_tmp) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs", "sigma.rdifs", "covariance")

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:12] <- dif_stat_tmp
          dif_stat[item_num, 13] <- i
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
      if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:12] <- dif_stat_tmp
        dif_stat[item_num, 13] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
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
      with_purify$moments <- mmt_df
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "rdif"
  rst$call <- cl
  rst
}

#' @describeIn rdif An object created by the function \code{\link{est_irt}}.
#'
#' @export
#'
rdif.est_irt <- function(x, score = NULL, group, focal.name, alpha = 0.05, missing = NA, purify = FALSE,
                         purify.by = c("rdifrs", "rdifr", "rdifs"), max.iter = 10, min.resp = NULL, method = "ML",
                         range = c(-5, 5), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1, verbose = TRUE, ...) {
  # match.call
  cl <- match.call()

  # extract information from an object
  data <- x$data
  D <- x$scale.D
  x <- x$par.est

  ## ----------------------------------
  ## (1) prepare DIF analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  # if(any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
  #   stop("The current version only supports dichotomous response data.", call.=FALSE)
  # }

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
  dif_rst <- rdif_one(x = x, data = data, score = score, group = group, focal.name = focal.name, D = D, alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- data.frame(
    id = x$id,
    dif_rst$moments$rdifr[, c(1, 3)],
    dif_rst$moments$rdifs[, c(1, 3)],
    dif_rst$covariance, stringsAsFactors = FALSE
  )
  names(no_purify$moments) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs", "sigma.rdifs", "covariance")
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
        id = rep(NA_character_, nrow(x)), rdifr = NA, z.rdifr = NA,
        rdifs = NA, z.rdifs = NA, rdifrs = NA, p.rdifr = NA, p.rdifs = NA, p.rdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA, stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.rdifr = NA, sigma.rdifr = NA,
        mu.rdifs = NA, sigma.rdifs = NA, covariance = NA, n.iter = NA, stringsAsFactors = FALSE
      )

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_df_tmp <- no_purify$moments

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
        flag_max <-
          switch(purify.by,
            rdifr = which.max(abs(dif_stat_tmp$z.rdifr)),
            rdifs = which.max(abs(dif_stat_tmp$z.rdifs)),
            rdifrs = which.max(dif_stat_tmp$rdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:12] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, D = D, alpha = alpha
        )

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <- data.frame(
          id = dif_rst_tmp$dif_stat$id,
          dif_rst_tmp$moments$rdifr[, c(1, 3)],
          dif_rst_tmp$moments$rdifs[, c(1, 3)],
          dif_rst_tmp$covariance, stringsAsFactors = FALSE
        )
        names(mmt_df_tmp) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs", "sigma.rdifs", "covariance")

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:12] <- dif_stat_tmp
          dif_stat[item_num, 13] <- i
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
      if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:12] <- dif_stat_tmp
        dif_stat[item_num, 13] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
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
      with_purify$moments <- mmt_df
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "rdif"
  rst$call <- cl
  rst
}


#' @describeIn rdif An object created by the function \code{\link{est_item}}.
#'
#' @export
#'
rdif.est_item <- function(x, group, focal.name, alpha = 0.05, missing = NA, purify = FALSE,
                          purify.by = c("rdifrs", "rdifr", "rdifs"), max.iter = 10, min.resp = NULL, method = "ML",
                          range = c(-5, 5), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1, verbose = TRUE, ...) {
  # match.call
  cl <- match.call()

  # extract information from an object
  data <- x$data
  score <- x$score
  D <- x$scale.D
  x <- x$par.est

  ## ----------------------------------
  ## (1) prepare DIF analysis
  ## ----------------------------------
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  # if(any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
  #   stop("The current version only supports dichotomous response data.", call.=FALSE)
  # }

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
  dif_rst <- rdif_one(x = x, data = data, score = score, group = group, focal.name = focal.name, D = D, alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL)
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <- data.frame(
    id = x$id,
    dif_rst$moments$rdifr[, c(1, 3)],
    dif_rst$moments$rdifs[, c(1, 3)],
    dif_rst$covariance, stringsAsFactors = FALSE
  )
  names(no_purify$moments) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs", "sigma.rdifs", "covariance")
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
        id = rep(NA_character_, nrow(x)), rdifr = NA, z.rdifr = NA,
        rdifs = NA, z.rdifs = NA, rdifrs = NA, p.rdifr = NA, p.rdifs = NA, p.rdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA, stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.rdifr = NA, sigma.rdifr = NA,
        mu.rdifs = NA, sigma.rdifs = NA, covariance = NA, n.iter = NA, stringsAsFactors = FALSE
      )

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_df_tmp <- no_purify$moments

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
        flag_max <-
          switch(purify.by,
            rdifr = which.max(abs(dif_stat_tmp$z.rdifr)),
            rdifs = which.max(abs(dif_stat_tmp$z.rdifs)),
            rdifrs = which.max(dif_stat_tmp$rdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:12] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, 13] <- i - 1
        mmt_df[del_item, 1:6] <- mmt_df_tmp[flag_max, ]
        mmt_df[del_item, 7] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- rowSums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, D = D, alpha = alpha
        )

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <- data.frame(
          id = dif_rst_tmp$dif_stat$id,
          dif_rst_tmp$moments$rdifr[, c(1, 3)],
          dif_rst_tmp$moments$rdifs[, c(1, 3)],
          dif_rst_tmp$covariance, stringsAsFactors = FALSE
        )
        names(mmt_df_tmp) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs", "sigma.rdifs", "covariance")

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:12] <- dif_stat_tmp
          dif_stat[item_num, 13] <- i
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
      if (max.iter == n_iter & !is.null(dif_item_tmp)) {
        warning("The iteration reached out the maximum number of iteration before purification is completed.", call. = FALSE)
        complete <- FALSE

        # add flagged DIF item at the last iteration
        dif_item <- c(dif_item, item_num[dif_item_tmp])

        # add the DIF statistics for rest of items
        dif_stat[item_num, 1:12] <- dif_stat_tmp
        dif_stat[item_num, 13] <- i
        mmt_df[item_num, 1:6] <- mmt_df_tmp
        mmt_df[item_num, 7] <- i
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
      with_purify$moments <- mmt_df
      with_purify$dif_item <- sort(dif_item)
      with_purify$n.iter <- n_iter
      with_purify$score <- score_puri
      with_purify$complete <- complete
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <- cbind(no_purify$moments, n.iter = 0)
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "rdif"
  rst$call <- cl
  rst
}


# This function conducts one iteration of DIF analysis using the IRT residual based statistics
#' @import dplyr
rdif_one <- function(x, data, score, group, focal.name, D = 1, alpha = 0.05) {
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

  # compute the residual-based DIF statistic
  rdifr <- colMeans(resid_foc, na.rm = TRUE) - colMeans(resid_ref, na.rm = TRUE) # the difference of mean raw residuals (DMRR)
  rdifs <- colMeans(resid_foc^2, na.rm = TRUE) - colMeans(resid_ref^2, na.rm = TRUE) # the difference of mean squared residuals (DMSR)

  # compute the means and variances of the two statistics for the hypothesis testing
  moments <- resid_moments(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats = cats
  )
  moments_rdifr <- moments$rdifr
  moments_rdifs <- moments$rdifs
  covar <- moments$covariance
  poolsd <- moments$poolsd

  # compute the chi-square statistics
  chisq <- c()
  for (i in 1:nitem) {
    if (i %in% all_miss) {
      chisq[i] <- NaN
    } else {
      # create a var-covariance matrix between rdifr and rdifs
      cov_mat <- array(NA, c(2, 2))

      # replace NAs with the analytically computed covariance
      cov_mat[col(cov_mat) != row(cov_mat)] <- covar[i]

      # replace NAs with the analytically computed variances
      diag(cov_mat) <- c(moments_rdifr[i, 2], moments_rdifs[i, 2])

      # create a vector of mean rdifr and mean rdifs
      mu_vec <- cbind(moments_rdifr[i, 1], moments_rdifs[i, 1])

      # create a vector of rdifr and rdifs
      est_mu_vec <- cbind(rdifr[i], rdifs[i])

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

  # standardize the two statistics of rdifr and rdifs
  z_stat_rdifr <- (rdifr - moments_rdifr$mu) / moments_rdifr$sigma
  z_stat_rdifs <- (rdifs - moments_rdifs$mu) / moments_rdifs$sigma

  # calculate p-values for all three statistics
  p_rdifr <- round(2 * stats::pnorm(q = abs(z_stat_rdifr), mean = 0, sd = 1, lower.tail = FALSE), 4)
  p_rdifs <- round(2 * stats::pnorm(q = abs(z_stat_rdifs), mean = 0, sd = 1, lower.tail = FALSE), 4)
  p_rdifrs <- round(stats::pchisq(chisq, df = 2, lower.tail = FALSE), 4)

  # compute three effect size for rdifr and rdifs
  rdif_stats <- list(x = rdifr, y = rdifs)
  n_total <- n_foc + n_ref

  # Cohen's d
  # efs_cohen <-
  #   purrr::map2(.x=rdif_stats, .y=poolsd$cohen, .f=function(x, y) round(x / y, 4)) %>%
  #   data.frame() %>%
  #   dplyr::rename_all(.f=function(x) c("co_rdifr", "co_rdifs"))

  # Hedge's g
  efs_hedge <-
    purrr::map2(
      .x = rdif_stats, .y = poolsd$hedge,
      .f = function(x, y) {
        round((x / y) * (1 - (3 / (4 * (n_foc + n_ref) - 9))), 4)
      }
    ) %>%
    data.frame() %>%
    dplyr::rename_all(.f = function(x) c("he_rdifr", "he_rdifs"))

  # Glass's delta
  efs_glass <-
    purrr::map2(.x = rdif_stats, .y = poolsd$glass, .f = function(x, y) round(x / y, 4)) %>%
    data.frame() %>%
    dplyr::rename_all(.f = function(x) c("gl_rdifr", "gl_rdifs"))

  # combine all effect size
  effect_size <- data.frame(id = x$id, efs_hedge, efs_glass, stringsAsFactors = FALSE)

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id, rdifr = round(rdifr, 4), z.rdifr = round(z_stat_rdifr, 4),
      rdifs = round(rdifs, 4), z.rdifs = round(z_stat_rdifs, 4), rdifrs = round(chisq, 4),
      p.rdifr = p_rdifr, p.rdifs = p_rdifs, p.rdifrs = p_rdifrs,
      n.ref = n_ref, n.foc = n_foc, n.total = n_total, stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # find the flagged items
  dif_item_rdifr <- as.numeric(which(p_rdifr <= alpha))
  dif_item_rdifs <- as.numeric(which(p_rdifs <= alpha))
  dif_item_rdifrs <- which(p_rdifrs <= alpha)
  if (length(dif_item_rdifr) == 0) dif_item_rdifr <- NULL
  if (length(dif_item_rdifs) == 0) dif_item_rdifs <- NULL
  if (length(dif_item_rdifrs) == 0) dif_item_rdifrs <- NULL

  # summarize the results
  rst <- list(
    dif_stat = stat_df, effect_size = effect_size,
    dif_item = list(rdifr = dif_item_rdifr, rdifs = dif_item_rdifs, rdifrs = dif_item_rdifrs),
    moments = list(rdifr = moments_rdifr, rdifs = moments_rdifs), covariance = covar, alpha = alpha
  )

  # return the results
  rst
}


# This function computes the mean and variance of the IRT based residual statistics
# and it computes the covariance between rdifr and rdifs
resid_moments <- function(p_ref, p_foc, n_ref, n_foc, resp_ref, resp_foc, cats) {
  # check the number of items
  nitem <- length(cats)

  # check the number of rows and columns in the two groups response data
  nrow_ref <- nrow(resp_ref)
  nrow_foc <- nrow(resp_foc)

  # create an empty list of two matrices to contain for the first and second moments of raw residual and squared residuals
  mu_ref <- purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = nrow_ref, ncol = nitem))
  mu_foc <- purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = nrow_foc, ncol = nitem))
  mu2_ref <- mu_ref
  mu2_foc <- mu_foc

  # create an empty matrix to contain the expectation of cube of raw residuals at each theta for all items
  mu_resid3_ref <- matrix(NA, nrow = nrow_ref, ncol = nitem)
  mu_resid3_foc <- matrix(NA, nrow = nrow_foc, ncol = nitem)

  # compute the first and second moments of raw residual and squared residuals for all items
  for (i in 1:nitem) {
    # compute the expected residuals for each score category
    Emat_ref <- matrix(0:(cats[i] - 1), nrow = nrow_ref, ncol = cats[i], byrow = TRUE)
    Emat_foc <- matrix(0:(cats[i] - 1), nrow = nrow_foc, ncol = cats[i], byrow = TRUE)
    exp_resid_ref <- Emat_ref - matrix(rowSums(Emat_ref * p_ref[[i]]), nrow = nrow_ref, ncol = cats[i], byrow = FALSE)
    exp_resid_foc <- Emat_foc - matrix(rowSums(Emat_foc * p_foc[[i]]), nrow = nrow_foc, ncol = cats[i], byrow = FALSE)

    # replace NA values into the missing data location
    exp_resid_ref[is.na(resp_ref[, i]), ] <- NA
    exp_resid_foc[is.na(resp_foc[, i]), ] <- NA

    # compute the expected values of raw and squared residuals to be used to compute the mean and variance parameters
    value_foc <- value_ref <- vector("list", 2)
    value_ref[[1]] <- exp_resid_ref
    value_foc[[1]] <- exp_resid_foc
    value_ref[[2]] <- exp_resid_ref^2
    value_foc[[2]] <- exp_resid_foc^2
    resid3_ref <- exp_resid_ref^3
    resid3_foc <- exp_resid_foc^3

    # compute the first and second moments across all thetas for the reference and focal groups
    for (j in 1:2) {
      mu_ref[[j]][, i] <- Rfast::rowsums(value_ref[[j]] * p_ref[[i]], na.rm = FALSE)
      mu2_ref[[j]][, i] <- Rfast::rowsums((value_ref[[j]])^2 * p_ref[[i]], na.rm = FALSE)
      mu_foc[[j]][, i] <- Rfast::rowsums(value_foc[[j]] * p_foc[[i]], na.rm = FALSE)
      mu2_foc[[j]][, i] <- Rfast::rowsums((value_foc[[j]])^2 * p_foc[[i]], na.rm = FALSE)
    }

    # compute the expectation (first moment) of the cube of raw residuals across all thetas for the reference and focal groups
    mu_resid3_ref[, i] <- Rfast::rowsums(resid3_ref * p_ref[[i]], na.rm = FALSE)
    mu_resid3_foc[, i] <- Rfast::rowsums(resid3_foc * p_foc[[i]], na.rm = FALSE)
  }

  # compute the variances across all thetas for the reference and focal groups
  var_ref <- purrr::map2(
    .x = mu2_ref, .y = mu_ref,
    .f = function(x, y) {
      x - y^2
    }
  )
  var_foc <- purrr::map2(
    .x = mu2_foc, .y = mu_foc,
    .f = function(x, y) {
      x - y^2
    }
  )

  # compute the sum of variances for each group
  # use V(aX - bY) = a^2 * V(X) + b^2 * V(Y)
  const_ref <- purrr::map(
    .x = var_ref,
    .f = function(x) {
      matrix((1 / n_ref^2), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    }
  )
  const_foc <- purrr::map(
    .x = var_foc,
    .f = function(x) {
      matrix((1 / n_foc^2), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    }
  )
  var_mu_ref <- purrr::map2(
    .x = const_ref, .y = var_ref,
    .f = function(x, y) {
      Rfast::colsums(x * y, na.rm = TRUE)
    }
  )
  var_mu_foc <- purrr::map2(
    .x = const_foc, .y = var_foc,
    .f = function(x, y) {
      Rfast::colsums(x * y, na.rm = TRUE)
    }
  )

  # compute the final mean and variance
  # use E(X - Y) = E(X) - E(Y)
  # use V(X - Y) = V(X) + V(Y)
  mu <- purrr::map2(
    .x = mu_foc, .y = mu_ref,
    .f = function(x, y) {
      colMeans(x, na.rm = TRUE) - colMeans(y, na.rm = TRUE)
    }
  )
  mu[[1]] <- round(mu[[1]], digits = 10)
  sigma2 <- purrr::map2(.x = var_mu_foc, .y = var_mu_ref, .f = function(x, y) x + y)
  sigma <- purrr::map(.x = sigma2, .f = function(x) sqrt(x))

  # compute two pooled SDs to compute the three effect size: Cohen's d, Hedge's g, Glass's delta
  # poolsd_cohen <-
  #   suppressWarnings(purrr::map2(.x=var_mu_foc, .y=var_mu_ref, .f=function(x, y) sqrt((x + y) / 2)))
  poolsd_hedge <-
    suppressWarnings(purrr::map2(
      .x = var_mu_foc, .y = var_mu_ref,
      .f = function(x, y) {
        sqrt((x * (n_foc - 1) + y * (n_ref - 1)) / (n_foc + n_foc - 2))
      }
    ))
  poolsd_glass <-
    suppressWarnings(purrr::map(.x = var_mu_ref, .f = function(x) sqrt(x)))
  poolsd <- list(hedge = poolsd_hedge, glass = poolsd_glass)

  # compute the covariance between rdifr and rdifs
  covar <- Rfast::colsums(mu_resid3_foc, na.rm = TRUE) * (1 / n_foc^2) + Rfast::colsums(mu_resid3_ref, na.rm = TRUE) * (1 / n_ref^2)

  # return the results
  rst <- purrr::pmap(
    .l = list(x = mu, y = sigma2, z = sigma),
    .f = function(x, y, z) data.frame(mu = x, sigma2 = y, sigma = z)
  )
  names(rst) <- c("rdifr", "rdifs")
  rst$covariance <- covar
  rst$poolsd <- poolsd
  rst
}
