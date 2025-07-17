#' IRT Residual-Based Differential Item Functioning (RDIF) Detection Framework
#'
#' This function computes three RDIF statistics for each item: \eqn{RDIF_{R}},
#' \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}} (Lim & Choe, 2023; Lim, et al., 2022).
#' \eqn{RDIF_{R}} primarily captures differences in raw residuals between two
#' groups, which are typically associated with uniform DIF. \eqn{RDIF_{S}}
#' primarily captures differences in squared residuals, which are typically
#' associated with nonuniform DIF. \eqn{RDIF_{RS}} jointly considers both types
#' of differences and is capable of detecting both uniform and nonuniform DIF.
#'
#' @inheritParams irtfit
#' @inheritParams est_score
#' @param score A numeric vector containing examinees' ability estimates (theta
#'   values). If not provided, [irtQ::rdif()] will estimate ability parameters
#'   internally before computing the RDIF statistics. See [irtQ::est_score()]
#'   for more information on scoring methods. Default is `NULL`.
#' @param group A numeric or character vector indicating examinees' group
#'   membership. The length of the vector must match the number of rows in the
#'   response data matrix.
#' @param focal.name A single numeric or character value specifying the focal
#'   group. For instance, given `group = c(0, 1, 0, 1, 1)` and '1' indicating
#'   the focal group, set `focal.name = 1`.
#' @param item.skip A numeric vector of item indices to exclude from DIF analysis.
#'  If `NULL`, all items are included. Useful for omitting specific items based on
#'  prior insights.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for hypothesis testing using the RDIF statistics. Default is `0.05`.
#' @param missing  A value indicating missing responses in the data set. Default
#'   is `NA`.
#' @param purify Logical. Indicates whether to apply a purification procedure.
#'   Default is `FALSE`.
#' @param purify.by A character string specifying which RDIF statistic is used
#'   to perform the purification. Available options are "rdifrs" for
#'   \eqn{RDIF_{RS}}, "rdifr" for \eqn{RDIF_{R}}, and "rdifs" for
#'   \eqn{RDIF_{S}}.
#' @param max.iter A positive integer specifying the maximum number of
#'   iterations allowed for the purification process. Default is `10`.
#' @param min.resp A positive integer specifying the minimum number of valid
#'   item responses required from an examinee in order to compute an ability
#'   estimate. Default is `NULL`. See **Details** for more information.
#' @param method A character string indicating the scoring method to use.
#'   Available options are:
#'   - `"ML"`: Maximum likelihood estimation
#'   - `"WL"`: Weighted likelihood estimation (Warm, 1989)
#'   - `"MAP"`: Maximum a posteriori estimation (Hambleton et al., 1991)
#'   - `"EAP"`: Expected a posteriori estimation (Bock & Mislevy, 1982)
#'
#'   Default is `"ML"`.
#' @param range A numeric vector of length two specifying the lower and upper
#'   bounds of the ability scale. This is used for the following scoring
#'   methods: `"ML"`, `"WL"`, and `"MAP"`. Default is `c(-5, 5)`.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal prior distribution. These values are used
#'   to generate the Gaussian quadrature points and weights. Ignored if `method`
#'   is `"ML"` or `"WL"`. Default is `c(0, 1)`.
#' @param nquad An integer indicating the number of Gaussian quadrature points
#'   to be generated from the normal prior distribution. Used only when `method`
#'   is `"EAP"`. Ignored for `"ML"`, `"WL"`, and `"MAP"`. Default is 41.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and their corresponding weights (in the second
#'   column) for the latent variable prior distribution. The weights and points
#'   can be conveniently generated using the function [irtQ::gen.weight()].
#'
#'   If `NULL` and `method = "EAP"`, default quadrature values are generated
#'   based on the `norm.prior` and `nquad` arguments. Ignored if `method` is
#'   `"ML"`, `"WL"`, or `"MAP"`.
#' @param ncore An integer specifying the number of logical CPU cores to use for
#'   parallel processing. Default is `1`. See [irtQ::est_score()] for details.
#' @param verbose Logical. If `TRUE`, progress messages from the purification
#'   procedure will be displayed; if `FALSE`, the messages will be suppressed.
#'   Default is `TRUE`.
#' @param ... Additional arguments passed to the [irtQ::est_score()] function.
#'
#' @details The RDIF framework (Lim & Choe, 2023; Lim et al., 2022) consists of
#'   three IRT residual-based statistics: \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and
#'   \eqn{RDIF_{RS}}. Under the null hypothesis that a test contains no DIF
#'   items, \eqn{RDIF_{R}} and \eqn{RDIF_{S}} asymptotically follow standard
#'   normal distributions. \eqn{RDIF_{RS}} is based on a bivariate normal
#'   distribution of the \eqn{RDIF_{R}} and \eqn{RDIF_{S}} statistics, and under
#'   the null hypothesis, it asymptotically follows a \eqn{\chi^{2}}
#'   distribution with 2 degrees of freedom. See Lim et al. (2022) for more
#'   details about the RDIF framework.
#'
#'   The [irtQ::rdif()] function computes all three RDIF statistics:
#'   \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}}. The current version of
#'   [irtQ::rdif()] supports both dichotomous and polytomous item response data.
#'   Note that for polytomous items, net DIF are assessed. To evaluate global
#'   DIF for polytomous items, use [irtQ::crdif()] function.
#'
#'   To compute the RDIF statistics, the [irtQ::rdif()] function requires:
#'   (1) item parameter estimates obtained from aggregate data (regardless
#'   of group membership), (2) examinees' ability estimates (e.g., ML), and
#'   (3) examinees' item response data. Note that the ability estimates must
#'   be based on the aggregate-data item parameters. The item parameter estimates
#'   should be provided in the `x` argument, the ability estimates in the `score`
#'   argument, and the response data in the `data` argument. If ability
#'   estimates are not provided (i.e., `score = NULL`), [irtQ::rdif()] will
#'   estimate them automatically using the scoring method specified via the
#'   `method` argument (e.g., `method = "ML"`).
#'
#'   The `group` argument should be a vector containing exactly two distinct
#'   values (either numeric or character), representing the reference and focal
#'   groups. Its length must match the number of rows in the response data,
#'   where each element corresponds to an examinee. Once `group` is specified, a
#'   single numeric or character value must be provided in the `focal.name`
#'   argument to indicate which level in `group` represents the focal group.
#'
#'   Similar to other DIF detection approaches, the RDIF framework supports an
#'   iterative purification process. When `purify = TRUE`, purification is
#'   conducted using one of the RDIF statistics specified in the `purify.by`
#'   argument (e.g., `purify.by = "rdifrs"`). At each iteration, examinees'
#'   ability estimates are recalculated based on the set of purified items using
#'   the scoring method specified in the `method` argument. The purification
#'   process continues until no additional DIF items are identified or the
#'   maximum number of iterations specified in `max.iter` is reached. See Lim et
#'   al. (2022) for more details on the purification procedure.
#'
#'   Scoring based on a small number of item responses can lead to large
#'   standard errors, potentially reducing the accuracy of DIF detection in the
#'   RDIF framework. The `min.resp` argument can be used to exclude examinees
#'   with insufficient response data from scoring, especially during the
#'   purification process. For example, if `min.resp` is not NULL (e.g.,
#'   `min.resp = 5`), examinees who responded to fewer than five items will have
#'   all their responses treated as missing (i.e., NA). As a result, their
#'   ability estimates will also be missing and will not be used in the
#'   computation of RDIF statistics. If `min.resp = NULL`, a score will be
#'   computed for any examinee with at least one valid item response.
#'
#' @return This function returns a list containing four main components:
#'
#' \item{no_purify}{A list of sub-objects containing the results of DIF analysis
#' without applying a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{dif_stat}{A data frame summarizing the RDIF analysis results for all
#'    items. The columns include: item ID, \eqn{RDIF_{R}} statistic, standardized
#'    \eqn{RDIF_{R}}, \eqn{RDIF_{S}} statistic, standardized \eqn{RDIF_{S}},
#'    \eqn{RDIF_{RS}} statistic, p-values for \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and
#'    \eqn{RDIF_{RS}}, sample sizes for the reference and focal groups, and total
#'    sample size. Note that \eqn{RDIF_{RS}} does not have a standardized value
#'    because it is a \eqn{\chi^{2}}-based statistic.}
#'     \item{moments}{A data frame reporting the first and second moments of the
#'     RDIF statistics. The columns include: item ID, mean and standard deviation
#'     of \eqn{RDIF_{R}}, mean and standard deviation of \eqn{RDIF_{S}}, and the
#'     covariance between \eqn{RDIF_{R}} and \eqn{RDIF_{S}}.}
#'     \item{dif_item}{A list of three numeric vectors identifying items flagged
#'     as DIF by each RDIF statistic: \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and
#'     \eqn{RDIF_{RS}}.}
#'     \item{score}{A numeric vector of ability estimates used to compute the RDIF
#'     statistics.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether the purification procedure
#' was applied.}
#'
#' \item{with_purify}{A list of sub-objects containing the results of DIF analysis
#' with a purification procedure. The sub-objects include:
#'   \describe{
#'     \item{purify.by}{A character string indicating the RDIF statistic used for
#'     purification. Possible values are "rdifr", "rdifs", and "rdifrs",
#'     corresponding to \eqn{RDIF_{R}}, \eqn{RDIF_{S}}, and \eqn{RDIF_{RS}},
#'     respectively.}
#'     \item{dif_stat}{A data frame reporting the RDIF analysis results for
#'     all items across the final iteration. Same structure as in \code{no_purify},
#'     with one additional column indicating the iteration number in which each
#'     result was obtained.}
#'     \item{moments}{A data frame reporting the moments of RDIF statistics
#'     across the final iteration. Includes the same columns as in
#'     \code{no_purify}, with an additional column for the iteration number.}
#'     \item{dif_item}{A list of three numeric vectors identifying DIF items
#'     flagged by each RDIF statistic.}
#'     \item{n.iter}{An integer indicating the total number of iterations
#'     performed during the purification process.}
#'     \item{score}{A numeric vector of purified ability estimates used to
#'     compute the final RDIF statistics.}
#'     \item{complete}{A logical value indicating whether the purification
#'     process converged. If FALSE, the maximum number of iterations was reached
#'     without convergence.}
#'   }
#' }
#'
#' \item{alpha}{A numeric value indicating the significance level (\eqn{\alpha})
#' used in hypothesis testing for RDIF statistics.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::est_item()], [irtQ::simdat()],
#'   [irtQ::shape_df()], [irtQ::est_score()]
#'
#' @references Lim, H., & Choe, E. M. (2023). Detecting differential item
#'   functioning in CAT using IRT residual DIF approach.
#'  *Journal of Educational Measurement, 60*(4), 626-650. \doi{doi:10.1111/jedm.12366}.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{doi:10.1111/jedm.12313}.
#'
#' @examples
#' \donttest{
#' # Load required package
#' library("dplyr")
#'
#' ## Uniform DIF detection
#' ###############################################
#' # (1) Generate data with known uniform DIF
#' ###############################################
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select 36 non-DIF items using the 3PLM model
#' par_nstd <-
#'   bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
#'   dplyr::filter(.data$model == "3PLM") %>%
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>%
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#'
#' # Generate 4 new DIF items for the reference group
#' difpar_ref <-
#'   shape_df(
#'     par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = 0.15),
#'     item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
#'   )
#'
#' # Add uniform DIF by shifting the b-parameters for the focal group
#' difpar_foc <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + rep(0.7, 4))
#'
#' # Combine the DIF and non-DIF items for both reference and focal groups
#' # Therefor, the first 4 items exhibit uniform DIF
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc <- rbind(difpar_foc, par_nstd)
#'
#' # Generate true ability values
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
#' #     from the combined response data
#' ###############################################
#'
#' # Estimate item parameters
#' est_mod <- est_irt(data = data, D = 1, model = "3PLM")
#' est_par <- est_mod$par.est
#'
#' # Estimate ability parameters using ML
#' score <- est_score(x = est_par, data = data, method = "ML")$est.theta
#'
#' ###############################################
#' # (3) Perform DIF analysis
#' ###############################################
#'
#' # Define group membership: 1 = focal group
#' group <- c(rep(0, 500), rep(1, 500))
#'
#' # (a)-1 Compute RDIF statistics with provided ability scores
#' #       (no purification)
#' dif_nopuri_1 <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05
#' )
#' print(dif_nopuri_1)
#'
#' # (a)-2 Compute RDIF statistics without providing ability scores
#' #       (no purification)
#' dif_nopuri_2 <- rdif(
#'   x = est_par, data = data, score = NULL,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   method = "ML"
#' )
#' print(dif_nopuri_2)
#'
#' # (b)-1 Compute RDIF statistics with purification based on RDIF(R)
#' dif_puri_r <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "rdifr"
#' )
#' print(dif_puri_r)
#'
#' # (b)-2 Compute RDIF statistics with purification based on RDIF(S)
#' dif_puri_s <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "rdifs"
#' )
#' print(dif_puri_s)
#'
#' # (b)-3 Compute RDIF statistics with purification based on RDIF(RS)
#' dif_puri_rs <- rdif(
#'   x = est_par, data = data, score = score,
#'   group = group, focal.name = 1, D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "rdifrs"
#' )
#' print(dif_puri_rs)
#'
#' }
#'
#' @export
rdif <- function(x, ...) UseMethod("rdif")

#' @describeIn rdif Default method for computing the three RDIF statistics using
#' a data frame `x` that contains item metadata
#'
#' @export
rdif.default <- function(x,
                         data,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         D = 1,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("rdifrs", "rdifr", "rdifs"),
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
      nquad = nquad, weights = weights, ncore = ncore, ...)$est.theta
  }

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <-
    rdif_one(x = x, data = data, score = score, group = group,
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
  no_purify$moments <- data.frame(
    id = x$id,
    dif_rst$moments$rdifr[, c(1, 3)],
    dif_rst$moments$rdifs[, c(1, 3)],
    dif_rst$covariance, stringsAsFactors = FALSE
  )
  names(no_purify$moments) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs",
                                "sigma.rdifs", "covariance")
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
        rdifr = NA, z.rdifr = NA,
        rdifs = NA, z.rdifs = NA,
        rdifrs = NA, p.rdifr = NA,
        p.rdifs = NA, p.rdifrs = NA,
        n.ref = NA, n.foc = NA, n.total = NA, n.iter = NA,
        stringsAsFactors = FALSE
      )
    mmt_df <-
      data.frame(
        id = rep(NA_character_, nrow(x)), mu.rdifr = NA, sigma.rdifr = NA,
        mu.rdifs = NA, sigma.rdifs = NA, covariance = NA, n.iter = NA,
        stringsAsFactors = FALSE
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
        score_puri <-
          est_score(
            x = x_puri, data = data_puri, D = D, method = method,
            range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
            ncore = ncore, ...)$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
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
        names(mmt_df_tmp) <- c("id", "mu.rdifr", "sigma.rdifr", "mu.rdifs",
                               "sigma.rdifs", "covariance")

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

#' @describeIn rdif An object created by the function [irtQ::est_irt()].
#'
#' @export
#'
rdif.est_irt <- function(x,
                         score = NULL,
                         group,
                         focal.name,
                         item.skip = NULL,
                         alpha = 0.05,
                         missing = NA,
                         purify = FALSE,
                         purify.by = c("rdifrs", "rdifr", "rdifs"),
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
  dif_rst <-
    rdif_one(x = x, data = data, score = score, group = group,
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
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
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


#' @describeIn rdif An object created by the function [irtQ::est_item()].
#'
#' @export
#'
rdif.est_item <- function(x,
                          group,
                          focal.name,
                          item.skip = NULL,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("rdifrs", "rdifr", "rdifs"),
                          max.iter = 10,
                          min.resp = NULL,
                          method = "ML",
                          range = c(-5, 5),
                          norm.prior =
                            c(0, 1),
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
  dif_rst <-
    rdif_one(x = x, data = data, score = score, group = group,
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
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- rdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, item.skip = item.skip.puri, D = D,
          alpha = alpha
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
rdif_one <- function(x,
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
  # poolsd <- moments$poolsd

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
  # rdif_stats <- list(x = rdifr, y = rdifs)
  n_total <- n_foc + n_ref

  # Cohen's d
  # efs_cohen <-
  #   purrr::map2(.x=rdif_stats, .y=poolsd$cohen, .f=function(x, y) round(x / y, 4)) %>%
  #   data.frame() %>%
  #   dplyr::rename_all(.f=function(x) c("co_rdifr", "co_rdifs"))

  # Hedge's g
  # efs_hedge <-
  #   purrr::map2(
  #     .x = rdif_stats, .y = poolsd$hedge,
  #     .f = function(x, y) {
  #       round((x / y) * (1 - (3 / (4 * (n_foc + n_ref) - 9))), 4)
  #     }
  #   ) %>%
  #   data.frame() %>%
  #   dplyr::rename_all(.f = function(x) c("he_rdifr", "he_rdifs"))

  # Glass's delta
  # efs_glass <-
  #   purrr::map2(.x = rdif_stats, .y = poolsd$glass, .f = function(x, y) round(x / y, 4)) %>%
  #   data.frame() %>%
  #   dplyr::rename_all(.f = function(x) c("gl_rdifr", "gl_rdifs"))

  # combine all effect size
  # effect_size <- data.frame(id = x$id, efs_hedge, efs_glass, stringsAsFactors = FALSE)

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      rdifr = round(rdifr, 4), z.rdifr = round(z_stat_rdifr, 4),
      rdifs = round(rdifs, 4), z.rdifs = round(z_stat_rdifs, 4),
      rdifrs = round(chisq, 4),
      p.rdifr = p_rdifr, p.rdifs = p_rdifs, p.rdifrs = p_rdifrs,
      n.ref = n_ref, n.foc = n_foc, n.total = n_total, stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # when there are items that should be skipped for the DIF analysis
  # insert NAs to the corresponding results of the items
  if (!is.null(item.skip)) {
    stat_df[item.skip, 2:9] <- NA
    p_rdifr[item.skip] <- NA
    p_rdifs[item.skip] <- NA
    p_rdifrs[item.skip] <- NA
    moments_rdifr[item.skip, ] <- NA
    moments_rdifs[item.skip, ] <- NA
    covar[item.skip] <- NA
  }

  # find the flagged items
  dif_item_rdifr <- as.numeric(which(p_rdifr <= alpha))
  dif_item_rdifs <- as.numeric(which(p_rdifs <= alpha))
  dif_item_rdifrs <- which(p_rdifrs <= alpha)
  if(length(dif_item_rdifr) == 0) dif_item_rdifr <- NULL
  if(length(dif_item_rdifs) == 0) dif_item_rdifs <- NULL
  if(length(dif_item_rdifrs) == 0) dif_item_rdifrs <- NULL

  # summarize the results
  rst <- list(
    dif_stat = stat_df,
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
