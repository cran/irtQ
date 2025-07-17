#' Generalized IRT residual-based DIF detection framework for multiple groups
#' (GRDIF)
#'
#' This function computes three GRDIF statistics, \eqn{GRDIF_{R}},
#' \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}, for analyzing differential item
#' functioning (DIF) among multiple groups (Lim et al., 2024). They are
#' specialized to capture uniform DIF, nonuniform DIF, and mixed DIF,
#' respectively.
#'
#' @inheritParams rdif
#' @param score A numeric vector containing examinees' ability estimates (theta
#'   values). If not provided, [irtQ::grdif()] will estimate ability parameters
#'   internally before computing the GRDIF statistics. See [irtQ::est_score()]
#'   for more information on scoring methods. Default is `NULL`.
#' @param focal.name A numeric or character vector specifying the levels
#'   associated with the focal groups. For example, consider `group = c(0, 0, 1,
#'   2, 2, 3, 3)`, where '1', '2', and '3' indicate three distinct focal groups
#'   and '0' represents the reference group. In this case, set `focal.name =
#'   c(1, 2, 3)`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for hypothesis testing using the GRDIF statistics. Default is `0.05`.
#' @param purify.by A character string specifying which GRDIF statistic is used
#'   to perform the purification. Available options are "grdifrs" for
#'   \eqn{GRDIF_{RS}}, "grdifr" for \eqn{GRDIF_{R}}, and "grdifs" for
#'   \eqn{GRDIF_{S}}.
#' @param post.hoc A logical value indicating whether to perform post-hoc RDIF
#'   analyses for all possible pairwise group comparisons on items flagged as
#'   statistically significant. Default is TRUE. See details below.
#' @param ... Additional arguments passed to the [irtQ::est_score()] function.
#'
#' @details
#' The GRDIF framework (Lim et al., 2024) is a generalized version of the RDIF
#' detection framework, designed to assess DIF across multiple groups. The
#' framework includes three statistics: \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and
#' \eqn{GRDIF_{RS}}, which are tailored to detect uniform, nonuniform, and mixed
#' DIF, respectively. Under the null hypothesis that the test contains no DIF
#' items, the statistics \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}
#' asymptotically follow \eqn{\chi^{2}} distributions with G-1, G-1, and 2(G-1)
#' degrees of freedom, respectively, where *G* represents the number of
#' groups being compared. For more information on the GRDIF framework, see
#' Lim et al. (2024).
#'
#' The [irtQ::grdif()] function computes all three GRDIF statistics:
#' \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}. It supports both
#' dichotomous and polytomous item response data. To compute these statistics,
#' [irtQ::grdif()] requires: (1) item parameter estimates obtained from aggregate data
#' (regardless of group membership); (2) examinees' ability estimates (e.g., ML);
#' and (3) their item response data. Note that ability estimates must be computed
#' using the item parameters estimated from the aggregate data. The item parameters
#' should be provided via the `x` argument, the ability estimates via the `score`
#' argument, and the response data via the `data` argument. If ability estimates
#' are not supplied (`score = NULL`), [irtQ::grdif()] automatically computes them
#' using the scoring method specified in the `method` argument (e.g., `method = "ML"`).
#'
#' The `group` argument accepts a vector of numeric or character values,
#' indicating the group membership of examinees. The vector may contain multiple
#' distinct values, where one represents the reference group and the others
#' represent focal groups. Its length must match the number of rows in the
#' response data, with each value corresponding to an examinee’s group
#' membership. Once `group` is specified, a numeric or character vector must be
#' supplied via the `focal.name` argument to define which group(s) in `group`
#' represent the focal groups. The reference group is defined as the group not
#' included in `focal.name`.
#'
#' Similar to the original RDIF framework for two-group comparisons, the GRDIF
#' framework supports an iterative purification process. When `purify = TRUE`,
#' purification is conducted based on the GRDIF statistic specified in the
#' `purify.by` argument (e.g., `purify.by = "grdifrs"`). During each iteration,
#' examinees' latent abilities are re-estimated using only the purified items,
#' with the scoring method determined by the `method` argument. The process
#' continues until no additional DIF items are flagged or until the number of
#' iterations reaches the specified `max.iter` limit. For details on the
#' purification procedure, see Lim et al. (2022).
#'
#' Scoring based on a limited number of items can lead to large standard errors,
#' which may compromise the effectiveness of DIF detection within the GRDIF
#' framework. The `min.resp` argument can be used to exclude ability estimates
#' with substantial standard errors, especially during the purification process.
#' For example, if `min.resp` is not NULL (e.g., `min.resp = 5`), examinees whose
#' total number of responses falls below the specified threshold will have their
#' responses treated as missing values (i.e., NA). Consequently, their ability
#' estimates will also be missing and will not be used when computing the GRDIF
#' statistics. If `min.resp = NULL`, an examinee's score will be computed as long
#' as at least one response is available.
#'
#' The `post.hoc` argument enables post-hoc RDIF analyses across all possible
#' pairwise group comparisons for items flagged as statistically significant.
#' For instance, consider four groups of examinees: A, B, C, and D. If
#' `post.hoc = TRUE`, the [irtQ::grdif()] function will perform pairwise RDIF analyses
#' for each flagged item across all group pairs (A-B, A-C, A-D, B-C, B-D, and C-D).
#' This provides a more detailed understanding of which specific group pairs
#' exhibit DIF. Note that when purification is enabled (i.e., `purify = TRUE`),
#' post-hoc RDIF analyses are conducted for each flagged item at each iteration
#' of the purification process.
#'
#' @return This function returns a list containing four main components:
#'
#' \item{no_purify}{A list of sub-objects presenting the results of DIF analysis
#' without a purification procedure. These include:
#'   \describe{
#'     \item{dif_stat}{A data frame summarizing the results of the three GRDIF
#'     statistics for all evaluated items. Columns include item ID,
#'     \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}} statistics; their
#'     corresponding p-values; the sample size of the reference group;
#'     the sample sizes of the focal groups; and the total sample size.}
#'     \item{moments}{A list of three data frames reporting the moments of
#'     mean raw residuals (MRRs) and mean squared residuals (MSRs) across all
#'     compared groups. The first contains means, the second variances, and
#'     the third covariances of MRRs and MSRs.}
#'     \item{dif_item}{A list of three numeric vectors indicating the items
#'     flagged as potential DIF items by each of the GRDIF statistics
#'     (\eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}).}
#'     \item{score}{A numeric vector of ability estimates used to compute
#'     the GRDIF statistics.}
#'     \item{post.hoc}{A list of three data frames containing post-hoc RDIF
#'     analysis results for all possible pairwise group comparisons. Each
#'     data frame corresponds to the results for items flagged by
#'     \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}}, respectively.}
#'   }
#' }
#'
#' \item{purify}{A logical value indicating whether the purification process
#' was applied.}
#'
#' \item{with_purify}{A list of sub-objects presenting the results of DIF
#' analysis with a purification procedure. These include:
#'   \describe{
#'     \item{purify.by}{A character string indicating which GRDIF statistic
#'      was used for purification: "grdifr", "grdifs", or "grdifrs",
#'      corresponding to \eqn{GRDIF_{R}}, \eqn{GRDIF_{S}}, and \eqn{GRDIF_{RS}},
#'      respectively.}
#'     \item{dif_stat}{A data frame summarizing the GRDIF results across
#'     iterations. Columns include item ID, the three GRDIF statistics and
#'     their p-values, sample size of the reference group, sample sizes of
#'     the focal groups, total sample size, and the iteration number at which
#'     each statistic was computed.}
#'     \item{moments}{A list of three data frames showing the MRR and MSR
#'     moments across iterations. The final column in each data frame indicates
#'     the iteration in which the statistics were computed.}
#'     \item{n.iter}{The total number of iterations executed during the
#'     purification process.}
#'     \item{score}{A numeric vector of the final purified ability estimates
#'     used for computing GRDIF statistics.}
#'     \item{post.hoc}{A data frame containing the post-hoc RDIF analysis
#'     results for flagged items across all possible pairwise group comparisons,
#'     updated at each iteration.}
#'     \item{complete}{A logical value indicating whether the purification
#'     process was completed. If `FALSE`, the process reached the maximum number
#'     of iterations without convergence.}
#'   }
#' }
#'
#' \item{alpha}{The significance level (\eqn{\alpha}) used for hypothesis
#' testing of the GRDIF statistics.}
#'
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::rdif()] [irtQ::est_irt()], [irtQ::est_item()],
#' [irtQ::simdat()], [irtQ::shape_df()], [irtQ::est_score()]
#'
#' @references Lim, H., & Choe, E. M. (2023). Detecting differential item
#'   functioning in CAT using IRT residual DIF approach.
#'  *Journal of Educational Measurement, 60*(4), 626-650. \doi{doi:10.1111/jedm.12366}.
#'
#'   Lim, H., Choe, E. M., & Han, K. T. (2022). A residual-based differential
#'   item functioning detection framework in item response theory. *Journal of
#'   Educational Measurement, 59*(1), 80-104. \doi{doi:10.1111/jedm.12313}.
#'
#'   Lim, H., Zhu, D., Choe, E. M., & Han, K. T. (2024). Detecting
#'   differential item functioning among multiple groups using IRT residual DIF
#'   framework. *Journal of Educational Measurement, 61*(4), 656-681.
#'
#' @examples
#' \donttest{
#' # Load required library
#' library("dplyr")
#'
#' ## Uniform DIF detection for four groups (1 reference, 3 focal)
#' ########################################################
#' # (1) Manipulate uniform DIF for all three focal groups
#' ########################################################
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select 36 non-DIF items modeled under 3PLM
#' par_nstd <-
#'   bring.flexmirt(file = flex_sam, "par")$Group1$full_df %>%
#'   dplyr::filter(.data$model == "3PLM") %>%
#'   dplyr::filter(dplyr::row_number() %in% 1:36) %>%
#'   dplyr::select(1:6)
#' par_nstd$id <- paste0("nondif", 1:36)
#'
#' # Generate four new items on which uniform DIF will be imposed
#' difpar_ref <-
#'   shape_df(
#'     par.drm = list(a = c(0.8, 1.5, 0.8, 1.5), b = c(0.0, 0.0, -0.5, -0.5), g = .15),
#'     item.id = paste0("dif", 1:4), cats = 2, model = "3PLM"
#'   )
#'
#' # Introduce DIF by shifting b-parameters differently for each focal group
#' difpar_foc1 <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(0.7, 0.7, 0, 0))
#' difpar_foc2 <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(0, 0, 0.7, 0.7))
#' difpar_foc3 <-
#'   difpar_ref %>%
#'   dplyr::mutate_at(.vars = "par.2", .funs = function(x) x + c(-0.4, -0.4, -0.5, -0.5))
#'
#' # Combine the 4 DIF and 36 non-DIF items for all four groups
#' # Therefore, the first four items contain uniform DIF across all focal groups
#' par_ref <- rbind(difpar_ref, par_nstd)
#' par_foc1 <- rbind(difpar_foc1, par_nstd)
#' par_foc2 <- rbind(difpar_foc2, par_nstd)
#' par_foc3 <- rbind(difpar_foc3, par_nstd)
#'
#' # Generate true abilities from different distributions
#' set.seed(128)
#' theta_ref <- rnorm(500, 0.0, 1.0)
#' theta_foc1 <- rnorm(500, -1.0, 1.0)
#' theta_foc2 <- rnorm(500, 1.0, 1.0)
#' theta_foc3 <- rnorm(500, 0.5, 1.0)
#'
#' # Simulate response data for each group
#' resp_ref <- irtQ::simdat(par_ref, theta = theta_ref, D = 1)
#' resp_foc1 <- irtQ::simdat(par_foc1, theta = theta_foc1, D = 1)
#' resp_foc2 <- irtQ::simdat(par_foc2, theta = theta_foc2, D = 1)
#' resp_foc3 <- irtQ::simdat(par_foc3, theta = theta_foc3, D = 1)
#' data <- rbind(resp_ref, resp_foc1, resp_foc2, resp_foc3)
#'
#' ########################################################
#' # (2) Estimate item and ability parameters
#' #     using aggregated data
#' ########################################################
#'
#' # Estimate item parameters
#' est_mod <- irtQ::est_irt(data = data, D = 1, model = "3PLM")
#' est_par <- est_mod$par.est
#'
#' # Estimate ability parameters using MLE
#' score <- irtQ::est_score(x = est_par, data = data, method = "ML")$est.theta
#'
#' ########################################################
#' # (3) Conduct DIF analysis
#' ########################################################
#'
#' # Create a group membership vector:
#' # 0 = reference group; 1, 2, 3 = focal groups
#' group <- c(rep(0, 500), rep(1, 500), rep(2, 500), rep(3, 500))
#'
#' # (a) Compute GRDIF statistics without purification,
#' #     and perform post-hoc pairwise comparisons for flagged items
#' dif_nopuri <- grdif(
#'   x = est_par, data = data, score = score, group = group,
#'   focal.name = c(1, 2, 3), D = 1, alpha = 0.05,
#'   purify = FALSE, post.hoc = TRUE
#' )
#' print(dif_nopuri)
#'
#' # Display post-hoc pairwise comparison results
#' print(dif_nopuri$no_purify$post.hoc)
#'
#' # (b) Compute GRDIF statistics with purification
#' #     based on GRDIF_R, including post-hoc comparisons
#' dif_puri_r <- grdif(
#'   x = est_par, data = data, score = score, group = group,
#'   focal.name = c(1, 2, 3), D = 1, alpha = 0.05,
#'   purify = TRUE, purify.by = "grdifr", post.hoc = TRUE
#' )
#' print(dif_puri_r)
#'
#' # Display post-hoc results before purification
#' print(dif_puri_r$no_purify$post.hoc)
#'
#' # Display post-hoc results after purification
#' print(dif_puri_r$with_purify$post.hoc)
#' }
#'
#' @export
grdif <- function(x, ...) UseMethod("grdif")

#' @describeIn grdif Default method to compute the three GRDIF statistics for
#' multiple-group data using a data frame `x` that contains item metadata.
#'
#' @import dplyr
#' @export
grdif.default <- function(x,
                          data,
                          score = NULL,
                          group,
                          focal.name,
                          D = 1,
                          alpha = 0.05,
                          missing = NA,
                          purify = FALSE,
                          purify.by = c("grdifrs", "grdifr", "grdifs"),
                          max.iter = 10,
                          min.resp = NULL,
                          post.hoc = TRUE,
                          method = "ML",
                          range = c(-4, 4),
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
      # count the number of responses for each question
      n_resp <- Rfast::rowsums(!is.na(data))

      # find the questions that have less than the minimum number of responses
      loc_less <- which(n_resp < min.resp & n_resp > 0)

      # replace the responses for those questions with NA
      data[loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <- grdif_one(x = x, data = data, score = score, group = group, focal.name = focal.name, D = D, alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(
    dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL,
    post.hoc = list(by.grdifr = NULL, by.grdifs = NULL, by.grdifrs = NULL)
  )
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, post.hoc = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <-
    list(
      mu = data.frame(id = x$id, dif_rst$moments$mu, stringsAsFactors = FALSE),
      var = data.frame(id = x$id, dif_rst$moments$var, stringsAsFactors = FALSE),
      covariance = data.frame(id = x$id, dif_rst$moments$covariance, stringsAsFactors = FALSE)
    )
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    n_df <-
      no_purify$dif_stat %>%
      dplyr::select(dplyr::contains("n.")) %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        grdifr = NA,
        grdifs = NA,
        grdifrs = NA,
        p.grdifr = NA,
        p.grdifs = NA,
        p.grdifrs = NA,
        n_df,
        n.iter = NA,
        stringsAsFactors = FALSE
      )
    mu_empty <-
      dif_rst$moments$mu %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    var_empty <-
      dif_rst$moments$var %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    covar_empty <-
      dif_rst$moments$covariance %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    mmt_df <-
      list(
        mu = data.frame(id = rep(NA_character_, nrow(x)), mu_empty, n.iter = NA, stringsAsFactors = FALSE),
        var = data.frame(id = x$id, var_empty, n.iter = NA, stringsAsFactors = FALSE),
        covariance = data.frame(id = x$id, covar_empty, n.iter = NA, stringsAsFactors = FALSE)
      )
    ncol.stat <- ncol(dif_stat)
    ncol.mu <- ncol.var <- ncol(mmt_df$mu)
    ncol.covar <- ncol(mmt_df$covariance)

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data
    score_puri <- score

    # create an empty list to contain the post-hoc list
    if (post.hoc) {
      post_df_pury <- list()
    }

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

        # a flagged item which has the largest significant GDIF statistic
        flag_max <-
          switch(purify.by,
                 grdifr = which.max(dif_stat_tmp$grdifr),
                 grdifs = which.max(dif_stat_tmp$grdifs),
                 grdifrs = which.max(dif_stat_tmp$grdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:(ncol.stat - 1)] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, ncol.stat] <- i - 1
        mmt_df$mu[del_item, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu[flag_max, ]
        mmt_df$var[del_item, 1:(ncol.var - 1)] <- mmt_df_tmp$var[flag_max, ]
        mmt_df$covariance[del_item, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance[flag_max, ]
        mmt_df$mu[del_item, ncol.mu] <- i - 1
        mmt_df$var[del_item, ncol.var] <- i - 1
        mmt_df$covariance[del_item, ncol.covar] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # post-hoc DIF analysis one flagged item by one flagged item
        if (post.hoc) {
          # select the flagged items for the post-hoc DIF analyses
          item4post <- flag_max

          # item metadata for the flagged items
          x_post <- x_puri[flag_max, ]

          # find all pairs of the groups
          pair.g <- utils::combn(x = sort(unique(group)), m = 2, simplify = FALSE)

          # post-hoc DIF analyses
          post.stats <- list()
          for (p in 1:length(pair.g)) {
            # select a pair of two groups
            pair.tmp <- pair.g[[p]]
            group.pair <- paste(pair.tmp, collapse = " & ")

            # find the location of examinees for the two compared groups
            loc_g1 <- which(group == pair.tmp[1])
            loc_g2 <- which(group == pair.tmp[2])
            loc.tmp <- c(loc_g1, loc_g2)
            focal.name.tmp <- pair.tmp[2]

            # select the group variables for the two groups
            group.tmp <- group[loc.tmp]

            # select the response data for the two groups
            data.tmp <- data_puri[c(loc_g1, loc_g2), item4post, drop = FALSE]

            # select the thetas for the two groups
            score.tmp <- score_puri[c(loc_g1, loc_g2)]

            # DIF analysis
            dif_post.tmp <-
              rdif_one(
                x = x_post, data = data.tmp, score = score.tmp, group = group.tmp,
                focal.name = focal.name.tmp, D = D, alpha = alpha
              )

            # add the stats to the list
            post.stats[[p]] <-
              dif_post.tmp$dif_stat %>%
              dplyr::mutate(group.pair = group.pair) %>%
              dplyr::relocate("group.pair", .after = "id")
          }

          # re-organize the stat list
          post_df_pury[[i]] <-
            purrr::map_df(.x = post.stats, ~ {
              .x
            }) %>%
            dplyr::mutate(n_iter = i - 1)
        }

        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- Rfast::rowsums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- grdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, D = D, alpha = alpha
        )

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <-
          list(
            mu = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$mu, stringsAsFactors = FALSE),
            var = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$var, stringsAsFactors = FALSE),
            covariance = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$covariance, stringsAsFactors = FALSE)
          )

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
          dif_stat[item_num, ncol.stat] <- i
          mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
          mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
          mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
          mmt_df$mu[item_num, ncol.mu] <- i
          mmt_df$var[item_num, ncol.var] <- i
          mmt_df$covariance[item_num, ncol.covar] <- i

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
        dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
        dif_stat[item_num, ncol.stat] <- i
        mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
        mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
        mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
        mmt_df$mu[item_num, ncol.mu] <- i
        mmt_df$var[item_num, ncol.var] <- i
        mmt_df$covariance[item_num, ncol.covar] <- i
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
      if (post.hoc) {
        with_purify$post.hoc <-
          post_df_pury %>%
          dplyr::bind_rows()
      }
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <-
        purrr::map(.x = no_purify$moments, ~ {
          cbind(.x, n.iter = 0)
        })
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }


  # post-hoc DIF analysis for the non-purified DIF analysis results
  if (post.hoc) {
    # select the flagged items for the post-hoc DIF analyses
    item4post <- no_purify$dif_item

    for (i in 1:3) {
      if (length(item4post[[i]]) != 0) {
        # item metadata for the flagged items
        x_post <- x[item4post[[i]], ]

        # find all pairs of the groups
        pair.g <- utils::combn(x = sort(unique(group)), m = 2, simplify = FALSE)

        # post-hoc DIF analyses
        post.stats <- list()
        for (p in 1:length(pair.g)) {
          # p <- 1
          # select a pair of two groups
          pair.tmp <- pair.g[[p]]
          group.pair <- paste(pair.tmp, collapse = " & ")

          # find the location of examinees for the two compared groups
          loc_g1 <- which(group == pair.tmp[1])
          loc_g2 <- which(group == pair.tmp[2])
          loc.tmp <- c(loc_g1, loc_g2)
          focal.name.tmp <- pair.tmp[2]

          # select the group variables for the two groups
          group.tmp <- group[loc.tmp]

          # select the response data for the two groups
          data.tmp <- data[c(loc_g1, loc_g2), item4post[[i]], drop = FALSE]

          # select the thetas for the two groups
          score.tmp <- score[c(loc_g1, loc_g2)]

          # DIF analysis
          dif_post.tmp <-
            rdif_one(
              x = x_post, data = data.tmp, score = score.tmp, group = group.tmp,
              focal.name = focal.name.tmp, D = D, alpha = alpha
            )

          # add the stats to the list
          post.stats[[p]] <-
            dif_post.tmp$dif_stat %>%
            dplyr::mutate(group.pair = group.pair) %>%
            dplyr::relocate("group.pair", .after = "id")
        }

        # re-organize the stat list
        post_df_nopury <-
          purrr::map(
            .x = 1:length(item4post[[i]]),
            .f = function(i) {
              purrr::map_df(.x = post.stats, ~ {
                .x[i, ]
              })
            }
          ) %>%
          dplyr::bind_rows()
        no_purify$post.hoc[[i]] <- post_df_nopury
      }
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "grdif"
  rst$call <- cl
  rst
}

#' @describeIn grdif An object created by the function [irtQ::est_irt()].
#' @import dplyr
#' @export
#'
grdif.est_irt <- function(x, score = NULL, group, focal.name, alpha = 0.05, missing = NA, purify = FALSE,
                          purify.by = c("grdifrs", "grdifr", "grdifs"), max.iter = 10, min.resp = NULL, post.hoc = TRUE,
                          method = "ML", range = c(-4, 4), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1,
                          verbose = TRUE, ...) {
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
      # count the number of responses for each question
      n_resp <- Rfast::rowsums(!is.na(data))

      # find the questions that have less than the minimum number of responses
      loc_less <- which(n_resp < min.resp & n_resp > 0)

      # replace the responses for those questions with NA
      data[loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <- grdif_one(x = x, data = data, score = score, group = group, focal.name = focal.name, D = D, alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(
    dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL,
    post.hoc = list(by.grdifr = NULL, by.grdifs = NULL, by.grdifrs = NULL)
  )
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, post.hoc = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <-
    list(
      mu = data.frame(id = x$id, dif_rst$moments$mu, stringsAsFactors = FALSE),
      var = data.frame(id = x$id, dif_rst$moments$var, stringsAsFactors = FALSE),
      covariance = data.frame(id = x$id, dif_rst$moments$covariance, stringsAsFactors = FALSE)
    )
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    n_df <-
      no_purify$dif_stat %>%
      dplyr::select(dplyr::contains("n.")) %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        grdifr = NA,
        grdifs = NA,
        grdifrs = NA,
        p.grdifr = NA,
        p.grdifs = NA,
        p.grdifrs = NA,
        n_df,
        n.iter = NA,
        stringsAsFactors = FALSE
      )
    mu_empty <-
      dif_rst$moments$mu %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    var_empty <-
      dif_rst$moments$var %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    covar_empty <-
      dif_rst$moments$covariance %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    mmt_df <-
      list(
        mu = data.frame(id = rep(NA_character_, nrow(x)), mu_empty, n.iter = NA, stringsAsFactors = FALSE),
        var = data.frame(id = x$id, var_empty, n.iter = NA, stringsAsFactors = FALSE),
        covariance = data.frame(id = x$id, covar_empty, n.iter = NA, stringsAsFactors = FALSE)
      )
    ncol.stat <- ncol(dif_stat)
    ncol.mu <- ncol.var <- ncol(mmt_df$mu)
    ncol.covar <- ncol(mmt_df$covariance)

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data
    score_puri <- score

    # create an empty list to contain the post-hoc list
    if (post.hoc) {
      post_df_pury <- list()
    }

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

        # a flagged item which has the largest significant GDIF statistic
        flag_max <-
          switch(purify.by,
                 grdifr = which.max(dif_stat_tmp$grdifr),
                 grdifs = which.max(dif_stat_tmp$grdifs),
                 grdifrs = which.max(dif_stat_tmp$grdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:(ncol.stat - 1)] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, ncol.stat] <- i - 1
        mmt_df$mu[del_item, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu[flag_max, ]
        mmt_df$var[del_item, 1:(ncol.var - 1)] <- mmt_df_tmp$var[flag_max, ]
        mmt_df$covariance[del_item, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance[flag_max, ]
        mmt_df$mu[del_item, ncol.mu] <- i - 1
        mmt_df$var[del_item, ncol.var] <- i - 1
        mmt_df$covariance[del_item, ncol.covar] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # post-hoc DIF analysis one flagged item by one flagged item
        if (post.hoc) {
          # select the flagged items for the post-hoc DIF analyses
          item4post <- flag_max

          # item metadata for the flagged items
          x_post <- x_puri[flag_max, ]

          # find all pairs of the groups
          pair.g <- utils::combn(x = sort(unique(group)), m = 2, simplify = FALSE)

          # post-hoc DIF analyses
          post.stats <- list()
          for (p in 1:length(pair.g)) {
            # select a pair of two groups
            pair.tmp <- pair.g[[p]]
            group.pair <- paste(pair.tmp, collapse = " & ")

            # find the location of examinees for the two compared groups
            loc_g1 <- which(group == pair.tmp[1])
            loc_g2 <- which(group == pair.tmp[2])
            loc.tmp <- c(loc_g1, loc_g2)
            focal.name.tmp <- pair.tmp[2]

            # select the group variables for the two groups
            group.tmp <- group[loc.tmp]

            # select the response data for the two groups
            data.tmp <- data_puri[c(loc_g1, loc_g2), item4post, drop = FALSE]

            # select the thetas for the two groups
            score.tmp <- score_puri[c(loc_g1, loc_g2)]

            # DIF analysis
            dif_post.tmp <-
              rdif_one(
                x = x_post, data = data.tmp, score = score.tmp, group = group.tmp,
                focal.name = focal.name.tmp, D = D, alpha = alpha
              )

            # add the stats to the list
            post.stats[[p]] <-
              dif_post.tmp$dif_stat %>%
              dplyr::mutate(group.pair = group.pair) %>%
              dplyr::relocate("group.pair", .after = "id")
          }

          # re-organize the stat list
          post_df_pury[[i]] <-
            purrr::map_df(.x = post.stats, ~ {
              .x
            }) %>%
            dplyr::mutate(n_iter = i - 1)
        }

        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- Rfast::rowsums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- grdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, D = D, alpha = alpha
        )

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <-
          list(
            mu = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$mu, stringsAsFactors = FALSE),
            var = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$var, stringsAsFactors = FALSE),
            covariance = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$covariance, stringsAsFactors = FALSE)
          )

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
          dif_stat[item_num, ncol.stat] <- i
          mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
          mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
          mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
          mmt_df$mu[item_num, ncol.mu] <- i
          mmt_df$var[item_num, ncol.var] <- i
          mmt_df$covariance[item_num, ncol.covar] <- i

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
        dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
        dif_stat[item_num, ncol.stat] <- i
        mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
        mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
        mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
        mmt_df$mu[item_num, ncol.mu] <- i
        mmt_df$var[item_num, ncol.var] <- i
        mmt_df$covariance[item_num, ncol.covar] <- i
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
      if (post.hoc) {
        with_purify$post.hoc <-
          post_df_pury %>%
          dplyr::bind_rows()
      }
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <-
        purrr::map(.x = no_purify$moments, ~ {
          cbind(.x, n.iter = 0)
        })
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }


  # post-hoc DIF analysis for the non-purified DIF analysis results
  if (post.hoc) {
    # select the flagged items for the post-hoc DIF analyses
    item4post <- no_purify$dif_item

    for (i in 1:3) {
      if (length(item4post[[i]]) != 0) {
        # item metadata for the flagged items
        x_post <- x[item4post[[i]], ]

        # find all pairs of the groups
        pair.g <- utils::combn(x = sort(unique(group)), m = 2, simplify = FALSE)

        # post-hoc DIF analyses
        post.stats <- list()
        for (p in 1:length(pair.g)) {
          # p <- 1
          # select a pair of two groups
          pair.tmp <- pair.g[[p]]
          group.pair <- paste(pair.tmp, collapse = " & ")

          # find the location of examinees for the two compared groups
          loc_g1 <- which(group == pair.tmp[1])
          loc_g2 <- which(group == pair.tmp[2])
          loc.tmp <- c(loc_g1, loc_g2)
          focal.name.tmp <- pair.tmp[2]

          # select the group variables for the two groups
          group.tmp <- group[loc.tmp]

          # select the response data for the two groups
          data.tmp <- data[c(loc_g1, loc_g2), item4post[[i]], drop = FALSE]

          # select the thetas for the two groups
          score.tmp <- score[c(loc_g1, loc_g2)]

          # DIF analysis
          dif_post.tmp <-
            rdif_one(
              x = x_post, data = data.tmp, score = score.tmp, group = group.tmp,
              focal.name = focal.name.tmp, D = D, alpha = alpha
            )

          # add the stats to the list
          post.stats[[p]] <-
            dif_post.tmp$dif_stat %>%
            dplyr::mutate(group.pair = group.pair) %>%
            dplyr::relocate("group.pair", .after = "id")
        }

        # re-organize the stat list
        post_df_nopury <-
          purrr::map(
            .x = 1:length(item4post[[i]]),
            .f = function(i) {
              purrr::map_df(.x = post.stats, ~ {
                .x[i, ]
              })
            }
          ) %>%
          dplyr::bind_rows()
        no_purify$post.hoc[[i]] <- post_df_nopury
      }
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "grdif"
  rst$call <- cl
  rst
}


#' @describeIn grdif An object created by the function [irtQ::est_item()].
#' @import dplyr
#' @export
#'
grdif.est_item <- function(x, group, focal.name, alpha = 0.05, missing = NA, purify = FALSE,
                           purify.by = c("grdifrs", "grdifr", "grdifs"), max.iter = 10, min.resp = NULL, post.hoc = TRUE,
                           method = "ML", range = c(-4, 4), norm.prior = c(0, 1), nquad = 41, weights = NULL, ncore = 1,
                           verbose = TRUE, ...) {
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
      # count the number of responses for each question
      n_resp <- Rfast::rowsums(!is.na(data))

      # find the questions that have less than the minimum number of responses
      loc_less <- which(n_resp < min.resp & n_resp > 0)

      # replace the responses for those questions with NA
      data[loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...
    )$est.theta
  }

  # a) when no purification is set
  # do only one iteration of DIF analysis
  dif_rst <- grdif_one(x = x, data = data, score = score, group = group, focal.name = focal.name, D = D, alpha = alpha)

  # create two empty lists to contain the results
  no_purify <- list(
    dif_stat = NULL, moments = NULL, dif_item = NULL, score = NULL,
    post.hoc = list(by.grdifr = NULL, by.grdifs = NULL, by.grdifrs = NULL)
  )
  with_purify <- list(
    purify.by = NULL, dif_stat = NULL, moments = NULL,
    dif_item = NULL, n.iter = NULL, score = NULL, post.hoc = NULL, complete = NULL
  )

  # record the first DIF detection results into the no purification list
  no_purify$dif_stat <- dif_rst$dif_stat
  no_purify$dif_item <- dif_rst$dif_item
  no_purify$moments <-
    list(
      mu = data.frame(id = x$id, dif_rst$moments$mu, stringsAsFactors = FALSE),
      var = data.frame(id = x$id, dif_rst$moments$var, stringsAsFactors = FALSE),
      covariance = data.frame(id = x$id, dif_rst$moments$covariance, stringsAsFactors = FALSE)
    )
  no_purify$score <- score

  # when purification is used
  if (purify) {
    # verify the criterion for purification
    purify.by <- match.arg(purify.by)

    # create an empty vector and empty data frames
    # to contain the detected DIF items, statistics, and moments
    dif_item <- NULL
    n_df <-
      no_purify$dif_stat %>%
      dplyr::select(dplyr::contains("n.")) %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    dif_stat <-
      data.frame(
        id = rep(NA_character_, nrow(x)),
        grdifr = NA,
        grdifs = NA,
        grdifrs = NA,
        p.grdifr = NA,
        p.grdifs = NA,
        p.grdifrs = NA,
        n_df,
        n.iter = NA,
        stringsAsFactors = FALSE
      )
    mu_empty <-
      dif_rst$moments$mu %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    var_empty <-
      dif_rst$moments$var %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    covar_empty <-
      dif_rst$moments$covariance %>%
      dplyr::mutate_all(~ {
        NA_integer_
      })
    mmt_df <-
      list(
        mu = data.frame(id = rep(NA_character_, nrow(x)), mu_empty, n.iter = NA, stringsAsFactors = FALSE),
        var = data.frame(id = x$id, var_empty, n.iter = NA, stringsAsFactors = FALSE),
        covariance = data.frame(id = x$id, covar_empty, n.iter = NA, stringsAsFactors = FALSE)
      )
    ncol.stat <- ncol(dif_stat)
    ncol.mu <- ncol.var <- ncol(mmt_df$mu)
    ncol.covar <- ncol(mmt_df$covariance)

    # extract the first DIF analysis results
    # and check if at least one DIF item is detected
    dif_item_tmp <- dif_rst$dif_item[[purify.by]]
    dif_stat_tmp <- dif_rst$dif_stat
    mmt_df_tmp <- no_purify$moments

    # copy the response data and item meta data
    x_puri <- x
    data_puri <- data
    score_puri <- score

    # create an empty list to contain the post-hoc list
    if (post.hoc) {
      post_df_pury <- list()
    }

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

        # a flagged item which has the largest significant GDIF statistic
        flag_max <-
          switch(purify.by,
            grdifr = which.max(dif_stat_tmp$grdifr),
            grdifs = which.max(dif_stat_tmp$grdifs),
            grdifrs = which.max(dif_stat_tmp$grdifrs)
          )

        # check an item that is deleted
        del_item <- item_num[flag_max]

        # add the deleted item as the DIF item
        dif_item <- c(dif_item, del_item)

        # add the DIF statistics and moments for the detected DIF item
        dif_stat[del_item, 1:(ncol.stat - 1)] <- dif_stat_tmp[flag_max, ]
        dif_stat[del_item, ncol.stat] <- i - 1
        mmt_df$mu[del_item, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu[flag_max, ]
        mmt_df$var[del_item, 1:(ncol.var - 1)] <- mmt_df_tmp$var[flag_max, ]
        mmt_df$covariance[del_item, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance[flag_max, ]
        mmt_df$mu[del_item, ncol.mu] <- i - 1
        mmt_df$var[del_item, ncol.var] <- i - 1
        mmt_df$covariance[del_item, ncol.covar] <- i - 1

        # refine the leftover items
        item_num <- item_num[-flag_max]

        # post-hoc DIF analysis one flagged item by one flagged item
        if (post.hoc) {
          # select the flagged items for the post-hoc DIF analyses
          item4post <- flag_max

          # item metadata for the flagged items
          x_post <- x_puri[flag_max, ]

          # find all pairs of the groups
          pair.g <- utils::combn(x = sort(unique(group)), m = 2, simplify = FALSE)

          # post-hoc DIF analyses
          post.stats <- list()
          for (p in 1:length(pair.g)) {
            # select a pair of two groups
            pair.tmp <- pair.g[[p]]
            group.pair <- paste(pair.tmp, collapse = " & ")

            # find the location of examinees for the two compared groups
            loc_g1 <- which(group == pair.tmp[1])
            loc_g2 <- which(group == pair.tmp[2])
            loc.tmp <- c(loc_g1, loc_g2)
            focal.name.tmp <- pair.tmp[2]

            # select the group variables for the two groups
            group.tmp <- group[loc.tmp]

            # select the response data for the two groups
            data.tmp <- data_puri[c(loc_g1, loc_g2), item4post, drop = FALSE]

            # select the thetas for the two groups
            score.tmp <- score_puri[c(loc_g1, loc_g2)]

            # DIF analysis
            dif_post.tmp <-
              rdif_one(
                x = x_post, data = data.tmp, score = score.tmp, group = group.tmp,
                focal.name = focal.name.tmp, D = D, alpha = alpha
              )

            # add the stats to the list
            post.stats[[p]] <-
              dif_post.tmp$dif_stat %>%
              dplyr::mutate(group.pair = group.pair) %>%
              dplyr::relocate("group.pair", .after = "id")
          }

          # re-organize the stat list
          post_df_pury[[i]] <-
            purrr::map_df(.x = post.stats, ~ {
              .x
            }) %>%
            dplyr::mutate(n_iter = i - 1)
        }

        # remove the detected DIF item data which has the largest statistic from the item metadata
        x_puri <- x_puri[-flag_max, ]

        # remove the detected DIF item data which has the largest statistic from the response data
        data_puri <- data_puri[, -flag_max]

        # if min.resp is not NULL, find the examinees who have the number of responses
        # less than specified value (e.g., 5). Then, replace their all responses with NA
        if (!is.null(min.resp)) {
          n_resp <- Rfast::rowsums(!is.na(data_puri))
          loc_less <- which(n_resp < min.resp & n_resp > 0)
          data_puri[loc_less, ] <- NA
        }

        # compute the updated ability estimates after deleting the detected DIF item data
        score_puri <- est_score(
          x = x_puri, data = data_puri, D = D, method = method, range = range, norm.prior = norm.prior,
          nquad = nquad, weights = weights, ncore = ncore, ...
        )$est.theta

        # do DIF analysis using the updated ability estimates
        dif_rst_tmp <- grdif_one(
          x = x_puri, data = data_puri, score = score_puri, group = group,
          focal.name = focal.name, D = D, alpha = alpha
        )

        # extract the first DIF analysis results
        # and check if at least one DIF item is detected
        dif_item_tmp <- dif_rst_tmp$dif_item[[purify.by]]
        dif_stat_tmp <- dif_rst_tmp$dif_stat
        mmt_df_tmp <-
          list(
            mu = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$mu, stringsAsFactors = FALSE),
            var = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$var, stringsAsFactors = FALSE),
            covariance = data.frame(id = dif_rst_tmp$dif_stat$id, dif_rst_tmp$moments$covariance, stringsAsFactors = FALSE)
          )

        # check if a further DIF item is flagged
        if (is.null(dif_item_tmp)) {
          # add no additional DIF item
          dif_item <- dif_item

          # add the DIF statistics for rest of items
          dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
          dif_stat[item_num, ncol.stat] <- i
          mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
          mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
          mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
          mmt_df$mu[item_num, ncol.mu] <- i
          mmt_df$var[item_num, ncol.var] <- i
          mmt_df$covariance[item_num, ncol.covar] <- i

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
        dif_stat[item_num, 1:(ncol.stat - 1)] <- dif_stat_tmp
        dif_stat[item_num, ncol.stat] <- i
        mmt_df$mu[item_num, 1:(ncol.mu - 1)] <- mmt_df_tmp$mu
        mmt_df$var[item_num, 1:(ncol.var - 1)] <- mmt_df_tmp$var
        mmt_df$covariance[item_num, 1:(ncol.covar - 1)] <- mmt_df_tmp$covariance
        mmt_df$mu[item_num, ncol.mu] <- i
        mmt_df$var[item_num, ncol.var] <- i
        mmt_df$covariance[item_num, ncol.covar] <- i
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
      if (post.hoc) {
        with_purify$post.hoc <-
          post_df_pury %>%
          dplyr::bind_rows()
      }
    } else {
      # in case when no DIF item is detected from the first DIF analysis results
      with_purify$purify.by <- purify.by
      with_purify$dif_stat <- cbind(no_purify$dif_stat, n.iter = 0)
      with_purify$moments <-
        purrr::map(.x = no_purify$moments, ~ {
          cbind(.x, n.iter = 0)
        })
      with_purify$n.iter <- 0
      with_purify$complete <- TRUE
    }
  }


  # post-hoc DIF analysis for the non-purified DIF analysis results
  if (post.hoc) {
    # select the flagged items for the post-hoc DIF analyses
    item4post <- no_purify$dif_item

    for (i in 1:3) {
      if (length(item4post[[i]]) != 0) {
        # item metadata for the flagged items
        x_post <- x[item4post[[i]], ]

        # find all pairs of the groups
        pair.g <- utils::combn(x = sort(unique(group)), m = 2, simplify = FALSE)

        # post-hoc DIF analyses
        post.stats <- list()
        for (p in 1:length(pair.g)) {
          # p <- 1
          # select a pair of two groups
          pair.tmp <- pair.g[[p]]
          group.pair <- paste(pair.tmp, collapse = " & ")

          # find the location of examinees for the two compared groups
          loc_g1 <- which(group == pair.tmp[1])
          loc_g2 <- which(group == pair.tmp[2])
          loc.tmp <- c(loc_g1, loc_g2)
          focal.name.tmp <- pair.tmp[2]

          # select the group variables for the two groups
          group.tmp <- group[loc.tmp]

          # select the response data for the two groups
          data.tmp <- data[c(loc_g1, loc_g2), item4post[[i]], drop = FALSE]

          # select the thetas for the two groups
          score.tmp <- score[c(loc_g1, loc_g2)]

          # DIF analysis
          dif_post.tmp <-
            rdif_one(
              x = x_post, data = data.tmp, score = score.tmp, group = group.tmp,
              focal.name = focal.name.tmp, D = D, alpha = alpha
            )

          # add the stats to the list
          post.stats[[p]] <-
            dif_post.tmp$dif_stat %>%
            dplyr::mutate(group.pair = group.pair) %>%
            dplyr::relocate("group.pair", .after = "id")
        }

        # re-organize the stat list
        post_df_nopury <-
          purrr::map(
            .x = 1:length(item4post[[i]]),
            .f = function(i) {
              purrr::map_df(.x = post.stats, ~ {
                .x[i, ]
              })
            }
          ) %>%
          dplyr::bind_rows()
        no_purify$post.hoc[[i]] <- post_df_nopury
      }
    }
  }

  # summarize the results
  rst <- list(no_purify = no_purify, purify = purify, with_purify = with_purify, alpha = alpha)

  # return the DIF detection results
  class(rst) <- "grdif"
  rst$call <- cl
  rst
}


# This function conducts one iteration of GRDIF analysis using the IRT residual based statistics
grdif_one <- function(x, data, score, group, focal.name, D = 1, alpha = 0.05) {
  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # check the unique score categories
  cats <- elm_item$cats

  # check the number of items
  nitem <- length(cats)

  # count N of the focal groups
  n.fg <- length(unique(focal.name))

  ## ---------------------------------
  # compute the GRDIF statistics
  ## ---------------------------------
  # find the location of examinees for the reference and the focal groups
  loc_ref <- which(!group %in% focal.name)
  loc_foc <-
    purrr::map(.x = focal.name, ~ {
      which(group == .x)
    })

  # divide the response data into the multiple groups data
  resp_ref <- data[loc_ref, , drop = FALSE]
  resp_foc <- purrr::map(.x = loc_foc, ~ {
    data[.x, , drop = FALSE]
  })

  # count sample sizes
  n_ref <- Rfast::colsums(!is.na(resp_ref))
  n_foc <- purrr::map(.x = resp_foc, ~ {
    Rfast::colsums(!is.na(.x))
  })

  # check if an item has all missing data for either of multiple groups
  all_miss <- sort(unique(c(
    which(n_ref == 0),
    unlist(purrr::map(.x = n_foc, ~ {
      which(.x == 0)
    }))
  )))

  # divide the thetas into multiple groups data
  score_ref <- score[loc_ref]
  score_foc <- purrr::map(.x = loc_foc, ~ {
    score[.x]
  })

  # compute the model-predicted probability of answering correctly (a.k.a. model-expected item score)
  extscore_ref <- trace(elm_item = elm_item, theta = score_ref, D = D, tcc = TRUE)$icc
  extscore_foc <-
    purrr::map(.x = score_foc, ~ {
      trace(elm_item = elm_item, theta = .x, D = D, tcc = TRUE)$icc
    })

  # compute the model probability of score categories
  prob_ref <- trace(elm_item = elm_item, theta = score_ref, D = D, tcc = FALSE)$prob.cats
  prob_foc <-
    purrr::map(.x = score_foc, ~ {
      trace(elm_item = elm_item, theta = .x, D = D, tcc = FALSE)$prob.cats
    })

  # replace NA values into the missing data location
  extscore_ref[is.na(resp_ref)] <- NA
  for (g in 1:n.fg) {
    extscore_foc[[g]][is.na(resp_foc[[g]])] <- NA
  }

  # compute the raw residuals across all groups
  resid_ref <- resp_ref - extscore_ref
  resid_foc <-
    purrr::map2(.x = resp_foc, .y = extscore_foc, ~ {
      .x - .y
    })

  # compute the mean raw residuals (mrr) and mean squared residuals (msr) across all groups
  mrr_r <- colMeans(resid_ref, na.rm = TRUE)
  mrr_f <- purrr::map(.x = resid_foc, ~ {
    colMeans(.x, na.rm = TRUE)
  })
  msr_r <- colMeans(resid_ref^2, na.rm = TRUE)
  msr_f <- purrr::map(.x = resid_foc, ~ {
    colMeans(.x^2, na.rm = TRUE)
  })

  # compute the means, variances, and covariances of mrr and msr across all groups
  moments <- resid_moments_mg(
    p_ref = prob_ref, p_foc = prob_foc, n_ref = n_ref, n_foc = n_foc,
    resp_ref = resp_ref, resp_foc = resp_foc, cats = cats, n.fg = n.fg
  )
  mu_rdif <- moments$mu
  var_rdif <- moments$var
  covar_rdif <- moments$covar

  # create a design matrix for grdif_r and grdif_s
  mat.tmp1 <- cbind(rep(-1, n.fg))
  mat.tmp2 <- diag(n.fg)
  dmat_r <- dmat_s <- cbind(mat.tmp1, mat.tmp2)

  # create a design matrix for grdif_rs
  mat.tmp1 <- cbind(rep(c(-1, 0), n.fg))
  mat.tmp2 <- cbind(rep(c(0, -1), n.fg))
  mat.tmp3 <- diag(n.fg * 2)
  dmat_rs <- cbind(mat.tmp1, mat.tmp2, mat.tmp3)

  # compute the chi-square statistics
  chisq_r <- c()
  chisq_s <- c()
  chisq_rs <- c()
  for (i in 1:nitem) {
    if (i %in% all_miss) {
      chisq[i] <- NaN
    } else {
      # vectors of mrr and msr estimates
      mu_mrr_est <- c(mrr_r[i], purrr::map_dbl(.x = mrr_f, ~ {
        .x[i]
      }))
      mu_msr_est <- c(msr_r[i], purrr::map_dbl(.x = msr_f, ~ {
        .x[i]
      }))

      # vectors of mrr and msr parameters
      mu_mrr_par <- c(
        mu_rdif$ref[[1]][i],
        purrr::map_dbl(.x = purrr::map(.x = mu_rdif$foc, ~ {
          .x[[1]]
        }), ~ {
          .x[i]
        })
      )
      mu_msr_par <- c(
        mu_rdif$ref[[2]][i],
        purrr::map_dbl(.x = purrr::map(.x = mu_rdif$foc, ~ {
          .x[[2]]
        }), ~ {
          .x[i]
        })
      )

      # a vector including both mrr and msr estimates & parameters
      mu_est <- mu_par <- rep(NA, (n.fg + 1) * 2)
      mu_est[2 * (1:(n.fg + 1)) - 1] <- mu_mrr_est
      mu_est[2 * (1:(n.fg + 1))] <- mu_msr_est
      mu_par[2 * (1:(n.fg + 1)) - 1] <- mu_mrr_par
      mu_par[2 * (1:(n.fg + 1))] <- mu_msr_par

      # covariance matrices for grdif_r, grdif_s, and grdif_rs
      var_mrr <- c(
        var_rdif$ref[[1]][i],
        purrr::map_dbl(.x = purrr::map(.x = var_rdif$foc, ~ {
          .x[[1]]
        }), ~ {
          .x[i]
        })
      )
      var_msr <- c(
        var_rdif$ref[[2]][i],
        purrr::map_dbl(.x = purrr::map(.x = var_rdif$foc, ~ {
          .x[[2]]
        }), ~ {
          .x[i]
        })
      )
      covar <- c(covar_rdif$ref[i], purrr::map_dbl(.x = covar_rdif$foc, ~ {
        .x[i]
      }))
      cov_rdifr <- diag(var_mrr)
      cov_rdifs <- diag(var_msr)
      covlist <- list()
      for (g in 1:(n.fg + 1)) {
        # create a var-covariance matrix between mrr and msr
        cov_mat <- array(NA, c(2, 2))

        # replace NAs with the analytically computed covariance
        cov_mat[col(cov_mat) != row(cov_mat)] <- covar[g]

        # replace NAs with the analytically computed variances
        diag(cov_mat) <- c(var_mrr[g], var_msr[g])

        # add the cov-matrix
        covlist[[g]] <- cov_mat
      }
      cov_rdifrs <-
        Matrix::bdiag(covlist) %>%
        data.matrix()

      # multiply the design matrix with statistic vectors
      vec_r <- dmat_r %*% (cbind(mu_mrr_est) - cbind(mu_mrr_par))
      vec_s <- dmat_s %*% (cbind(mu_msr_est) - cbind(mu_msr_par))
      vec_rs <- dmat_rs %*% (cbind(mu_est) - cbind(mu_par))

      # multiply the design matrix with the covariance matrix
      dcov_r <- dmat_r %*% cov_rdifr %*% t(dmat_r)
      dcov_s <- dmat_s %*% cov_rdifs %*% t(dmat_s)
      dcov_rs <- dmat_rs %*% cov_rdifrs %*% t(dmat_rs)
      # dcov_r <- tcrossprod(x = (dmat_r %*% cov_rdifr), y = dmat_r)
      # dcov_s <- tcrossprod(x = (dmat_s %*% cov_rdifs), y = dmat_s)
      # dcov_rs <- tcrossprod(x = (dmat_rs %*% cov_rdifrs), y = dmat_rs)

      # find the inverse of the multiple covariance matrix
      inv_dcov_r <- suppressWarnings(tryCatch(
        {
          solve(dcov_r, tol = 1e-200)
        },
        error = function(e) {
          NULL
        }
      ))
      if (is.null(inv_dcov_r)) {
        inv_dcov_r <- suppressWarnings(tryCatch(
          {
            solve(dcov_r + 1e-15, tol = 1e-200)
          },
          error = function(e) {
            NULL
          }
        ))
        if (is.null(inv_dcov_r)) {
          inv_dcov_r <- suppressWarnings(tryCatch(
            {
              solve(dcov_r + 1e-10, tol = 1e-200)
            },
            error = function(e) {
              NULL
            }
          ))
        }
      }
      inv_dcov_s <- suppressWarnings(tryCatch(
        {
          solve(dcov_s, tol = 1e-200)
        },
        error = function(e) {
          NULL
        }
      ))
      if (is.null(inv_dcov_s)) {
        inv_dcov_s <- suppressWarnings(tryCatch(
          {
            solve(dcov_s + 1e-15, tol = 1e-200)
          },
          error = function(e) {
            NULL
          }
        ))
        if (is.null(inv_dcov_s)) {
          inv_dcov_s <- suppressWarnings(tryCatch(
            {
              solve(dcov_s + 1e-10, tol = 1e-200)
            },
            error = function(e) {
              NULL
            }
          ))
        }
      }
      inv_dcov_rs <- suppressWarnings(tryCatch(
        {
          solve(dcov_rs, tol = 1e-200)
        },
        error = function(e) {
          NULL
        }
      ))
      if (is.null(inv_dcov_rs)) {
        inv_dcov_rs <- suppressWarnings(tryCatch(
          {
            solve(dcov_rs + 1e-15, tol = 1e-200)
          },
          error = function(e) {
            NULL
          }
        ))
        if (is.null(inv_dcov_rs)) {
          inv_dcov_rs <- suppressWarnings(tryCatch(
            {
              solve(dcov_rs + 1e-10, tol = 1e-200)
            },
            error = function(e) {
              NULL
            }
          ))
        }
      }

      # compute the chi-square statistic
      chisq_r[i] <- as.numeric(t(vec_r) %*% inv_dcov_r %*% vec_r)
      chisq_s[i] <- as.numeric(t(vec_s) %*% inv_dcov_s %*% vec_s)
      chisq_rs[i] <- as.numeric(t(vec_rs) %*% inv_dcov_rs %*% vec_rs)
      # chisq_r[i] <- as.numeric(crossprod(x=vec_r, y=(inv_dcov_r %*% vec_r)))
      # chisq_s[i] <- as.numeric(crossprod(x=vec_s, y=(inv_dcov_s %*% vec_s)))
      # chisq_rs[i] <- as.numeric(crossprod(x=vec_rs, y=(inv_dcov_rs %*% vec_rs)))
    }
  }

  # calculate p-values for all three statistics
  p_grdifr <- round(stats::pchisq(chisq_r, df = n.fg, lower.tail = FALSE), 4)
  p_grdifs <- round(stats::pchisq(chisq_s, df = n.fg, lower.tail = FALSE), 4)
  p_grdifrs <- round(stats::pchisq(chisq_rs, df = n.fg * 2, lower.tail = FALSE), 4)

  # create a data.frame of sample sizes
  names(n_foc) <- paste0("n.foc", 1:n.fg)
  n_df <-
    data.frame(n.ref = n_ref, n_foc) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n.total = sum(dplyr::c_across(cols = dplyr::everything())))

  # create a data frame to contain the results
  stat_df <-
    data.frame(
      id = x$id,
      grdifr = round(chisq_r, 4),
      grdifs = round(chisq_s, 4),
      grdifrs = round(chisq_rs, 4),
      p.grdifr = p_grdifr,
      p.grdifs = p_grdifs,
      p.grdifrs = p_grdifrs,
      n_df,
      stringsAsFactors = FALSE
    )
  rownames(stat_df) <- NULL

  # find the flagged items
  dif_item_grdifr <- which(p_grdifr <= alpha)
  dif_item_grdifs <- which(p_grdifs <= alpha)
  dif_item_grdifrs <- which(p_grdifrs <= alpha)
  if (length(dif_item_grdifr) == 0) dif_item_grdifr <- NULL
  if (length(dif_item_grdifs) == 0) dif_item_grdifs <- NULL
  if (length(dif_item_grdifrs) == 0) dif_item_grdifrs <- NULL

  # data frame of the moments
  mu_df <- data.frame(mu_rdif)
  var_df <- data.frame(var_rdif)
  covar_df <- data.frame(covar_rdif)
  names(mu_df) <- c(
    paste0(c("mrr.", "msr."), "ref"),
    paste0(
      rep(c("mrr.", "msr."), n.fg),
      rep(paste0("foc", 1:n.fg), each = 2)
    )
  )
  names(var_df) <- c(
    paste0(c("mrr.", "msr."), "ref"),
    paste0(
      rep(c("mrr.", "msr."), n.fg),
      rep(paste0("foc", 1:n.fg), each = 2)
    )
  )
  names(covar_df) <- c("ref", paste0("foc", 1:n.fg))
  mmt_df <- list(mu = mu_df, var = var_df, covariance = covar_df)

  # summarize the results
  rst <- list(
    dif_stat = stat_df,
    dif_item = list(
      grdifr = dif_item_grdifr,
      grdifs = dif_item_grdifs,
      grdifrs = dif_item_grdifrs
    ),
    moments = mmt_df, alpha = alpha
  )

  # return the results
  rst
}

# This function computes the mean and variance of the IRT based residual statistics across all groups
# Also, it computes the covariance between mrr and msr across all groups
resid_moments_mg <- function(p_ref, p_foc, n_ref, n_foc, resp_ref, resp_foc, cats, n.fg) {
  # check the number of items
  nitem <- length(cats)

  # count the number of rows for the multiple group response data
  nrow_ref <- nrow(resp_ref)
  nrow_foc <- purrr::map(.x = resp_foc, ~ {
    nrow(.x)
  })

  # create an empty list of two matrices to contain the first moments of raw residuals and squared residuals
  mu_ref <- purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = nrow_ref, ncol = nitem))
  mu_foc <-
    purrr::map(
      .x = nrow_foc,
      ~ {
        purrr::map(.x = 1:2, .f = function(x) matrix(NA, nrow = .x, ncol = nitem))
      }
    )

  # create an empty list of two matrices to contain the second moments of raw residuals and squared residuals
  mu2_ref <- mu_ref
  mu2_foc <- mu_foc

  # create an empty matrix to contain the expectation of cube of raw residuals at each theta for all items
  mu_resid3_ref <- matrix(NA, nrow = nrow_ref, ncol = nitem)
  mu_resid3_foc <- purrr::map(.x = nrow_foc, ~ {
    matrix(NA, nrow = .x, ncol = nitem)
  })

  # compute the first and second moments of raw residual and squared residuals for all items
  for (i in 1:nitem) {
    # compute the expected residuals for each score category
    Emat_ref <- matrix(0:(cats[i] - 1), nrow = nrow_ref, ncol = cats[i], byrow = TRUE)
    Emat_foc <-
      purrr::map(
        .x = nrow_foc,
        ~ {
          matrix(0:(cats[i] - 1), nrow = .x, ncol = cats[i], byrow = TRUE)
        }
      )
    exp_resid_ref <-
      Emat_ref - matrix(rowSums(Emat_ref * p_ref[[i]]), nrow = nrow_ref, ncol = cats[i], byrow = FALSE)
    exp_resid_foc <-
      purrr::map(.x = 1:n.fg, ~ {
        Emat_foc[[.x]] - matrix(rowSums(Emat_foc[[.x]] * p_foc[[.x]][[i]]),
          nrow = nrow_foc[[.x]], ncol = cats[i], byrow = FALSE
        )
      })

    # replace NA values into the missing data location
    exp_resid_ref[is.na(resp_ref[, i]), ] <- NA
    for (g in 1:n.fg) {
      exp_resid_foc[[g]][is.na(resp_foc[[g]][, i]), ] <- NA
    }

    # compute the expected values of raw and squared residuals to be used to compute the mean and variance parameters
    value_ref <- vector("list", 2)
    value_ref[[1]] <- exp_resid_ref
    value_ref[[2]] <- exp_resid_ref^2
    resid3_ref <- exp_resid_ref^3
    value_foc <- purrr::map(.x = 1:n.fg, ~ {
      vector("list", 2)
    })
    resid3_foc <- vector("list", n.fg)
    for (g in 1:n.fg) {
      value_foc[[g]][[1]] <- exp_resid_foc[[g]]
      value_foc[[g]][[2]] <- exp_resid_foc[[g]]^2
      resid3_foc[[g]] <- exp_resid_foc[[g]]^3
    }

    # compute the first and second moments across all thetas for the reference and focal groups
    for (j in 1:2) {
      mu_ref[[j]][, i] <- rowSums(value_ref[[j]] * p_ref[[i]], na.rm = FALSE)
      mu2_ref[[j]][, i] <- rowSums((value_ref[[j]])^2 * p_ref[[i]], na.rm = FALSE)
      for (g in 1:n.fg) {
        mu_foc[[g]][[j]][, i] <- rowSums(value_foc[[g]][[j]] * p_foc[[g]][[i]], na.rm = FALSE)
        mu2_foc[[g]][[j]][, i] <- rowSums((value_foc[[g]][[j]])^2 * p_foc[[g]][[i]], na.rm = FALSE)
      }
    }

    # compute the expectation (first moment) of the cube of raw residuals across all thetas for the reference and focal groups
    mu_resid3_ref[, i] <- rowSums(resid3_ref * p_ref[[i]], na.rm = FALSE)
    for (g in 1:n.fg) {
      mu_resid3_foc[[g]][, i] <- rowSums(resid3_foc[[g]] * p_foc[[g]][[i]], na.rm = FALSE)
    }
  }

  # compute the variances across all thetas for the reference and focal groups
  var_ref <- purrr::map2(.x = mu2_ref, .y = mu_ref, .f = function(x, y) x - y^2)
  var_foc <-
    purrr::map(
      .x = 1:n.fg,
      ~ {
        purrr::map2(.x = mu2_foc[[.x]], .y = mu_foc[[.x]], .f = function(x, y) x - y^2)
      }
    )

  # compute the variances of MRR and MSR, respectively, across all groups
  # use V(aX - bY) = a^2 * V(X) + b^2 * V(Y)
  const_ref <- purrr::map(
    .x = var_ref,
    ~ {
      matrix((1 / n_ref^2), nrow = nrow(.x), ncol = ncol(.x), byrow = TRUE)
    }
  )
  const_foc <-
    purrr::map(
      .x = 1:n.fg,
      ~ {
        purrr::map(
          .x = var_foc[[.x]],
          .f = function(x) matrix((1 / n_foc[[.x]]^2), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
        )
      }
    )
  var_rdif_ref <- purrr::map2(.x = const_ref, .y = var_ref, ~ {
    colSums(.x * .y, na.rm = TRUE)
  })
  var_rdif_foc <-
    purrr::map(
      .x = 1:n.fg,
      ~ {
        purrr::map2(.x = const_foc[[.x]], .y = var_foc[[.x]], .f = function(x, y) colSums(x * y, na.rm = TRUE))
      }
    )
  names(var_rdif_ref) <- c("mrr", "msr")
  names(var_rdif_foc) <- paste0("foc", 1:n.fg)
  var_rdif_foc <-
    purrr::map(.x = var_rdif_foc, .f = function(x) {
      stats::setNames(x, c("mrr", "msr"))
    })

  # compute the means and variances of MRR and MSR, respectively, across all group
  mu_rdif_ref <- purrr::map(.x = mu_ref, ~ {
    colMeans(x = .x, na.rm = TRUE)
  })
  mu_rdif_foc <- purrr::map(.x = 1:n.fg, .f = function(x) {
    purrr::map(.x = mu_foc[[x]], ~ {
      colMeans(x = .x, na.rm = TRUE)
    })
  })
  names(mu_rdif_ref) <- c("mrr", "msr")
  names(mu_rdif_foc) <- paste0("foc", 1:n.fg)
  mu_rdif_foc <-
    purrr::map(
      .x = mu_rdif_foc,
      .f = function(x) {
        list(mrr = round(x[[1]], 10), msr = x[[2]])
      }
    )
  mu_rdif_ref[[1]] <- round(mu_rdif_ref[[1]], 10)

  # compute the covariance between mrr and msr across all groups
  covar_ref <- colSums(mu_resid3_ref, na.rm = TRUE) * (1 / n_ref^2)
  covar_foc <-
    purrr::map2(.x = mu_resid3_foc, .y = n_foc, ~ {
      colSums(.x, na.rm = TRUE) * (1 / .y^2)
    })
  names(covar_foc) <- paste0("foc", 1:n.fg)

  # return the results
  mu_rdif <- list(ref = mu_rdif_ref, foc = mu_rdif_foc)
  var_rdif <- list(ref = var_rdif_ref, foc = var_rdif_foc)
  covar_rdif <- list(ref = covar_ref, foc = covar_foc)
  rst <- list(mu = mu_rdif, var = var_rdif, covar = covar_rdif)
  rst
}
