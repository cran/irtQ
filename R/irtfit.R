#' Traditional IRT Item Fit Statistics
#'
#' This function computes traditional IRT item fit statistics, including the
#' \eqn{\chi^{2}} fit statistic (e.g., Bock, 1960; Yen, 1981),
#' the log-likelihood ratio \eqn{\chi^{2}} fit statistic (\eqn{G^{2}}; McKinley
#' & Mills, 1985), and the infit and outfit statistics (Ames et al., 2015). It
#' also returns contingency tables used to compute the \eqn{\chi^{2}} and
#' \eqn{G^{2}} statistics.
#'
#' @inheritParams est_score
#' @param x A data frame containing item metadata (e.g., item parameters,
#'   number of categories, IRT model types, etc.); or an object of class
#'   `est_irt` obtained from [irtQ::est_irt()], or `est_item` from
#'   [irtQ::est_item()].
#'
#'   See [irtQ::est_irt()] or [irtQ::simdat()] for more details about the item
#'   metadata. This data frame can be easily created using the
#'   [irtQ::shape_df()] function.
#' @param score A numeric vector containing examinees' ability estimates
#'   (theta values).
#' @param group.method A character string specifying the method used to group
#'   examinees along the ability scale when computing the \eqn{\chi^{2}} and
#'   \eqn{G^{2}} fit statistics. Available options are:
#'    - `"equal.width"`: Divides the ability scale into intervals of equal width.
#'    - `"equal.freq"`: Divides the examinees into groups with (approximately)
#'    equal numbers of examinees.
#'
#'   Note that `"equal.freq"` does not guarantee exactly equal group sizes
#'   due to ties in ability estimates. Default is `"equal.width"`.
#'   The number of groups and the range of the ability scale are controlled by
#'   the `n.width` and `range.score` arguments, respectively.
#' @param n.width An integer specifying the number of intervals (groups) into which
#' the ability scale is divided for computing the fit statistics. Default is 10.
#' @param loc.theta A character string indicating the point on the ability scale
#'  at which the expected category probabilities are calculated for each group.
#'  Available options are:
#'  - `"average"`: Uses the average ability estimate of examinees within each group.
#'  - `"middle"`: Uses the midpoint of each group's ability interval.
#'
#'  Default is `"average"`.
#' @param range.score A numeric vector of length two specifying the lower and upper
#'   bounds of the ability scale. Ability estimates below the lower bound or above the
#'   upper bound are truncated to the respective bound. If `NULL`, the observed minimum
#'   and maximum of the `score` vector are used. Note that this range restriction is
#'   independent of the grouping method specified in `group.method`.
#'   Default is `NULL`.
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha}) for
#'   the hypothesis tests of the \eqn{\chi^{2}} and \eqn{G^{2}} item fit statistics.
#'   Default is `0.05`.
#' @param overSR A numeric threshold used to identify ability groups (intervals)
#'   whose standardized residuals exceed the specified value. This is used to
#'   compute the proportion of misfitting groups per item. Default is 2.
#' @param min.collapse An integer specifying the minimum expected frequency required
#'   for a cell in the contingency table. Neighboring groups will be merged if any
#'   expected cell frequency falls below this threshold when computing the
#'   \eqn{\chi^{2}} and \eqn{G^{2}} statistics. Default is 1.
#' @param pcm.loc A vector of integers indicating the locations (indices) of
#'   partial credit model (PCM) items for which slope parameters are fixed.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' To compute the \eqn{\chi^2} and \eqn{G^2} item fit statistics, the `group.method`
#' argument determines how the ability scale is divided into groups:
#' - `"equal.width"`: Examinees are grouped based on intervals of equal width a
#'  long the ability scale.
#' - `"equal.freq"`: Examinees are grouped such that each group contains
#'  (approximately) the same number of individuals.
#'
#' Note that `"equal.freq"` does not guarantee *exactly* equal frequencies across
#' all groups, since grouping is based on quantiles.
#'
#' When dividing the ability scale into intervals to compute the \eqn{\chi^2}
#' and \eqn{G^2} fit statistics, the intervals should be:
#' - **Wide enough** to ensure that each group contains a sufficient number of
#'  examinees (to avoid unstable estimates),
#' - **Narrow enough** to ensure that examinees within each group are relatively
#'  homogeneous in ability (Hambleton et al., 1991).
#'
#' If you want to divide the ability scale into a number of groups other than
#' the default of 10, specify the desired number using the `n.width` argument.
#' For reference:
#' - Yen (1981) used 10 fixed-width groups,
#' - Bock (1960) allowed for flexibility in the number of groups.
#'
#' Regarding degrees of freedom (*df*):
#' - The \eqn{\chi^2} statistic is approximately chi-square distributed with
#'   degrees of freedom equal to the number of ability groups minus the number
#'   of item parameters (Ames et al., 2015).
#' - The \eqn{G^2} statistic is approximately chi-square distributed with
#'   degrees of freedom equal to the number of ability groups (Ames et al., 2015;
#'   Muraki & Bock, 2003).
#'
#'   Note that if `"DRM"` is specified for an item in the item metadata set,
#'   the item is treated as a `"3PLM"` when computing the degrees of freedom
#'   for the \eqn{\chi^{2}} fit statistic.
#'
#' Note that infit and outfit statistics should be interpreted with caution when
#' applied to non-Rasch models. The returned object—particularly the contingency
#' tables—can be passed to [irtQ::plot.irtfit()] to generate raw and standardized
#' residual plots (Hambleton et al., 1991).
#'
#' @return This function returns an object of class `irtfit`, which includes
#' the following components:
#'
#'   \item{fit_stat}{A data frame containing the results of three IRT item fit statistics—
#'   \eqn{\chi^{2}}, \eqn{G^{2}}, infit, and outfit—for all evaluated items.
#'   Each row corresponds to one item, and the columns include:
#'   the item ID; \eqn{\chi^{2}} statistic; \eqn{G^{2}} statistic;
#'   degrees of freedom for \eqn{\chi^{2}} and \eqn{G^{2}};
#'   critical values and p-values for both statistics;
#'   outfit and infit values;
#'   the number of examinees used to compute these statistics;
#'   and the proportion of ability groups (prior to cell collapsing) that have
#'   standardized residuals greater than the threshold specified in the `overSR` argument.}
#'
#'   \item{contingency.fitstat}{A list of contingency tables used to compute the
#'   \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics for all items. Note that the
#'   cell-collapsing strategy is applied to these tables to ensure sufficient
#'   expected frequencies.}
#'
#'   \item{contingency.plot}{A list of contingency tables used to generate raw
#'   and standardized residual plots (Hambleton et al., 1991) via the
#'   [irtQ::plot.irtfit()]. Note that these tables are based on the original,
#'   uncollapsed groupings.}
#'
#'   \item{individual.info}{A list of data frames containing individual residuals
#'   and corresponding variance values. This information is used to compute infit
#'   and outfit statistics.}
#'
#'   \item{item_df}{A data frame containing the item metadata provided in the
#'   argument `x`.}
#'
#'   \item{ancillary}{A list of ancillary information used during the item fit
#'   analysis.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::plot.irtfit()], [irtQ::shape_df()], [irtQ::est_irt()],
#' [irtQ::est_item()]
#'
#' @references Ames, A. J., & Penfield, R. D. (2015). An NCME Instructional
#' Module on Item-Fit Statistics for Item Response Theory Models.
#' *Educational Measurement: Issues and Practice, 34*(3), 39-48.
#'
#' Bock, R.D. (1960), *Methods and applications of optimal scaling*. Chapel
#' Hill, NC: L.L. Thurstone Psychometric Laboratory.
#'
#' Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).*Fundamentals of
#' item response theory*. Newbury Park, CA: Sage.
#'
#' McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit
#' statistics.
#' *Applied Psychological Measurement, 9*, 49-57.
#'
#' Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test
#' scoring for rating scale data (Computer Software). Chicago, IL: Scientific
#' Software International. URL http://www.ssicentral.com
#'
#' Wells, C. S., & Bolt, D. M. (2008). Investigation of a nonparametric
#' procedure for assessing goodness-of-fit in item response theory. *Applied
#' Measurement in Education, 21*(1), 22-40.
#'
#' Yen, W. M. (1981). Using simulation results to choose a latent trait model.
#' *Applied Psychological Measurement, 5*, 245-262.
#'
#' @examples
#' \donttest{
#' ## Example 1
#' ## Use the simulated CAT data
#' # Identify items with more than 10,000 responses
#' over10000 <- which(colSums(simCAT_MX$res.dat, na.rm = TRUE) > 10000)
#'
#' # Select items with more than 10,000 responses
#' x <- simCAT_MX$item.prm[over10000, ]
#'
#' # Extract response data for the selected items
#' data <- simCAT_MX$res.dat[, over10000]
#'
#' # Extract examinees' ability estimates
#' score <- simCAT_MX$score
#'
#' # Compute item fit statistics
#' fit1 <- irtfit(
#'   x = x, score = score, data = data, group.method = "equal.width",
#'   n.width = 10, loc.theta = "average", range.score = NULL, D = 1, alpha = 0.05,
#'   missing = NA, overSR = 2
#' )
#'
#' # View the fit statistics
#' fit1$fit_stat
#'
#' # View the contingency tables used to compute fit statistics
#' fit1$contingency.fitstat
#'
#'
#' ## Example 2
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select the first two dichotomous items and the last polytomous item
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # Generate ability values from a standard normal distribution
#' set.seed(10)
#' score <- rnorm(1000, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' # Compute item fit statistics
#' fit2 <- irtfit(
#'   x = x, score = score, data = data, group.method = "equal.freq",
#'   n.width = 11, loc.theta = "average", range.score = c(-4, 4), D = 1, alpha = 0.05
#' )
#'
#' # View the fit statistics
#' fit2$fit_stat
#'
#' # View the contingency tables used to compute fit statistics
#' fit2$contingency.fitstat
#'
#' # Plot raw and standardized residuals for the first item (dichotomous)
#' plot(x = fit2, item.loc = 1, type = "both", ci.method = "wald",
#'      show.table = TRUE, ylim.sr.adjust = TRUE)
#'
#' # Plot raw and standardized residuals for the third item (polytomous)
#' plot(x = fit2, item.loc = 3, type = "both", ci.method = "wald",
#'      show.table = FALSE, ylim.sr.adjust = TRUE)
#' }
#'
#' @export
irtfit <- function(x, ...) UseMethod("irtfit")

#' @describeIn irtfit Default method for computing traditional IRT item fit
#' statistics using a data frame `x` that contains item metadata.
#' @import dplyr
#' @export
#'
irtfit.default <- function(x,
                           score,
                           data,
                           group.method = c("equal.width", "equal.freq"),
                           n.width = 10,
                           loc.theta = "average",
                           range.score = NULL,
                           D = 1,
                           alpha = 0.05,
                           missing = NA,
                           overSR = 2,
                           min.collapse = 1,
                           pcm.loc = NULL,
                           ...) {

  # match.call
  cl <- match.call()

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # create a vector of PCM item indicators
  pcm.lg <- logical(nrow(x))
  pcm.lg[pcm.loc] <- TRUE

  # transform scores to a vector form
  if (is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check if there are items which have zero or one response frequency
  n.score <- Rfast::colsums(!is.na(data))
  if (all(n.score %in% c(0L, 1L))) {
    stop("Every item has frequency of zero or one for the item response data. Each item must have more than two item responses.", call. = FALSE)
  }

  if (any(n.score %in% c(0L, 1L))) {
    del_item <- which(n.score %in% c(0L, 1L))

    # delete the items which have no frequency of scores from the data set
    x <- x[-del_item, ]
    data <- data[, -del_item]
    pcm.lg <- pcm.lg[-del_item]

    # warning message
    memo <- paste0(
      paste0("item ", del_item, collapse = ", "),
      " is/are excluded in the analysis. Because the item(s) has/have frequency of zero or one for the item response data."
    )
    warning(memo, call. = FALSE)
  }

  # restrict the range of scores if required
  if (!is.null(range.score)) {
    score <- ifelse(score < range.score[1], range.score[1], score)
    score <- ifelse(score > range.score[2], range.score[2], score)
  } else {
    tmp.val <- max(ceiling(abs(range(score, na.rm = TRUE))))
    range.score <- c(-tmp.val, tmp.val)
  }

  # compute item fit statistics and obtain contingency tables across all items
  fits <-
    purrr::map(1:nrow(x), .f = function(i) {
      itemfit(
        x_item = x[i, ], score = score, resp = data[, i], group.method = group.method,
        n.width = n.width, loc.theta = loc.theta, D = D, alpha = alpha, overSR = overSR,
        min.collapse = min.collapse, is.pcm = pcm.lg[i]
      )
    })

  # extract fit statistics
  fit_stat <-
    purrr::map(fits, .f = function(i) i$fit.stats) %>%
    do.call(what = "rbind")
  fit_stat <- data.frame(id = x$id, fit_stat)

  # extract the contingency tables used to compute the fit statistics
  contingency.fitstat <-
    purrr::map(fits, .f = function(i) i$contingency.fitstat)
  names(contingency.fitstat) <- x$id

  # extract the contingency tables to be used to draw residual plots
  contingency.plot <-
    purrr::map(fits, .f = function(i) i$contingency.plot)
  names(contingency.plot) <- x$id

  # extract the individual residuals and variances
  individual.info <-
    purrr::map(fits, .f = function(i) i$individual.info)
  names(individual.info) <- x$id

  # return results
  rst <- list(
    fit_stat = fit_stat, contingency.fitstat = contingency.fitstat,
    contingency.plot = contingency.plot,
    item_df = x, individual.info = individual.info,
    ancillary = list(range.score = range.score, alpha = alpha, overSR = overSR, scale.D = D)
  )
  class(rst) <- "irtfit"
  rst$call <- cl

  rst
}


#' @describeIn irtfit An object created by the function [irtQ::est_item()].
#' @import dplyr
#' @export
#'
irtfit.est_item <- function(x,
                            group.method = c("equal.width", "equal.freq"),
                            n.width = 10,
                            loc.theta = "average",
                            range.score = NULL,
                            alpha = 0.05,
                            missing = NA,
                            overSR = 2,
                            min.collapse = 1,
                            pcm.loc = NULL,
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

  # create a vector of PCM item indicators
  pcm.lg <- logical(nrow(x))
  pcm.lg[pcm.loc] <- TRUE

  # transform scores to a vector form
  if (is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check if there are items which have zero or one response frequency
  n.score <- Rfast::colsums(!is.na(data))
  if (all(n.score %in% c(0L, 1L))) {
    stop("Every item has frequency of zero or one for the item response data. Each item must have more than two item responses.", call. = FALSE)
  }

  if (any(n.score %in% c(0L, 1L))) {
    del_item <- which(n.score %in% c(0L, 1L))

    # delete the items which have no frequency of scores from the data set
    x <- x[-del_item, ]
    data <- data[, -del_item]
    pcm.lg <- pcm.lg[-del_item]

    # warning message
    memo <- paste0(
      paste0("item ", del_item, collapse = ", "),
      " is/are excluded in the analysis. Because the item(s) has/have frequency of zero or one for the item response data."
    )
    warning(memo, call. = FALSE)
  }

  # restrict the range of scores if required
  if (!is.null(range.score)) {
    score <- ifelse(score < range.score[1], range.score[1], score)
    score <- ifelse(score > range.score[2], range.score[2], score)
  } else {
    tmp.val <- max(ceiling(abs(range(score, na.rm = TRUE))))
    range.score <- c(-tmp.val, tmp.val)
  }

  # compute item fit statistics and obtain contingency tables across all items
  fits <-
    purrr::map(1:nrow(x), .f = function(i) {
      itemfit(
        x_item = x[i, ], score = score, resp = data[, i], group.method = group.method,
        n.width = n.width, loc.theta = loc.theta, D = D, alpha = alpha, overSR = overSR,
        min.collapse = min.collapse, is.pcm = pcm.lg[i]
      )
    })

  # extract fit statistics
  fit_stat <-
    purrr::map(fits, .f = function(i) i$fit.stats) %>%
    do.call(what = "rbind")
  fit_stat <- data.frame(id = x$id, fit_stat)

  # extract the contingency tables used to compute the fit statistics
  contingency.fitstat <-
    purrr::map(fits, .f = function(i) i$contingency.fitstat)
  names(contingency.fitstat) <- x$id

  # extract the contingency tables to be used to draw residual plots
  contingency.plot <-
    purrr::map(fits, .f = function(i) i$contingency.plot)
  names(contingency.plot) <- x$id

  # extract the individual residuals and variances
  individual.info <-
    purrr::map(fits, .f = function(i) i$individual.info)
  names(individual.info) <- x$id

  # return results
  rst <- list(
    fit_stat = fit_stat, contingency.fitstat = contingency.fitstat,
    contingency.plot = contingency.plot,
    item_df = x, individual.info = individual.info,
    ancillary = list(range.score = range.score, alpha = alpha, overSR = overSR, scale.D = D)
  )
  class(rst) <- "irtfit"
  rst$call <- cl

  rst
}


#' @describeIn irtfit An object created by the function [irtQ::est_irt()].
#' @import dplyr
#' @export
#'
irtfit.est_irt <- function(x,
                           score,
                           group.method = c("equal.width", "equal.freq"),
                           n.width = 10,
                           loc.theta = "average",
                           range.score = NULL,
                           alpha = 0.05,
                           missing = NA,
                           overSR = 2,
                           min.collapse = 1,
                           pcm.loc = NULL,
                           ...) {

  # match.call
  cl <- match.call()

  # extract information from an object
  data <- x$data
  D <- x$scale.D
  x <- x$par.est

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # create a vector of PCM item indicators
  pcm.lg <- logical(nrow(x))
  pcm.lg[pcm.loc] <- TRUE

  # transform scores to a vector form
  if (is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check if there are items which have zero or one response frequency
  n.score <- Rfast::colsums(!is.na(data))
  if (all(n.score %in% c(0L, 1L))) {
    stop("Every item has frequency of zero or one for the item response data. Each item must have more than two item responses.", call. = FALSE)
  }

  if (any(n.score %in% c(0L, 1L))) {
    del_item <- which(n.score %in% c(0L, 1L))

    # delete the items which have no frequency of scores from the data set
    x <- x[-del_item, ]
    data <- data[, -del_item]
    pcm.lg <- pcm.lg[-del_item]

    # warning message
    memo <- paste0(
      paste0("item ", del_item, collapse = ", "),
      " is/are excluded in the analysis. Because the item(s) has/have frequency of zero or one for the item response data."
    )
    warning(memo, call. = FALSE)
  }

  # restrict the range of scores if required
  if (!is.null(range.score)) {
    score <- ifelse(score < range.score[1], range.score[1], score)
    score <- ifelse(score > range.score[2], range.score[2], score)
  } else {
    tmp.val <- max(ceiling(abs(range(score, na.rm = TRUE))))
    range.score <- c(-tmp.val, tmp.val)
  }

  # compute item fit statistics and obtain contingency tables across all items
  fits <-
    purrr::map(1:nrow(x), .f = function(i) {
      itemfit(
        x_item = x[i, ], score = score, resp = data[, i], group.method = group.method,
        n.width = n.width, loc.theta = loc.theta, D = D, alpha = alpha, overSR = overSR,
        min.collapse = min.collapse, is.pcm = pcm.lg[i]
      )
    })

  # extract fit statistics
  fit_stat <-
    purrr::map(fits, .f = function(i) i$fit.stats) %>%
    do.call(what = "rbind")
  fit_stat <- data.frame(id = x$id, fit_stat)

  # extract the contingency tables used to compute the fit statistics
  contingency.fitstat <-
    purrr::map(fits, .f = function(i) i$contingency.fitstat)
  names(contingency.fitstat) <- x$id

  # extract the contingency tables to be used to draw residual plots
  contingency.plot <-
    purrr::map(fits, .f = function(i) i$contingency.plot)
  names(contingency.plot) <- x$id

  # extract the individual residuals and variances
  individual.info <-
    purrr::map(fits, .f = function(i) i$individual.info)
  names(individual.info) <- x$id

  # return results
  rst <- list(
    fit_stat = fit_stat, contingency.fitstat = contingency.fitstat,
    contingency.plot = contingency.plot,
    item_df = x, individual.info = individual.info,
    ancillary = list(range.score = range.score, alpha = alpha, overSR = overSR, scale.D = D)
  )
  class(rst) <- "irtfit"
  rst$call <- cl

  rst
}



#' @importFrom janitor adorn_totals adorn_percentages
#' @importFrom tibble rownames_to_column remove_rownames column_to_rownames
#' @import dplyr
itemfit <- function(x_item, score, resp, group.method = c("equal.width", "equal.freq"),
                    n.width = 10, loc.theta = "average", D = 1, alpha = 0.05, overSR = 2, min.collapse = 1,
                    is.pcm = FALSE) {
  # break down the item metadata into several elements
  elm_item <- breakdown(x_item)

  # number of score categories
  cats <- elm_item$cats

  # delete missing responses and corresponding thetas
  na_lg <- is.na(resp)
  resp <- resp[!na_lg]
  score <- score[!na_lg]

  # assign factor levels to the item responses
  resp <- factor(resp, levels = (seq_len(cats) - 1))

  # compute cut scores to divide score groups
  group.method <- match.arg(group.method)
  cutscore <- switch(group.method,
    equal.width = seq(from = min(score), to = max(score), length.out = n.width + 1),
    equal.freq = stats::quantile(score, probs = seq(0, 1, length.out = n.width + 1), type = 9)
  )

  # when there are the same cutscores,
  # use only unique cutscores
  cutscore <- unique(cutscore)

  # assign score group variable to each score
  intv <- cut(score, breaks = cutscore, right = FALSE, include.lowest = TRUE, dig.lab = 7)

  # create a contingency table for the frequencies of score points
  obs.freq <-
    as.data.frame.matrix(table(intv, resp)) %>%
    tibble::rownames_to_column(var = "interval") %>%
    janitor::adorn_totals(where = "col", name = "total")
  delrow.lg <- obs.freq$total == 0
  obs.freq <- obs.freq[!delrow.lg, ]

  # create a contingency table for the category proportions
  obs.prop <-
    obs.freq %>%
    janitor::adorn_percentages(denominator = "row")

  # find a theta point for each score group
  loc.theta <- tolower(loc.theta)
  if (loc.theta == "middle") {
    theta <- purrr::map_dbl(.x = 2:length(cutscore), .f = function(i) mean(c(cutscore[i - 1], cutscore[i])))
    theta <- theta[!delrow.lg]
  } else {
    theta <-
      data.frame(intv, score) %>%
      dplyr::group_by(.data$intv) %>%
      dplyr::summarize(ave = mean(.data$score), .groups = "drop") %>%
      dplyr::pull(2)
  }

  # compute expected probabilities of endorsing an answer to each score category
  exp.prob <-
    trace(elm_item = elm_item, theta = theta, D = D, tcc = FALSE)$prob.cats[[1]] %>%
    data.frame()
  colnames(exp.prob) <- 0:(cats - 1)

  ## -------------------------------------------------------------------------
  # find a z-score corresponding to significance level
  zscore <- stats::qnorm(1 - alpha)

  # compute raw residuals (rr)
  rr <- obs.prop[2:(cats + 1)] - exp.prob

  # compute the standard errors (se)
  se <- sqrt((exp.prob * (1 - exp.prob)) / obs.freq$total)

  # standardize the residuals (sr)
  sr <- rr / se

  # compute a proportion of groups (or intervals) that have the standardized
  # residuals greater than a specified criterion
  over_sr <- sum(abs(sr) > overSR)
  over_sr_prop <- round(over_sr / (cats * nrow(sr)), 3)

  # create a full contingency table to draw IRT residual plots
  ctg_tb <-
    data.frame(
      point = theta, obs.freq = obs.freq, obs.prop = obs.prop[, 2:(cats + 1)],
      exp.prob = exp.prob, raw.rsd = rr, se = se, std.rsd = sr
    ) %>%
    dplyr::relocate("obs.freq.interval", .before = "point") %>%
    dplyr::relocate("obs.freq.total", .after = "point") %>%
    dplyr::rename("interval" = "obs.freq.interval", "total" = "obs.freq.total")

  ## ------------------------------------------------------------------------------
  # collapsing the contingency tables to compute the chi-square fit statistics
  # check the number of expected frequency for all cells
  exp.freq <- exp.prob * obs.freq$total

  # collapse the expected and observed frequency tables
  ftable_info <-
    data.frame(exp.freq = exp.freq, obs.freq = obs.freq) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "obs.freq.interval")
  for (i in 1:cats) {
    ftable_info <- collapse_ftable(x = ftable_info, col = i, min.collapse = min.collapse)
  }

  # new contingency tables after collapsing
  exp.freq.cp <- dplyr::select(ftable_info, dplyr::contains("exp.freq"))
  obs.freq.cp <- dplyr::select(ftable_info, dplyr::contains("obs.freq"))
  freq.tot.cp <- obs.freq.cp$obs.freq.total
  exp.prob.cp <- exp.freq.cp / freq.tot.cp
  obs.prop.cp <- (obs.freq.cp / freq.tot.cp)[, 1:cats]
  colnames(exp.prob.cp) <- paste0("exp.prob.", 0:(cats - 1))
  colnames(obs.prop.cp) <- paste0("obs.prop.", 0:(cats - 1))

  ## -------------------------------------------------------------------------
  # create a contingency table to compute the item fit statistics
  # first, compute raw residuals
  rr.cp <- obs.prop.cp - exp.prob.cp
  colnames(rr.cp) <- paste0("raw.rsd.", 0:(cats - 1))

  # create a full contingency table for chi-square fit statistic
  ctg_tb.cp <-
    data.frame(obs.freq.cp, exp.freq.cp, obs.prop.cp, exp.prob.cp, rr.cp) %>%
    dplyr::relocate("obs.freq.total", .before = "obs.freq.0") %>%
    dplyr::rename("total" = "obs.freq.total")
  rownames(ctg_tb.cp) <- 1:nrow(ctg_tb.cp)

  # compute the chi-square statistic (X2)
  x2 <- sum(freq.tot.cp * (rr.cp^2 / exp.prob.cp), na.rm = TRUE)

  # compute the likelihood ratio chi-square fit statistic (G2)
  g2 <- 2 * sum(obs.freq.cp[, 1:cats] * log(obs.prop.cp / exp.prob.cp), na.rm = TRUE)

  # find the number of parameters for each item
  model <- x_item$model
  if (is.pcm) model <- "PCM"
  count_prm <- NA
  count_prm[model %in% "1PLM"] <- 1
  count_prm[model %in% "2PLM"] <- 2
  count_prm[model %in% c("3PLM", "DRM")] <- 3
  count_prm[model %in% "PCM"] <- x_item[model %in% "PCM", 2] - 1
  count_prm[model %in% "GPCM"] <- x_item[model %in% "GPCM", 2]
  count_prm[model %in% "GRM"] <- x_item[model %in% "GRM", 2]

  # find a critical value and compute the p values
  df.x2 <- nrow(exp.freq.cp) * (ncol(exp.freq.cp) - 1) - count_prm
  df.g2 <- nrow(exp.freq.cp) * (ncol(exp.freq.cp) - 1)
  crtval.x2 <- stats::qchisq(1 - alpha, df = df.x2, lower.tail = TRUE)
  crtval.g2 <- stats::qchisq(1 - alpha, df = df.g2, lower.tail = TRUE)
  pval.x2 <- 1 - stats::pchisq(x2, df = df.x2, lower.tail = TRUE)
  pval.g2 <- 1 - stats::pchisq(g2, df = df.g2, lower.tail = TRUE)

  ## ------------------------------------------------------------------------------
  # infit & outfit
  # individual expected probabilities for each score category
  indiv_exp.prob <- trace(
    elm_item = elm_item, theta = score, D = D,
    tcc = FALSE
  )$prob.cats[[1]]

  # individual observed proportion for each score category
  n.resp <- length(resp)
  indiv_obs.prob <-
    table(1:n.resp, resp) %>%
    as.data.frame.matrix()

  # a matrix of the expected scores for each category
  Emat <- matrix(0:(cats - 1), nrow(indiv_exp.prob), ncol(indiv_exp.prob), byrow = TRUE)

  # residuals
  resid <- rowSums(indiv_obs.prob * Emat) - rowSums(Emat * indiv_exp.prob)

  # variance
  Var <- rowSums((Emat - rowSums(Emat * indiv_exp.prob))^2 * indiv_exp.prob)

  # compute outfit & infit
  outfit <- sum(resid^2 / Var) / n.resp
  infit <- sum(resid^2) / sum(Var)

  ## ------------------------------------------------------------------------------
  # summary of fit statistics
  fitstats <- data.frame(
    X2 = round(x2, 3), G2 = round(g2, 3), df.X2 = df.x2, df.G2 = df.g2,
    crit.val.X2 = round(crtval.x2, 2), crit.val.G2 = round(crtval.g2, 2),
    p.X2 = round(pval.x2, 3), p.G2 = round(pval.g2, 3),
    outfit = round(outfit, 3), infit = round(infit, 3),
    N = n.resp, overSR.prop = over_sr_prop
  )

  ## ------------------------------------------------------------------------------
  # return results
  list(
    fit.stats = fitstats, contingency.fitstat = ctg_tb.cp, contingency.plot = ctg_tb,
    individual.info = data.frame(resid = resid, Var = Var)
  )
}
