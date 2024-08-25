#' Pseudo-count D2 method
#'
#' @description This function calculates the Pseudo-count \eqn{D^{2}} statistic
#' to evaluate item parameter drift, as described by Cappaert et al. (2018) and Stone (2000).
#' The Pseudo-count \eqn{D^{2}} statistic is designed to detect item parameter drift efficiently
#' without requiring item recalibration, making it especially valuable in computerized adaptive
#' testing (CAT) environments. This method compares observed and expected response frequencies
#' across quadrature points, which represent latent ability levels. The expected frequencies are
#' computed using the posterior distribution of each examinee's ability (Stone, 2000), providing
#' a robust and sensitive measure of item parameter drift, ensuring the stability and accuracy
#' of the test over time.
#'
#' @param x A data frame containing the metadata for the item bank, which includes
#' important information for each item such as the number of score categories and the
#' IRT model applied. See \code{\link{est_irt}}, \code{\link{irtfit}},
#' \code{\link{info}} or \code{\link{simdat}} for more detail about the item metadata.
#' @param data A matrix containing examinees' response data of the items in the
#' argument \code{x}. A row and column indicate the examinees and items, respectively.
#' @param D A scaling factor in IRT models to make the logistic function as close as
#' possible to the normal ogive function (if set to 1.7). Default is 1.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param Quadrature A numeric vector of two components specifying the number of quadrature
#' points (in the first component) and the symmetric minimum and maximum values of these points
#' (in the second component). For example, a vector of c(49, 6) indicates 49 rectangular
#' quadrature points over -6 and 6.
#' @param weights A two-column matrix or data frame containing the quadrature points
#' (in the first column) and the corresponding weights (in the second column) of the latent
#' variable prior distribution. If not NULL, the scale of the latent ability distribution will
#' be will be fixed to the scale of the provided quadrature points and weights.
#' The weights and quadrature points can be easily obtained using the function \code{\link{gen.weight}}.
#' If NULL, a normal prior density is used based on the information provided in the arguments
#' of \code{Quadrature}, \code{group.mean}, and \code{group.var}). Default is NULL.
#' @param group.mean A numeric value to set the mean of latent variable prior distribution
#' when \code{weights = NULL}. Default is 0.
#' @param group.var A positive numeric value to set the variance of latent variable
#' prior distribution when \code{weights = NULL}. Default is 1.
#'
#' @details The Pseudo-count \eqn{D^{2}} values are calculated by summing the
#' weighted squared differences between the observed and expected frequencies
#' for each score category across all items. The expected frequencies are determined
#' using the posterior distribution of each examinee's ability (Stone, 2000).
#'
#' The Pseudo-count \eqn{D^{2}} statistic is calculated as:
#' \deqn{
#' Pseudo-count D^{2} = \sum_{k=1}^{Q} \left( \frac{r_{0k} + r_{1k}}{N}\right) \left( \frac{r_{1k}}{r_{0k} + r_{1k}} - E_{1k} \right)^2
#' }
#'
#' where \eqn{r_{0k}} and \eqn{r_{1k}} are the pseudo-counts for the incorrect and correct responses
#' at each ability level \eqn{k}, \eqn{E_{1k}} is the expected proportion of correct responses at each ability level \eqn{k},
#' calculated using item parameters from the item bank, and \eqn{N} is the total count of examinees
#' who received each item
#'
#' The \code{\link{pcd2}} function is designed to be flexible and allows for detailed control over
#' the computation process through its various arguments:
#'
#' \describe{
#'   \item{x}{The metadata should include key information such as the number of score
#'   categories for each item, the IRT model applied (e.g., "1PLM", "2PLM", "3PLM"),
#'   and the item parameters (e.g., discrimination, difficulty, guessing). This
#'   data frame is crucial because it defines the structure of the items and
#'   the parameters used in the calculation of expected frequencies.}
#'
#'   \item{data}{The response matrix should be preprocessed to ensure that missing values
#'   are correctly coded (using the `missing` argument if needed). The matrix is
#'   the foundation for calculating both observed and expected frequencies for
#'   each item and score category. Ensure that the number of items in this matrix
#'   matches the number specified in the `x` argument.}
#'
#'   \item{Quadrature}{The quadrature points are used to approximate the latent ability
#'   distribution, which is essential for accurately estimating the posterior
#'   distribution of abilities. Adjust the number and range of quadrature points
#'   based on the precision required and the characteristics of the examinee population.}
#'
#'   \item{weights}{This argument allows you to provide a custom two-column matrix or data
#'   frame where the first column contains quadrature points and the second column contains
#'   the corresponding weights. This provides flexibility in defining the latent ability
#'   distribution. If `weights` is set to `NULL`, the function generates a normal prior
#'   distribution based on the values provided in `Quadrature`, `group.mean`, and
#'   `group.var`. Use this argument when you have a specific latent ability distribution
#'   you wish to apply, such as one derived from empirical data or an alternative theoretical
#'   distribution.}
#'
#'   \item{group.mean, group.var}{These numeric values define the mean and variance
#'   of the latent ability distribution when `weights` is not provided. The default values
#'   are `group.mean = 0` and `group.var = 1`, which assume a standard normal distribution.
#'   These values are used to generate quadrature points and weights internally if `weights`
#'   is `NULL`. Adjust `group.mean` and `group.var` to reflect the expected distribution
#'   of abilities in your examinee population, particularly if you suspect that the latent
#'   trait distribution deviates from normality.}
#'
#' }
#'
#' When setting these arguments, consider the specific characteristics of your item
#' bank, the distribution of abilities in your examinee population, and the computational
#' precision required. Properly configured, the function provides robust and accurate
#' Pseudo-count \eqn{D^{2}} statistics, enabling effective monitoring of item parameter
#' drift in CAT or other IRT-based testing environments.
#'
#' @return A data frame containing the Pseudo-count \eqn{D^{2}} statistic for each item
#' in the analysis. The data frame includes the following columns:
#'  \item{id}{The unique identifier for each item, corresponding to the item IDs provided in
#'   the \code{x} argument.}
#'  \item{pcd2}{The computed Pseudo-count \eqn{D^{2}} statistic for each item, which quantifies
#'   the degree of item parameter drift by comparing observed and expected response frequencies.}
##' \item{N}{The number of examinees whose responses were used to calculate the Pseudo-count
#'   \eqn{D^{2}} statistic for each item, reflecting the sample size involved in the computation.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Cappaert, K. J., Wen, Y., & Chang, Y. F. (2018). Evaluating CAT-adjusted
#' approaches for suspected item parameter drift detection. \emph{Measurement:
#' Interdisciplinary Research and Perspectives, 16}(4), 226-238.
#'
#' Stone, C. A. (2000). Monte Carlo based null distribution for an alternative
#' goodness-of-fit test statistic in IRT models. \emph{Journal of educational
#' measurement, 37}(1), 58-75.
#'
#' @examples
#' ## Compute the pseudo-count D2 statistics for the dichotomous items
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # select the first 30 3PLM item metadata to be examined
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[1:30, 1:6]
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(25)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # simulate the response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' # compute the pseudo-count D2 statistics
#' ps_d2 <- pcd2(x = x, data = data)
#' print(ps_d2)
#'
#' @import dplyr
#' @export
pcd2 <- function(x, data, D = 1, missing = NA, Quadrature = c(49, 6.0), weights = NULL,
                 group.mean = 0.0, group.var = 1.0) {

  # transform a data set to matrix
  data <- data.matrix(data)

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # stop when the model includes any polytomous model
  if (any(x$model %in% c("GRM", "GPCM")) | any(x$cats > 2)) {
    stop("The current version only supports dichotomous response data.", call. = FALSE)
  }

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

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # count the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call. = FALSE)
  }

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
                      sum((y / z) * ((x[, 1] / y) - h[, 1])^2)
                    })

  # return the results
  rst <-  data.frame(id = id, pcd2 = pc_d2, N = n.resp)
  rst


}


