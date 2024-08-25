#' S-X2 fit statistic
#'
#' @description This function computes \eqn{S-X^{2}} (Orlando & Thissen, 2000, 2003) item fit statistic.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...), an object
#' of class \code{\link{est_item}} obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}}
#' obtained from the function \code{\link{est_irt}}. See \code{\link{irtfit}}, \code{\link{info}}, or \code{\link{simdat}}
#' for more details about the item metadata. The data frame of item metadata can be easily obtained using the function \code{\link{shape_df}}.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param alpha A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test for \eqn{S-X^{2}} fit statistic. Default is .05.
#' @param min.collapse An integer value to indicate the minimum frequency of cells to be collapsed. Default is 1. See below for details.
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
#' c(0,1).
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 30.
#' @param weights A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
#' using the function \code{\link{gen.weight}}. If missing, default values are used (see the arguments of \code{norm.prior} and \code{nquad}).
#' @param pcm.loc A vector of integer values indicating the locations of partial credit model (PCM) items whose slope parameters are fixed
#' to certain values. Default is NULL.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Often, very small expected frequencies in the contingency tables used to compute \eqn{\chi^{2}} fit statistics could
#' compromise the accuracy of the \eqn{\chi^{2}} approximation for their distribution (Orlando & Thissen, 2000).
#' To avoid this problem, Orlando and Thissen (2000) used an algorithm of collapsing adjacent test score groups to maintain
#' a minimum expected category frequency of 1. However, if Orlando and Thissen's cell collapsing approach is applied to polytomous data,
#' too much information would be lost (Kang & Chen, 2008). Thus, Kang and Chen (2008) collapsed adjacent cells of item score categories
#' for a specific score group to ensure a minimum expected category frequency of 1. The same collapsing strategies were applied
#' in the function \code{\link{sx2_fit}}. If a minimum expected category frequency needs to be set to different number, you can specify
#' the minimum value in the argument \code{min.collapse}.
#'
#' Note that if "DRM" is specified for an item in the item metadata set, the item is considered as "3PLM" to compute degree of freedom of
#' the \eqn{S-X^{2}} fit statistic.
#'
#' Also, any missing responses in \code{data} are replaced with incorrect responses (i.e., 0s).
#'
#' @return This function returns a list. Within a list, several internal objects are contained such as:
#' \item{fit_stat}{A data frame containing the results of \eqn{S-X^{2}} fit statistics for all items.}
#' \item{item_df}{The item metadata specified in the argument \code{x}.}
#' \item{exp_freq}{A list containing the collapsed expected frequency tables for all items.}
#' \item{obs_freq}{A list containing the collapsed observed frequency tables for all items.}
#' \item{exp_prob}{A list containing the collapsed expected probability tables for all items.}
#' \item{obs_prop}{A list containing the collapsed observed proportion tables for all items.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{irtfit}}, \code{\link{info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{est_item}}
#'
#' @references
#' Kang, T., & Chen, T. T. (2008). Performance of the generalized S-X2 item fit index for polytomous IRT models.
#' \emph{Journal of Educational Measurement, 45}(4), 391-406.
#'
#' Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for dichotomous item response theory models.
#' \emph{Applied Psychological Measurement, 24}(1), 50-64.
#'
#' Orlando, M., & Thissen, D. (2003). Further investigation of the performance of S-X2: An item fit index for use with
#' dichotomous item response theory models. \emph{Applied Psychological Measurement, 27}(4), 289-298.
#'
#' @examples
#'
#' ## example 1: all five polytomous IRT models are GRM
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # select the item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # simulate the response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # compute fit statistics
#' fit1 <- sx2_fit(x = x, data = data, nquad = 30)
#'
#' # fit statistics
#' fit1$fit_stat
#' }
#'
#' ## example 2: first 39th and 40th items follows GRM and 53rd, 54th, and 55th items
#' ##            follow PCM (thus, the slope parameters are fixed to 1)
#' # replace the model names with GPCM and
#' # assign 1 to the slope parameters for the 53rd, 54th, and 55th items
#' x[53:55, 3] <- "GPCM"
#' x[53:55, 4] <- 1
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(25)
#' score <- rnorm(1000, mean = 0, sd = 1)
#'
#' # simulate the response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # compute fit statistics
#' fit2 <- sx2_fit(x = x, data = data, nquad = 30, pcm.loc = 53:55)
#'
#' # fit statistics
#' fit2$fit_stat
#' }
#'
#' @export
sx2_fit <- function(x, ...) UseMethod("sx2_fit")


#' @describeIn sx2_fit Default method to compute \eqn{S-X^{2}} fit statistics for a data frame \code{x} containing the item metadata.
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
sx2_fit.default <- function(x, data, D = 1, alpha = 0.05, min.collapse = 1, norm.prior = c(0, 1),
                            nquad = 30, weights, pcm.loc = NULL, ...) {
  ## ------------------------------------------------------------------------------------------------
  # check missing data
  # replace NAs with 0
  na.lg <- is.na(data)
  if (any(na.lg)) {
    data[na.lg] <- 0
    memo <- "Any missing responses are replaced with 0s. \n"
    warning(memo, call. = FALSE)
  }

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  ## ------------------------------------------------------------------
  ## 1. data preparation
  ## ------------------------------------------------------------------
  # save item metadata to another objective
  full_df <- x

  # count the number of items
  nitem <- nrow(x)

  # compute raw sum scores for all examinees
  rawscore <- Rfast::rowsums(data)

  # generate weights and thetas
  if (missing(weights)) {
    wts <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    wts <- data.frame(weights)
    nquad <- nrow(wts)
  }

  # select category variable
  cats <- elm_item$cats

  # total sum score
  t.score <- sum(cats - 1)

  # frequencies of raw sum scores
  score.freq <- as.numeric(table(factor(rawscore, levels = c(0:t.score))))

  ## ------------------------------------------------------------------
  ## 2. calculate the likelihoods
  ## ------------------------------------------------------------------
  # compute category probabilities for all items
  prob.cats <- trace(elm_item = elm_item, theta = wts[, 1], D = D, tcc = FALSE)$prob.cats

  # estimate likelihoods of getting raw sum scores using lord-wingersky algorithm
  lkhd <- t(lwRecurive(prob.cats = prob.cats, cats = cats, n.theta = nquad))

  # estimate likelihoods of getting raw sum scores except the examined item using lord-wingersky algorithm
  lkhd_noitem <-
    purrr::map(
      .x = 1:nitem,
      .f = function(i) {
        t(lwRecurive(prob.cats = prob.cats[-i], cats = cats[-i], n.theta = nquad))
      }
    )

  ## ------------------------------------------------------------------
  ## 3. prepare the contingency tables
  ## ------------------------------------------------------------------
  # contingency tables of the expected frequencies for all items
  exp_freq <-
    purrr::pmap(
      .l = list(x = cats, y = prob.cats, z = lkhd_noitem),
      .f = function(x, y, z) {
        expFreq(t.score,
          cats = x, prob.cats = y,
          lkhd_noitem = z, lkhd, wts, score.freq
        )
      }
    )

  # contingency tables of the observed frequencies for all items
  obs_freq <-
    purrr::map2(
      .x = data.frame(data), .y = cats,
      .f = function(x, y) {
        obsFreq(
          rawscore = rawscore, response = x,
          t.score = t.score, cats = y
        )
      }
    )

  ## ------------------------------------------------------------------
  ## 4. collapse cells in the contingency tables
  ## ------------------------------------------------------------------
  # cbind two lists of the expected frequency tables and the observed frequency tables
  ftable_info <- cbind(exp_freq, obs_freq)

  # collapse the expected and observed frequency tables
  # (1) DRM items
  if (!is.null(idx.drm)) {
    # select frequency tables of dichotomous items
    if (length(idx.drm) > 1) {
      ftable_info_drm <- ftable_info[idx.drm, ]
    } else {
      ftable_info_drm <- rbind(ftable_info[idx.drm, ])
    }

    # collapse the frequency tables for all dichotomous items
    for (i in 1:length(idx.drm)) {
      x <- data.frame(ftable_info_drm[i, ])

      # collapse based on the column of the incorrect responses
      tmp1 <- collapse_ftable(x = x, col = 1, min.collapse = min.collapse)

      # collapse based on the column of the correct responses
      tmp2 <- collapse_ftable(x = tmp1, col = 2, min.collapse = min.collapse)

      # replace the frequency tables with the collapsed frequency tables
      ftable_info_drm[i, ][[1]] <- tmp2[, 1:2]
      ftable_info_drm[i, ][[2]] <- tmp2[, 3:4]
    }

    ftable_info[idx.drm, ] <- ftable_info_drm
  }

  # (2) PRM items
  if (!is.null(idx.prm)) {
    # select frequency tables of polytomous items
    if (length(idx.prm) > 1) {
      ftable_info_plm <- ftable_info[idx.prm, ]
    } else {
      ftable_info_plm <- rbind(ftable_info[idx.prm, ])
    }

    # collapse the frequency tables for all polytomous items
    for (i in 1:length(idx.prm)) {
      # select the expected and observed frequency tables for the corresponding items
      exp_tmp <- data.frame(ftable_info_plm[i, 1])
      obs_tmp <- data.frame(ftable_info_plm[i, 2])

      # check if there are any rows a sum of the number of examinees across all score categories are zero
      # if so, delete them
      if (any(rowSums(exp_tmp) == 0L)) {
        exp_tmp <- exp_tmp[rowSums(exp_tmp) != 0L, ]
        obs_tmp <- obs_tmp[rowSums(obs_tmp) != 0L, ]
      }
      col.name <- colnames(exp_tmp)

      # create two empty list to contain the collapsed results
      exp_table <- vector("list", nrow(exp_tmp))
      obs_table <- vector("list", nrow(obs_tmp))

      # collapsing cells of each frequency table
      for (j in 1:nrow(exp_tmp)) {
        x <- data.frame(exp = as.numeric(exp_tmp[j, ]), obs = as.numeric(obs_tmp[j, ]))
        tmp <- collapse_ftable(x = x, col = 1, min.collapse = min.collapse)
        exp_table[[j]] <- tmp$exp
        obs_table[[j]] <- tmp$obs
      }

      ftable_info_plm[i, ][[1]] <-
        data.frame(bind.fill(exp_table, type = "rbind")) %>%
        stats::setNames(nm = col.name)
      ftable_info_plm[i, ][[2]] <- data.frame(bind.fill(obs_table, type = "rbind")) %>%
        stats::setNames(nm = col.name)
    }

    ftable_info[idx.prm, ] <- ftable_info_plm
  }

  # replace the two frequency tables with the collapsed two frequency tables
  exp_freq2 <-
    ftable_info[, 1] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "exp_freq.", replacement = "")
    })

  obs_freq2 <-
    ftable_info[, 2] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "obs_freq.", replacement = "")
    })

  ## ------------------------------------------------------------------
  ## 5. calculate the fit statistics
  ## ------------------------------------------------------------------
  # check the number of parameters for each item
  model <- full_df$model
  model[pcm.loc] <- "PCM"
  count_prm <- rep(NA, length(model))
  count_prm[model %in% "1PLM"] <- 1
  count_prm[model %in% "2PLM"] <- 2
  count_prm[model %in% c("3PLM", "DRM")] <- 3
  count_prm[model %in% "PCM"] <- full_df[model %in% "PCM", 2] - 1
  count_prm[model %in% "GPCM"] <- full_df[model %in% "GPCM", 2]
  count_prm[model %in% "GRM"] <- full_df[model %in% "GRM", 2]

  # compute the fit statistics for all items
  infoList <- list(exp_freq2, obs_freq2, as.list(count_prm))
  fitstat_list <- purrr::pmap(.l = infoList, .f = chisq_stat, crt.delta = 0.0, alpha = alpha)

  # make a data.frame for the fit statistics results
  fit_stat <- list(chisq = NULL, df = NULL, crit.val = NULL)
  for (i in 1:3) {
    fit_stat[[i]] <- purrr::map_dbl(fitstat_list, .f = function(x) x[[i]])
  }
  pval <- purrr::pmap_dbl(
    .l = list(x = fit_stat[[1]], y = fit_stat[[2]]),
    .f = function(x, y) 1 - stats::pchisq(q = x, df = y, lower.tail = TRUE)
  )
  fit_stat <- data.frame(id = full_df$id, fit_stat, p = round(pval, 3), stringsAsFactors = FALSE)
  fit_stat$chisq <- round(fit_stat$chisq, 3)
  fit_stat$crit.val <- round(fit_stat$crit.val, 3)
  rownames(fit_stat) <- NULL

  # extract the expected and observed proportion tables
  exp_prob <- purrr::map(fitstat_list, .f = function(x) x$exp.prop)
  obs_prop <- purrr::map(fitstat_list, .f = function(x) x$obs.prop)

  ## ------------------------------------------------------------------
  ## 3. return the results
  ## ------------------------------------------------------------------
  list(
    fit_stat = fit_stat, item_df = full_df, exp_freq = exp_freq2, obs_freq = obs_freq2,
    exp_prob = exp_prob, obs_prop = obs_prop
  )
}


#' @describeIn sx2_fit An object created by the function \code{\link{est_item}}.
#' @import dplyr
#' @export
sx2_fit.est_item <- function(x, alpha = 0.05, min.collapse = 1, norm.prior = c(0, 1),
                             nquad = 30, weights, pcm.loc = NULL, ...) {
  # extract information from an object
  data <- x$data
  D <- x$scale.D
  x <- x$par.est

  ## ------------------------------------------------------------------------------------------------
  # check missing data
  # replace NAs with 0
  na.lg <- is.na(data)
  if (any(na.lg)) {
    data[na.lg] <- 0
    memo <- "Any missing responses are replaced with 0s. \n"
    warning(memo, call. = FALSE)
  }

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  ## ------------------------------------------------------------------
  ## 1. data preparation
  ## ------------------------------------------------------------------
  # save item metadata to another objective
  full_df <- x

  # count the number of items
  nitem <- nrow(x)

  # compute raw sum scores for all examinees
  rawscore <- Rfast::rowsums(data)

  # generate weights and thetas
  if (missing(weights)) {
    wts <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    wts <- data.frame(weights)
    nquad <- nrow(wts)
  }

  # select category variable
  cats <- elm_item$cats

  # total sum score
  t.score <- sum(cats - 1)

  # frequencies of raw sum scores
  score.freq <- as.numeric(table(factor(rawscore, levels = c(0:t.score))))

  ## ------------------------------------------------------------------
  ## 2. calculate the likelihoods
  ## ------------------------------------------------------------------
  # compute category probabilities for all items
  prob.cats <- trace(elm_item = elm_item, theta = wts[, 1], D = D, tcc = FALSE)$prob.cats

  # estimate likelihoods of getting raw sum scores using lord-wingersky algorithm
  lkhd <- t(lwRecurive(prob.cats = prob.cats, cats = cats, n.theta = nquad))

  # estimate likelihoods of getting raw sum scores except the examined item using lord-wingersky algorithm
  lkhd_noitem <-
    purrr::map(
      .x = 1:nitem,
      .f = function(i) {
        t(lwRecurive(prob.cats = prob.cats[-i], cats = cats[-i], n.theta = nquad))
      }
    )

  ## ------------------------------------------------------------------
  ## 3. prepare the contingency tables
  ## ------------------------------------------------------------------
  # contingency tables of the expected frequencies for all items
  exp_freq <-
    purrr::pmap(
      .l = list(x = cats, y = prob.cats, z = lkhd_noitem),
      .f = function(x, y, z) {
        expFreq(t.score,
          cats = x, prob.cats = y,
          lkhd_noitem = z, lkhd, wts, score.freq
        )
      }
    )

  # contingency tables of the observed frequencies for all items
  obs_freq <-
    purrr::map2(
      .x = data.frame(data), .y = cats,
      .f = function(x, y) {
        obsFreq(
          rawscore = rawscore, response = x,
          t.score = t.score, cats = y
        )
      }
    )

  ## ------------------------------------------------------------------
  ## 4. collapse cells in the contingency tables
  ## ------------------------------------------------------------------
  # cbind two lists of the expected frequency tables and the observed frequency tables
  ftable_info <- cbind(exp_freq, obs_freq)

  # collapse the expected and observed frequency tables
  # (1) DRM items
  if (!is.null(idx.drm)) {
    # select frequency tables of dichotomous items
    if (length(idx.drm) > 1) {
      ftable_info_drm <- ftable_info[idx.drm, ]
    } else {
      ftable_info_drm <- rbind(ftable_info[idx.drm, ])
    }

    # collapse the frequency tables for all dichotomous items
    for (i in 1:length(idx.drm)) {
      x <- data.frame(ftable_info_drm[i, ])

      # collapse based on the column of the incorrect responses
      tmp1 <- collapse_ftable(x = x, col = 1, min.collapse = min.collapse)

      # collapse based on the column of the correct responses
      tmp2 <- collapse_ftable(x = tmp1, col = 2, min.collapse = min.collapse)

      # replace the frequency tables with the collapsed frequency tables
      ftable_info_drm[i, ][[1]] <- tmp2[, 1:2]
      ftable_info_drm[i, ][[2]] <- tmp2[, 3:4]
    }

    ftable_info[idx.drm, ] <- ftable_info_drm
  }

  # (2) PRM items
  if (!is.null(idx.prm)) {
    # select frequency tables of polytomous items
    if (length(idx.prm) > 1) {
      ftable_info_plm <- ftable_info[idx.prm, ]
    } else {
      ftable_info_plm <- rbind(ftable_info[idx.prm, ])
    }

    # collapse the frequency tables for all polytomous items
    for (i in 1:length(idx.prm)) {
      # select the expected and observed frequency tables for the corresponding items
      exp_tmp <- data.frame(ftable_info_plm[i, 1])
      obs_tmp <- data.frame(ftable_info_plm[i, 2])

      # check if there are any rows a sum of the number of examinees across all score categories are zero
      # if so, delete them
      if (any(rowSums(exp_tmp) == 0L)) {
        exp_tmp <- exp_tmp[rowSums(exp_tmp) != 0L, ]
        obs_tmp <- obs_tmp[rowSums(obs_tmp) != 0L, ]
      }
      col.name <- colnames(exp_tmp)

      # create two empty list to contain the collapsed results
      exp_table <- vector("list", nrow(exp_tmp))
      obs_table <- vector("list", nrow(obs_tmp))

      # collapsing cells of each frequency table
      for (j in 1:nrow(exp_tmp)) {
        x <- data.frame(exp = as.numeric(exp_tmp[j, ]), obs = as.numeric(obs_tmp[j, ]))
        tmp <- collapse_ftable(x = x, col = 1, min.collapse = min.collapse)
        exp_table[[j]] <- tmp$exp
        obs_table[[j]] <- tmp$obs
      }

      ftable_info_plm[i, ][[1]] <-
        data.frame(bind.fill(exp_table, type = "rbind")) %>%
        stats::setNames(nm = col.name)
      ftable_info_plm[i, ][[2]] <- data.frame(bind.fill(obs_table, type = "rbind")) %>%
        stats::setNames(nm = col.name)
    }

    ftable_info[idx.prm, ] <- ftable_info_plm
  }

  # replace the two frequency tables with the collapsed two frequency tables
  exp_freq2 <-
    ftable_info[, 1] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "exp_freq.", replacement = "")
    })

  obs_freq2 <-
    ftable_info[, 2] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "obs_freq.", replacement = "")
    })

  ## ------------------------------------------------------------------
  ## 5. calculate the fit statistics
  ## ------------------------------------------------------------------
  # check the number of parameters for each item
  model <- full_df$model
  model[pcm.loc] <- "PCM"
  count_prm <- rep(NA, length(model))
  count_prm[model %in% "1PLM"] <- 1
  count_prm[model %in% "2PLM"] <- 2
  count_prm[model %in% c("3PLM", "DRM")] <- 3
  count_prm[model %in% "PCM"] <- full_df[model %in% "PCM", 2] - 1
  count_prm[model %in% "GPCM"] <- full_df[model %in% "GPCM", 2]
  count_prm[model %in% "GRM"] <- full_df[model %in% "GRM", 2]

  # compute the fit statistics for all items
  infoList <- list(exp_freq2, obs_freq2, as.list(count_prm))
  fitstat_list <- purrr::pmap(.l = infoList, .f = chisq_stat, crt.delta = 0.0, alpha = alpha)

  # make a data.frame for the fit statistics results
  fit_stat <- list(chisq = NULL, df = NULL, crit.val = NULL)
  for (i in 1:3) {
    fit_stat[[i]] <- purrr::map_dbl(fitstat_list, .f = function(x) x[[i]])
  }
  pval <- purrr::pmap_dbl(
    .l = list(x = fit_stat[[1]], y = fit_stat[[2]]),
    .f = function(x, y) 1 - stats::pchisq(q = x, df = y, lower.tail = TRUE)
  )
  fit_stat <- data.frame(id = full_df$id, fit_stat, p = round(pval, 3), stringsAsFactors = FALSE)
  fit_stat$chisq <- round(fit_stat$chisq, 3)
  fit_stat$crit.val <- round(fit_stat$crit.val, 3)
  rownames(fit_stat) <- NULL

  # extract the expected and observed proportion tables
  exp_prob <- purrr::map(fitstat_list, .f = function(x) x$exp.prop)
  obs_prop <- purrr::map(fitstat_list, .f = function(x) x$obs.prop)

  ## ------------------------------------------------------------------
  ## 3. return the results
  ## ------------------------------------------------------------------
  list(
    fit_stat = fit_stat, item_df = full_df, exp_freq = exp_freq2, obs_freq = obs_freq2,
    exp_prob = exp_prob, obs_prop = obs_prop
  )
}

#' @describeIn sx2_fit An object created by the function \code{\link{est_irt}}.
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
sx2_fit.est_irt <- function(x, alpha = 0.05, min.collapse = 1, norm.prior = c(0, 1),
                            nquad = 30, weights, pcm.loc = NULL, ...) {
  # extract information from an object
  data <- x$data
  D <- x$scale.D
  x <- x$par.est

  ## ------------------------------------------------------------------------------------------------
  # check missing data
  # replace NAs with 0
  na.lg <- is.na(data)
  if (any(na.lg)) {
    data[na.lg] <- 0
    memo <- "Any missing responses are replaced with 0s. \n"
    warning(memo, call. = FALSE)
  }

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  ## ------------------------------------------------------------------
  ## 1. data preparation
  ## ------------------------------------------------------------------
  # save item metadata to another objective
  full_df <- x

  # count the number of items
  nitem <- nrow(x)

  # compute raw sum scores for all examinees
  rawscore <- Rfast::rowsums(data)

  # generate weights and thetas
  if (missing(weights)) {
    wts <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    wts <- data.frame(weights)
    nquad <- nrow(wts)
  }

  # select category variable
  cats <- elm_item$cats

  # total sum score
  t.score <- sum(cats - 1)

  # frequencies of raw sum scores
  score.freq <- as.numeric(table(factor(rawscore, levels = c(0:t.score))))

  ## ------------------------------------------------------------------
  ## 2. calculate the likelihoods
  ## ------------------------------------------------------------------
  # compute category probabilities for all items
  prob.cats <- trace(elm_item = elm_item, theta = wts[, 1], D = D, tcc = FALSE)$prob.cats

  # estimate likelihoods of getting raw sum scores using lord-wingersky algorithm
  lkhd <- t(lwRecurive(prob.cats = prob.cats, cats = cats, n.theta = nquad))

  # estimate likelihoods of getting raw sum scores except the examined item using lord-wingersky algorithm
  lkhd_noitem <-
    purrr::map(
      .x = 1:nitem,
      .f = function(i) {
        t(lwRecurive(prob.cats = prob.cats[-i], cats = cats[-i], n.theta = nquad))
      }
    )

  ## ------------------------------------------------------------------
  ## 3. prepare the contingency tables
  ## ------------------------------------------------------------------
  # contingency tables of the expected frequencies for all items
  exp_freq <-
    purrr::pmap(
      .l = list(x = cats, y = prob.cats, z = lkhd_noitem),
      .f = function(x, y, z) {
        expFreq(t.score,
          cats = x, prob.cats = y,
          lkhd_noitem = z, lkhd, wts, score.freq
        )
      }
    )

  # contingency tables of the observed frequencies for all items
  obs_freq <-
    purrr::map2(
      .x = data.frame(data), .y = cats,
      .f = function(x, y) {
        obsFreq(
          rawscore = rawscore, response = x,
          t.score = t.score, cats = y
        )
      }
    )

  ## ------------------------------------------------------------------
  ## 4. collapse cells in the contingency tables
  ## ------------------------------------------------------------------
  # cbind two lists of the expected frequency tables and the observed frequency tables
  ftable_info <- cbind(exp_freq, obs_freq)

  # collapse the expected and observed frequency tables
  # (1) DRM items
  if (!is.null(idx.drm)) {
    # select frequency tables of dichotomous items
    if (length(idx.drm) > 1) {
      ftable_info_drm <- ftable_info[idx.drm, ]
    } else {
      ftable_info_drm <- rbind(ftable_info[idx.drm, ])
    }

    # collapse the frequency tables for all dichotomous items
    for (i in 1:length(idx.drm)) {
      x <- data.frame(ftable_info_drm[i, ])

      # collapse based on the column of the incorrect responses
      tmp1 <- collapse_ftable(x = x, col = 1, min.collapse = min.collapse)

      # collapse based on the column of the correct responses
      tmp2 <- collapse_ftable(x = tmp1, col = 2, min.collapse = min.collapse)

      # replace the frequency tables with the collapsed frequency tables
      ftable_info_drm[i, ][[1]] <- tmp2[, 1:2]
      ftable_info_drm[i, ][[2]] <- tmp2[, 3:4]
    }

    ftable_info[idx.drm, ] <- ftable_info_drm
  }

  # (2) PRM items
  if (!is.null(idx.prm)) {
    # select frequency tables of polytomous items
    if (length(idx.prm) > 1) {
      ftable_info_plm <- ftable_info[idx.prm, ]
    } else {
      ftable_info_plm <- rbind(ftable_info[idx.prm, ])
    }

    # collapse the frequency tables for all polytomous items
    for (i in 1:length(idx.prm)) {
      # select the expected and observed frequency tables for the corresponding items
      exp_tmp <- data.frame(ftable_info_plm[i, 1])
      obs_tmp <- data.frame(ftable_info_plm[i, 2])

      # check if there are any rows a sum of the number of examinees across all score categories are zero
      # if so, delete them
      if (any(rowSums(exp_tmp) == 0L)) {
        exp_tmp <- exp_tmp[rowSums(exp_tmp) != 0L, ]
        obs_tmp <- obs_tmp[rowSums(obs_tmp) != 0L, ]
      }
      col.name <- colnames(exp_tmp)

      # create two empty list to contain the collapsed results
      exp_table <- vector("list", nrow(exp_tmp))
      obs_table <- vector("list", nrow(obs_tmp))

      # collapsing cells of each frequency table
      for (j in 1:nrow(exp_tmp)) {
        x <- data.frame(exp = as.numeric(exp_tmp[j, ]), obs = as.numeric(obs_tmp[j, ]))
        tmp <- collapse_ftable(x = x, col = 1, min.collapse = min.collapse)
        exp_table[[j]] <- tmp$exp
        obs_table[[j]] <- tmp$obs
      }

      ftable_info_plm[i, ][[1]] <-
        data.frame(bind.fill(exp_table, type = "rbind")) %>%
        stats::setNames(nm = col.name)
      ftable_info_plm[i, ][[2]] <- data.frame(bind.fill(obs_table, type = "rbind")) %>%
        stats::setNames(nm = col.name)
    }

    ftable_info[idx.prm, ] <- ftable_info_plm
  }

  # replace the two frequency tables with the collapsed two frequency tables
  exp_freq2 <-
    ftable_info[, 1] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "exp_freq.", replacement = "")
    })

  obs_freq2 <-
    ftable_info[, 2] %>%
    purrr::map(.f = function(x) {
      dplyr::rename_all(x, .funs = list("gsub"), pattern = "obs_freq.", replacement = "")
    })

  ## ------------------------------------------------------------------
  ## 5. calculate the fit statistics
  ## ------------------------------------------------------------------
  # check the number of parameters for each item
  model <- full_df$model
  model[pcm.loc] <- "PCM"
  count_prm <- rep(NA, length(model))
  count_prm[model %in% "1PLM"] <- 1
  count_prm[model %in% "2PLM"] <- 2
  count_prm[model %in% c("3PLM", "DRM")] <- 3
  count_prm[model %in% "PCM"] <- full_df[model %in% "PCM", 2] - 1
  count_prm[model %in% "GPCM"] <- full_df[model %in% "GPCM", 2]
  count_prm[model %in% "GRM"] <- full_df[model %in% "GRM", 2]

  # compute the fit statistics for all items
  infoList <- list(exp_freq2, obs_freq2, as.list(count_prm))
  fitstat_list <- purrr::pmap(.l = infoList, .f = chisq_stat, crt.delta = 0.0, alpha = alpha)

  # make a data.frame for the fit statistics results
  fit_stat <- list(chisq = NULL, df = NULL, crit.val = NULL)
  for (i in 1:3) {
    fit_stat[[i]] <- purrr::map_dbl(fitstat_list, .f = function(x) x[[i]])
  }
  pval <- purrr::pmap_dbl(
    .l = list(x = fit_stat[[1]], y = fit_stat[[2]]),
    .f = function(x, y) 1 - stats::pchisq(q = x, df = y, lower.tail = TRUE)
  )
  fit_stat <- data.frame(id = full_df$id, fit_stat, p = round(pval, 3), stringsAsFactors = FALSE)
  fit_stat$chisq <- round(fit_stat$chisq, 3)
  fit_stat$crit.val <- round(fit_stat$crit.val, 3)
  rownames(fit_stat) <- NULL

  # extract the expected and observed proportion tables
  exp_prob <- purrr::map(fitstat_list, .f = function(x) x$exp.prop)
  obs_prop <- purrr::map(fitstat_list, .f = function(x) x$obs.prop)

  ## ------------------------------------------------------------------
  ## 3. return the results
  ## ------------------------------------------------------------------
  list(
    fit_stat = fit_stat, item_df = full_df, exp_freq = exp_freq2, obs_freq = obs_freq2,
    exp_prob = exp_prob, obs_prop = obs_prop
  )
}
