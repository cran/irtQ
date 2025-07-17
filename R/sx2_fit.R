#' S-X2 Fit Statistic
#'
#' @description Computes the \eqn{S\text{-}X^2} item fit statistic proposed by
#'   Orlando and Thissen (2000, 2003). This statistic evaluates the fit of IRT
#'   models by comparing observed and expected item response frequencies across
#'   summed score groups.
#'
#' @inheritParams est_score
#' @param alpha A numeric value specifying the significance level (\eqn{\alpha})
#'   for the hypothesis test associated with the \eqn{S\text{-}X^2} statistic.
#'   Default is 0.05.
#' @param min.collapse An integer specifying the minimum expected frequency
#'   required per cell before adjacent cells are collapsed. Default is 1. See
#'   **Details**.
#' @param nquad An integer specifying the number of Gaussian quadrature points
#'   used to approximate the normal prior distribution. Default is 30.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (first column) and their corresponding weights (second column) for
#'   the latent ability distribution. If omitted, default values are generated
#'   using [irtQ::gen.weight()] according to the `norm.prior` and `nquad` arguments.
#' @param pcm.loc An optional integer vector indicating the row indices of items
#'   that follow the partial credit model (PCM), where slope parameters are
#'   fixed. Default is `NULL`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @details
#' The accuracy of the \eqn{\chi^{2}} approximation in item fit statistics can be
#' compromised when expected cell frequencies in contingency tables are too small
#' (Orlando & Thissen, 2000). To address this issue, Orlando and Thissen (2000)
#' proposed collapsing adjacent summed score groups to ensure a minimum expected
#' frequency of at least 1.
#'
#' However, applying this collapsing approach directly to polytomous item data
#' can result in excessive information loss (Kang & Chen, 2008). To mitigate this,
#' Kang and Chen (2008) instead collapsed adjacent response categories *within*
#' each summed score group, maintaining a minimum expected frequency of 1 per
#' category. The same collapsing strategies are implemented in [irtQ::sx2_fit()].
#' If a different minimum expected frequency is desired, it can be specified via
#' the `min.collapse` argument.
#'
#' When an item is labeled as "DRM" in the item metadata, it is treated as a 3PLM
#' item when computing the degrees of freedom for the \eqn{S\text{-}X^2} statistic.
#'
#' Additionally, any missing responses in the `data` are automatically replaced
#' with incorrect responses (i.e., 0s).
#'
#' @return
#' A list containing the following components:
#' \item{fit_stat}{A data frame summarizing the \eqn{S\text{-}X^2} fit statistics
#' for all items, including the chi-square value, degrees of freedom, critical
#' value, and p-value.}
#' \item{item_df}{A data frame containing the item metadata as specified in the
#' input argument `x`.}
#' \item{exp_freq}{A list of collapsed expected frequency tables for all items.}
#' \item{obs_freq}{A list of collapsed observed frequency tables for all items.}
#' \item{exp_prob}{A list of collapsed expected probability tables for all items.}
#' \item{obs_prop}{A list of collapsed observed proportion tables for all items.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::irtfit()], [irtQ::simdat()], [irtQ::shape_df()],
#' [irtQ::est_irt()], [irtQ::est_item()]
#'
#' @references Kang, T., & Chen, T. T. (2008). Performance of the generalized
#' S-X2 item fit index for polytomous IRT models.
#' *Journal of Educational Measurement, 45*(4), 391-406.
#'
#' Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for
#' dichotomous item response theory models.
#' *Applied Psychological Measurement, 24*(1), 50-64.
#'
#' Orlando, M., & Thissen, D. (2003). Further investigation of the performance
#' of S-X2: An item fit index for use with dichotomous item response theory
#' models. *Applied Psychological Measurement, 27*(4), 289-298.
#'
#' @examples
#' ## Example 1: All five polytomous IRT items follow the GRM
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select the item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # Compute fit statistics
#' fit1 <- sx2_fit(x = x, data = data, nquad = 30)
#'
#' # Display fit statistics
#' fit1$fit_stat
#' }
#'
#' ## Example 2: Items 39 and 40 follow the GRM, and items 53, 54, and 55
#' ##            follow the PCM (with slope parameters fixed to 1)
#' # Replace the model names with "GPCM" and
#' # set the slope parameters of items 53–55 to 1
#' x[53:55, 3] <- "GPCM"
#' x[53:55, 4] <- 1
#'
#' # Generate examinees' abilities from N(0, 1)
#' set.seed(25)
#' score <- rnorm(1000, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # Compute fit statistics
#' fit2 <- sx2_fit(x = x, data = data, nquad = 30, pcm.loc = 53:55)
#'
#' # Display fit statistics
#' fit2$fit_stat
#' }
#'
#' @export
sx2_fit <- function(x, ...) UseMethod("sx2_fit")


#' @describeIn sx2_fit Default method for computing \eqn{S\text{-}X^{2}} fit
#' statistics from a data frame `x` containing item metadata.
#' @importFrom Rfast rowsums
#' @import dplyr
#' @export
sx2_fit.default <- function(x,
                            data,
                            D = 1,
                            alpha = 0.05,
                            min.collapse = 1,
                            norm.prior = c(0, 1),
                            nquad = 30,
                            weights,
                            pcm.loc = NULL,
                            ...) {

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

  # transform a data set to matrix
  data <- data.matrix(data)

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


#' @describeIn sx2_fit An object created by the function [irtQ::est_item()].
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

  # transform a data set to matrix
  data <- data.matrix(data)

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

#' @describeIn sx2_fit An object created by the function [irtQ::est_irt()].
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

  # transform a data set to matrix
  data <- data.matrix(data)

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
