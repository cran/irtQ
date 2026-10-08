#' Test-Level Reliability Summary (Cronbach's Alpha)
#'
#' Computes a test-level classical test theory (CTT) reliability summary -
#' Cronbach's alpha (both the raw and standardized forms), the standard error
#' of measurement (SEM), and the average item difficulty and discrimination -
#' from scored item response data.
#'
#' @inheritParams ctt_item
#'
#' @details
#' Two forms of Cronbach's alpha are always computed and reported as separate
#' columns:
#'
#' Raw alpha (`alpha`) uses the standard variance-based formula,
#' `alpha = (k / (k - 1)) * (1 - sum(item variances) / total-score
#' variance)`, where k is the number of items. This formula is a general
#' reliability coefficient that applies unchanged to dichotomous and
#' polytomous item scores alike, and it reflects the reliability of the
#' actual (unweighted) total score obtained by simply summing the item
#' scores - the score most tests actually use for reporting and decisions.
#' In everyday terms, raw alpha asks: "if every item kept its own natural
#' scale and spread, how consistently do these items agree with each other?"
#'
#' Standardized alpha (`alpha_std`) instead first standardizes every item to
#' the same scale (unit variance) before combining them, using the
#' equivalent formula `alpha_std = (k * r_bar) / (1 + (k - 1) * r_bar)`,
#' where `r_bar` is the average pairwise correlation among all items. In
#' everyday terms, standardized alpha asks: "if every item counted equally
#' regardless of how much it happens to vary in this particular sample, how
#' consistently would these items agree with each other?" Because it removes
#' the influence of any single item's variance, standardized alpha is most
#' useful when items differ substantially in scale or format (e.g., a mix of
#' dichotomous and polytomous items with very different score ranges); when
#' all items share the same scale and format (as with a dichotomous
#' selected-response test scored 0/1), raw and standardized alpha are
#' typically close, and raw alpha remains the more directly interpretable of
#' the two since it matches the reliability of the score actually used in
#' practice. See Cronbach (1951) for the original derivation of coefficient
#' alpha, and Osburn (2000) for a discussion contrasting the raw
#' (covariance-based) and standardized (correlation-based) forms.
#'
#' The standard error of measurement (SEM) is computed as
#' `SEM = SD(total score) * sqrt(1 - alpha)` (using raw alpha), following the
#' standard CTT relationship between test reliability and measurement
#' precision.
#'
#' Average difficulty and average discrimination are the simple means of the
#' per-item difficulty and discrimination values computed by [ctt_item()]
#' (any `NA` per-item value, e.g. from a constant item, is excluded via
#' `na.rm = TRUE`). As in [ctt_item()], both the raw (uncorrected) and
#' corrected (item-excluded) item-total correlations are always averaged and
#' reported as separate columns. `ctt_alpha()` itself has no flagging step,
#' so `correct` has no further effect here beyond being passed through to
#' [ctt_item()] for internal consistency.
#'
#' Because `mean_difficulty`/`mean_discrimination_raw`/
#' `mean_discrimination_corrected` are averaged with `na.rm = TRUE`, they
#' inherit [ctt_item()]'s `cats`-auto-inference rule: when `cats` is not
#' supplied, it is inferred as the observed maximum score plus one, with a
#' minimum of two. An item that every examinee scores 0 on therefore enters
#' the average with a difficulty of 0 and is flagged as "difficulty too low"
#' by [ctt_item()]. For a polytomous item whose highest score category is not
#' observed in the sample, however, `cats` is inferred too small and the
#' difficulty is biased upward, so supply `cats` explicitly whenever the
#' maximum possible score may not have been observed.
#'
#' @return A one-row data frame containing:
#' \item{n_examinee}{number of examinees included (after listwise deletion).}
#' \item{n_item}{number of items.}
#' \item{alpha}{Cronbach's alpha, raw (covariance-based) form.}
#' \item{alpha_std}{Cronbach's alpha, standardized (correlation-based) form -
#'   equivalent to computing raw alpha after first standardizing every item
#'   to unit variance. See **Details** for when this differs meaningfully
#'   from `alpha`.}
#' \item{sem}{the standard error of measurement (based on raw alpha).}
#' \item{mean_difficulty}{average item difficulty.}
#' \item{mean_discrimination_raw}{average raw (uncorrected) item-total
#'   correlation across items.}
#' \item{mean_discrimination_corrected}{average corrected (item-excluded)
#'   item-total correlation across items.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [ctt_item()], [irtQ::score_resp()]
#'
#' @references
#'   Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
#'   tests. *Psychometrika, 16*(3), 297-334. \doi{10.1007/BF02310555}.
#'
#'   Osburn, H. G. (2000). Coefficient alpha and related internal consistency
#'   reliability coefficients. *Psychological Methods, 5*(3), 343-355.
#'   \doi{10.1037/1082-989X.5.3.343}.
#'
#' @keywords internal
ctt_alpha <- function(data, item.id = NULL, cats = NULL, correct = FALSE,
                       missing = NA) {

  # coerce to a plain data frame so column-wise access (data[[j]]) behaves
  # consistently for matrix/tibble/data.frame input alike
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  # number of items (columns); at least two are required for alpha to be
  # mathematically defined
  n_item <- ncol(data)
  if (n_item < 2L) {
    stop("`data` must contain at least two items to compute alpha.",
         call. = FALSE)
  }

  # validate item.id length up front (also re-validated inside ctt_item(),
  # but checking here gives a clearer error before any other work is done)
  if (!is.null(item.id) && length(item.id) != n_item) {
    stop("length(item.id) must equal ncol(data): one ID per item.",
         call. = FALSE)
  }

  # recode a user-specified missing-value sentinel to NA before analysis,
  # mirroring the `missing` argument convention used elsewhere in irtQ
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # listwise-delete any examinee with a remaining missing response, so alpha
  # and the total-score variance are computed on a common sample
  complete_rows <- stats::complete.cases(data)
  n_dropped <- sum(!complete_rows)
  if (n_dropped > 0L) {
    warning(n_dropped, " examinee(s) with missing item responses were ",
            "excluded listwise from ctt_alpha().", call. = FALSE)
    data <- data[complete_rows, , drop = FALSE]
  }

  # stop when a score or a `cats` value is not valid
  check_ctt_scores(data, cats, item.id)
  n_examinee <- nrow(data)

  # per-item variances and the total-score variance, used in the alpha
  # formula below
  item_var <- vapply(data, stats::var, numeric(1))  # variance of each item
  total <- rowSums(data)                             # total score per examinee
  total_var <- stats::var(total)                     # variance of total score

  # Cronbach's alpha; NA when the total score has zero variance (e.g., every
  # examinee has the same total), since the ratio is then undefined
  alpha <- if (total_var > 0) {
    (n_item / (n_item - 1)) * (1 - sum(item_var) / total_var)
  } else {
    NA_real_
  }

  # standard error of measurement: SD(total score) * sqrt(1 - alpha)
  sem <- if (!is.na(alpha)) stats::sd(total) * sqrt(1 - alpha) else NA_real_

  # standardized alpha: equivalent to raw alpha computed after standardizing
  # every item to unit variance, via the average-inter-item-correlation form
  # alpha_std = (k * r_bar) / (1 + (k - 1) * r_bar); r_bar is the mean of all
  # off-diagonal (item, item) correlations. A constant item yields NA
  # correlations with every other item, so na.rm = TRUE excludes those pairs;
  # if every pairwise correlation is NA (e.g., all items constant), r_bar and
  # therefore alpha_std are NaN, mirroring alpha's own NA-when-undefined
  # behavior above.
  # the zero-standard-deviation warning for a constant item is suppressed
  # because its NA correlations are handled by na.rm = TRUE below
  item_cor <- suppressWarnings(stats::cor(data))     # k x k item correlation matrix
  r_bar <- mean(item_cor[upper.tri(item_cor)], na.rm = TRUE)  # mean off-diagonal r
  alpha_std <- (n_item * r_bar) / (1 + (n_item - 1) * r_bar)

  # average item difficulty/discrimination, obtained by calling ctt_item()
  # on the same (already missing-recoded and listwise-deleted) data, so the
  # two functions always agree on the underlying per-item formulas; missing
  # is left at its default (NA) here since recoding already happened above
  item_stats <- ctt_item(data = data, item.id = item.id, cats = cats,
                          correct = correct, missing = NA, flag = FALSE)$item
  mean_difficulty <- mean(item_stats$difficulty, na.rm = TRUE)
  mean_discrimination_raw <- mean(item_stats$discrimination_raw, na.rm = TRUE)
  mean_discrimination_corrected <- mean(item_stats$discrimination_corrected,
                                        na.rm = TRUE)

  # assemble the one-row test-level summary, rounded for readable reporting
  data.frame(
    n_examinee = n_examinee,
    n_item = n_item,
    alpha = round(alpha, 3),
    alpha_std = round(alpha_std, 3),
    sem = round(sem, 3),
    mean_difficulty = round(mean_difficulty, 3),
    mean_discrimination_raw = round(mean_discrimination_raw, 3),
    mean_discrimination_corrected = round(mean_discrimination_corrected, 3)
  )
}
