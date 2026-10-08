#' Classical Test Theory Item Analysis
#'
#' Computes traditional classical test theory (CTT) item statistics -
#' difficulty, an item-total correlation (discrimination), and Cronbach's
#' alpha with the item removed - for dichotomous or polytomous item response
#' data, along with optional flagging of items that fall outside commonly
#' used quality thresholds.
#'
#' @param data A data frame or matrix of already-scored item responses, with
#'   examinees in rows and items in columns. Item scores must range from 0 to
#'   `cats[j] - 1` for each item j (0/1 for a dichotomous item; 0, 1, 2, ...
#'   for a polytomous/partial-credit item). An error is raised when a score
#'   is not a whole number between 0 and `cats[j] - 1`.
#' @param item.id A character vector of item identifiers, in the same order
#'   as the columns of `data`. If `NULL` (default), item IDs are generated
#'   automatically as `paste0("V", 1:ncol(data))`, following the convention
#'   used elsewhere in irtQ (e.g., [irtQ::shape_df()], [irtQ::est_irt()]).
#'   Note that this generated ID does not fall back to `colnames(data)`; pass
#'   `item.id` explicitly to label items using the column names of `data` or
#'   any other identifier scheme.
#' @param cats A numeric vector giving the number of score categories for
#'   each item (e.g., 2 for a dichotomous item), following the `cats`
#'   convention used elsewhere in irtQ (see, e.g., [irtQ::shape_df()]). If
#'   `NULL` (default), the number of categories for each item is inferred
#'   from the observed maximum score in `data` plus one, with a minimum of
#'   two (i.e., `max(max(data[, j]) + 1, 2)`); supply `cats` explicitly
#'   whenever the maximum possible score may not have been observed in the
#'   sample. Each value must be a whole number of at least 2.
#' @param correct Logical. Both the raw (uncorrected) item-total correlation -
#'   where an item (or, in [ctt_distr()], a response option/category) is
#'   correlated with the total score that includes its own contribution -
#'   and the corrected item-total correlation - excluding its own
#'   contribution - are always computed and reported as separate columns
#'   (see **Details**/**Value**). This argument only controls which of the
#'   two is used for flagging: if `TRUE`, flagging (in [ctt_item()]'s `flag`
#'   column, [ctt_alpha()]'s summary discrimination figure, and, in
#'   [ctt_distr()], the `crit.distractor` check) is based on the corrected
#'   value; if `FALSE` (default), on the raw value.
#' @param missing A value indicating missing responses in `data`, analogous
#'   to the `missing` argument in [irtQ::est_irt()] and [irtQ::score_resp()].
#'   Any cell equal to `missing` is recoded to `NA` before analysis. Default
#'   is `NA`. Examinees with any remaining missing item response are excluded
#'   listwise from all statistics computed by this function.
#' @param flag Logical. If `TRUE` (default), items are flagged when their
#'   difficulty or discrimination falls outside the thresholds given in
#'   `crit.p` and `crit.dis`.
#' @param crit.p A numeric vector of length two giving the lower and upper
#'   difficulty bounds used for flagging: difficulty below the first value is
#'   flagged as too difficult, and difficulty above the second value is
#'   flagged as too easy. Default is `c(0.10, 0.95)`.
#' @param crit.dis A single numeric value giving the minimum acceptable
#'   discrimination (item-total correlation); items strictly below this value
#'   are flagged as poorly discriminating. Default is `0.20`.
#'
#' @details
#' Difficulty for item j is defined generally as the mean observed score
#' divided by the item's maximum possible score,
#' `mean(data[, j], na.rm = TRUE) / (cats[j] - 1)`. For a dichotomous item
#' (`cats[j] = 2`), this reduces to the familiar proportion-correct difficulty
#' index. For a polytomous item, this expresses the average score as a
#' proportion of the maximum attainable score, so that difficulty remains
#' interpretable on the same 0-1 scale regardless of the number of score
#' categories.
#'
#' Discrimination for item j is the Pearson correlation between the item
#' score and the total score, which for a dichotomous item is mathematically
#' equivalent to the point-biserial correlation. Both the uncorrected (raw)
#' item-total correlation and the corrected (item-excluded) item-total
#' correlation are always computed and returned as separate columns; see,
#' e.g., Crocker and Algina (1986) for discussion of both conventions.
#' This mirrors how some software reports both side by side (e.g.,
#' `psych::alpha()` reports the raw item-total correlation as `raw.r` and
#' the corrected version as `r.drop`) rather than defaulting to one or the
#' other. The `correct` argument only selects which of the two feeds the
#' discrimination flagging criterion (see `crit.dis`).
#'
#' Alpha-with-item-removed for item j is Cronbach's alpha recomputed on the
#' remaining `ncol(data) - 1` items, using the same variance-based formula
#' used for an overall, test-level alpha (`k / (k - 1) * (1 - sum(item
#' variances) / total variance)`). A value of `NA` is returned wherever a
#' needed variance is zero (e.g., a constant item, or fewer than two items
#' remaining), since the relevant ratio is then undefined.
#'
#' @return A list with two elements:
#' \item{item}{A data frame with one row per item, containing the item
#'   label, number of score categories, difficulty, the raw (uncorrected)
#'   item-total correlation (`discrimination_raw`), the corrected
#'   (item-excluded) item-total correlation (`discrimination_corrected`),
#'   and alpha-with-item-removed, plus a `flag` column (a character
#'   description of which criteria were triggered, or `""` when none were;
#'   present only when `flag = TRUE`). Flagging is based on whichever of
#'   `discrimination_raw`/`discrimination_corrected` is selected by
#'   `correct`.}
#' \item{crit}{A list echoing the `crit.p` and `crit.dis` thresholds used, for
#'   reference in downstream reporting.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::score_resp()]
#'
#' @references
#' Crocker, L., & Algina, J. (1986). *Introduction to classical and modern
#' test theory*. Holt, Rinehart and Winston.
#'
#' @keywords internal
ctt_item <- function(data, item.id = NULL, cats = NULL, correct = FALSE,
                      missing = NA, flag = TRUE, crit.p = c(0.10, 0.95),
                      crit.dis = 0.20) {

  # coerce to a plain data frame so column-wise access (data[[j]]) behaves
  # consistently for matrix/tibble/data.frame input alike
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  # number of items (columns) in the scored response data; at least two are
  # required so that item-total correlations and alpha are well defined
  n_item <- ncol(data)
  if (n_item < 2L) {
    stop("`data` must contain at least two items to compute item-total ",
         "correlations and alpha.", call. = FALSE)
  }

  # item labels: use the supplied item.id if given (validating its length),
  # otherwise auto-generate V1, V2, ... following irtQ's own convention
  # (this does NOT fall back to colnames(data); see @param item.id)
  if (is.null(item.id)) {
    item_names <- paste0("V", seq_len(n_item))
  } else {
    if (length(item.id) != n_item) {
      stop("length(item.id) must equal ncol(data): one ID per item.",
           call. = FALSE)
    }
    item_names <- as.character(item.id)
  }

  # recode a user-specified missing-value sentinel to NA before analysis,
  # mirroring the `missing` argument convention used elsewhere in irtQ
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # listwise-delete any examinee with a remaining missing response, so that
  # every correlation and the total score are computed on a common sample
  complete_rows <- stats::complete.cases(data)
  n_dropped <- sum(!complete_rows)
  if (n_dropped > 0L) {
    warning(n_dropped, " examinee(s) with missing item responses were ",
            "excluded listwise from ctt_item().", call. = FALSE)
    data <- data[complete_rows, , drop = FALSE]
  }

  # stop when a score or a `cats` value is not valid
  check_ctt_scores(data, cats, item_names)

  # infer the number of score categories per item from the observed maximum
  # score, when not supplied explicitly (cats[j] = observed max score + 1,
  # with a minimum of two categories)
  if (is.null(cats)) {
    cats <- pmax(vapply(data, function(x) max(x) + 1, numeric(1)), 2)
  }

  # total score across all items (row sums), used for item-total correlations
  total <- rowSums(data)

  # pre-allocate per-item result vectors, filled in the loop below; both the
  # raw and corrected item-total correlations are always computed (see
  # @param correct), so both get their own vector regardless of `correct`
  difficulty <- numeric(n_item)
  discrimination_raw <- numeric(n_item)
  discrimination_corrected <- numeric(n_item)
  alpha_removed <- numeric(n_item)

  # loop over items; each item's statistics depend on the item itself and,
  # for discrimination/alpha, on the rest of the test
  for (j in seq_len(n_item)) {

    item_score <- data[[j]]     # this item's scored responses
    max_score <- cats[j] - 1    # maximum attainable score for this item

    # difficulty = mean score / max score; NA when the item has no possible
    # score range (max_score <= 0, i.e., a single-category/constant item)
    difficulty[j] <- if (max_score > 0) mean(item_score) / max_score else NA_real_

    # reference total scores for the raw (item-included) and corrected
    # (item-excluded) item-total correlations; both are always computed
    ref_total_raw <- total
    ref_total_corrected <- total - item_score

    # discrimination = Pearson correlation between item score and reference
    # total; NA when either vector is constant, since correlation is undefined
    sd_item_ok <- stats::sd(item_score) > 0
    discrimination_raw[j] <- if (sd_item_ok && stats::sd(ref_total_raw) > 0) {
      stats::cor(item_score, ref_total_raw)
    } else {
      NA_real_
    }
    discrimination_corrected[j] <- if (sd_item_ok && stats::sd(ref_total_corrected) > 0) {
      stats::cor(item_score, ref_total_corrected)
    } else {
      NA_real_
    }

    # alpha-with-item-removed: Cronbach's alpha recomputed using only the
    # remaining items (same variance-based formula as the test-level alpha)
    rest <- data[, -j, drop = FALSE]     # all items except item j
    k_rest <- ncol(rest)                  # number of remaining items
    if (k_rest >= 2L) {
      item_var_rest <- vapply(rest, stats::var, numeric(1))  # per-item variances
      total_var_rest <- stats::var(rowSums(rest))            # remaining total var
      alpha_removed[j] <- if (total_var_rest > 0) {
        (k_rest / (k_rest - 1)) * (1 - sum(item_var_rest) / total_var_rest)
      } else {
        NA_real_
      }
    } else {
      alpha_removed[j] <- NA_real_   # fewer than 2 remaining items: alpha undefined
    }
  }

  # discrimination values used for flagging: selected by `correct` (TRUE =
  # corrected/item-excluded, FALSE = raw/item-included), while both columns
  # are always reported in item_df below
  discrimination_for_flag <- if (correct) discrimination_corrected else discrimination_raw

  # assemble the per-item result data frame, rounding for readable reporting
  item_df <- data.frame(
    item = item_names,
    cats = unname(cats),
    difficulty = round(difficulty, 3),
    discrimination_raw = round(discrimination_raw, 3),
    discrimination_corrected = round(discrimination_corrected, 3),
    alpha_removed = round(alpha_removed, 3),
    stringsAsFactors = FALSE
  )

  # optional flagging based on the difficulty/discrimination thresholds;
  # each item's flag column lists every criterion it triggered, or "" if none
  if (flag) {
    flag_txt <- character(n_item)
    for (j in seq_len(n_item)) {
      msgs <- character(0)
      if (!is.na(difficulty[j]) && difficulty[j] < crit.p[1]) {
        msgs <- c(msgs, "difficulty too low")
      }
      if (!is.na(difficulty[j]) && difficulty[j] > crit.p[2]) {
        msgs <- c(msgs, "difficulty too high")
      }
      if (!is.na(discrimination_for_flag[j]) && discrimination_for_flag[j] < crit.dis) {
        msgs <- c(msgs, "discrimination too low")
      }
      flag_txt[j] <- paste(msgs, collapse = "; ")
    }
    item_df$flag <- flag_txt
  }

  # return the per-item statistics and the thresholds used to flag them
  list(item = item_df, crit = list(crit.p = crit.p, crit.dis = crit.dis))
}

# validates scored item responses and the `cats` vector for the CTT functions.
# `data` is a data frame with missing values already removed, `cats` is NULL
# or one number of categories per item, and `item_names` labels the items in
# the error messages. returns NULL invisibly and stops on the first problem
check_ctt_scores <- function(data, cats = NULL, item_names = NULL) {

  # number of items and default item labels for the messages
  n_item <- ncol(data)
  if (is.null(item_names)) {
    item_names <- paste0("V", seq_len(n_item))
  }

  # every column must hold numeric or logical scores
  not_numeric <- !vapply(data, function(x) is.numeric(x) || is.logical(x),
                         logical(1))
  if (any(not_numeric)) {
    stop("Item scores must be numeric. Non-numeric item(s): ",
         paste(item_names[not_numeric], collapse = ", "), ".", call. = FALSE)
  }

  # at least one examinee must remain after the missing-data handling
  if (nrow(data) == 0L) {
    stop("No examinee has complete item responses.", call. = FALSE)
  }

  # every score must be a finite whole number of at least 0
  not_whole <- vapply(data, function(x) {
    any(!is.finite(x) | x < 0 | abs(x - round(x)) > 1e-8)
  }, logical(1))
  if (any(not_whole)) {
    stop("Item scores must be whole numbers between 0 and cats - 1. ",
         "Invalid score(s) in item(s): ",
         paste(item_names[not_whole], collapse = ", "), ".", call. = FALSE)
  }

  # nothing more to check when `cats` is not supplied
  if (is.null(cats)) {
    return(invisible(NULL))
  }

  # one `cats` value per item, each a whole number of at least 2
  if (length(cats) != n_item) {
    stop("length(cats) must equal ncol(data): one value per item.",
         call. = FALSE)
  }
  if (!is.numeric(cats) || any(!is.finite(cats)) || any(cats < 2) ||
      any(abs(cats - round(cats)) > 1e-8)) {
    stop("`cats` must contain whole numbers of at least 2, one per item.",
         call. = FALSE)
  }

  # no score may exceed the maximum score implied by `cats`
  too_high <- vapply(seq_len(n_item), function(j) {
    max(data[[j]]) > cats[j] - 1 + 1e-8
  }, logical(1))
  if (any(too_high)) {
    stop("Item scores must be whole numbers between 0 and cats - 1. ",
         "Score(s) above cats - 1 in item(s): ",
         paste(item_names[too_high], collapse = ", "), ".", call. = FALSE)
  }

  invisible(NULL)
}
