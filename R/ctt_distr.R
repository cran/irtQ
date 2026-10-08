#' Option/Category Response Distribution and Distractor Analysis
#'
#' Computes, for each item, the percentage of examinees choosing each
#' response option (or, for already-scored polytomous data, each score
#' category), together with each option's/category's point-biserial
#' correlation with the total score. When an answer key is supplied, the
#' function operates in "selected-response" mode - treating `data` as raw,
#' unscored option responses (as consumed by [irtQ::score_resp()]) - and
#' additionally reports omission/double-marking rates and flags distractors
#' (incorrect options) whose correlation with the total score is positive, a
#' common indicator of a miskeyed or poorly written item. When no key is
#' supplied, the function operates in "scored-category" mode - treating
#' `data` as an already-scored response matrix (0, 1, 2, ..., cats - 1) - and
#' simply reports each score category's response proportion and
#' category-total correlation, with no notion of a "correct" option to flag
#' distractors against.
#'
#' @note Selected-response mode (`key` supplied) is restricted to
#'   **dichotomous, single-key items only**, exactly like
#'   [irtQ::score_resp()], which it calls internally to obtain the 0/1 scored
#'   responses and total score: a selected-response item has exactly one
#'   correct option and therefore cannot yield a partial-credit score by
#'   construction, so a mixed dichotomous/polytomous test form cannot be
#'   analyzed in this mode. Scored-category mode (`key = NULL`) has no such
#'   restriction: because `cats` is always a per-item vector, a single call
#'   fully supports a genuinely mixed dichotomous/polytomous test form (see
#'   **Examples**). If raw selected-response data and already-scored
#'   polytomous data need to be analyzed together, score the
#'   selected-response items first (via [irtQ::score_resp()] or this
#'   function's own selected-response mode), then combine the result with
#'   the polytomous scores and call `ctt_distr()` again in scored-category
#'   mode on the combined matrix.
#'
#' @inheritParams ctt_item
#' @param data In selected-response mode (`key` supplied), a data frame or
#'   matrix of raw item responses in the same format as the `data` argument
#'   of [irtQ::score_resp()] (a selected option number, a missing-response
#'   indicator, or a comma-separated string such as `"1,5"` for a
#'   double-marked response). In scored-category mode (`key = NULL`), a data
#'   frame or matrix of already-scored item responses with scores 0 to
#'   `cats[j] - 1` for each item j; an error is raised when a score is not a
#'   whole number between 0 and `cats[j] - 1`. Either way, examinees are rows
#'   and items are columns.
#' @param key `NULL` (default) for scored-category mode, or - to use
#'   selected-response mode - either a vector of correct options in item
#'   order, or a data frame with columns `item` and `key`, exactly as
#'   accepted by [irtQ::score_resp()]. As in [irtQ::score_resp()], each
#'   item's option-coding scheme (numeric, Latin-letter, or general label
#'   such as a Korean syllable label) is inferred independently from its own
#'   key value; see **Details**.
#' @param opt A vector giving the full set of possible response option
#'   values in selected-response mode (e.g., `1:5` for a five-option item,
#'   or `c("A", "B", "C", "D", "E")` for a letter-coded item), applied to
#'   every item. Only used when `key` is supplied; it is silently ignored in
#'   scored-category mode (`key = NULL`), which instead derives its category
#'   set from `cats`. If `NULL` (default), the option set is inferred as the
#'   sorted set of all single (non-blank, non-double-marked) response values
#'   observed anywhere in `data`, using each item's own option-coding scheme
#'   to decide what counts as a single response (see **Details**). Because
#'   `opt` is shared across every item, `ctt_distr()` assumes a test form
#'   uses one consistent option-coding scheme throughout; see **Details**
#'   for what happens if items actually differ.
#' @param cats A numeric vector giving the number of score categories per
#'   item, following the same convention as [ctt()]. Only used in
#'   scored-category mode (`key = NULL`); silently ignored when `key` is
#'   supplied. If `NULL` (default), inferred per item as the observed maximum
#'   score plus one, with a minimum of two. Each value must be a whole number
#'   of at least 2.
#' @param total An optional numeric vector of total scores, one per
#'   examinee, aligned with the rows of `data`. If `NULL` (default), the
#'   total score is computed internally: via [irtQ::score_resp()] in
#'   selected-response mode, or as `rowSums(data)` in scored-category mode
#'   (see **Details** for how `total` interacts with `missing`/listwise
#'   deletion in scored-category mode).
#' @param missing A value indicating missing responses, analogous to the
#'   `missing` argument in [irtQ::est_irt()] and [irtQ::score_resp()]. Its
#'   effect differs by mode (see **Details**): in selected-response mode, a
#'   cell equal to `missing` is recoded to `NA` and then treated as a blank
#'   (omitted) response, contributing to that item's omission rate, without
#'   removing the examinee from the analysis. In scored-category mode, a
#'   cell equal to `missing` is recoded to `NA` and the examinee's row is
#'   then excluded listwise, as in [ctt()]. Default is `NA`.
#' @param crit.distractor A single numeric value; in selected-response mode,
#'   a distractor (non-key option) is flagged when its point-biserial
#'   correlation with the total score exceeds this value. Default is `0`.
#'
#' @details
#' In selected-response mode, each item's option-coding scheme is inferred
#' from its own `key` value, using exactly the same three-way rule as
#' [irtQ::score_resp()]: a key value that parses as a number is compared
#' numerically; a key value consisting of one or more Latin letters is
#' compared as letters, case-insensitively, and a response must itself
#' consist only of Latin letters to count as a single valid option; any
#' other key value (e.g., a Korean syllable label, a Roman numeral written
#' in a non-Latin numeral script, or a circled-number symbol) is treated as
#' a general option label, under which any non-blank, non-double-marked
#' response is accepted as a single valid option and compared to the key
#' (and to `opt`) with a case-insensitive exact string match. This ensures
#' `ctt_distr()` interprets a given `data`/`key` pair identically to
#' [irtQ::score_resp()]. Because `opt` (the option universe used for every
#' item's frequency table) is a single vector shared across all items,
#' `ctt_distr()` is intended for a test form that uses one consistent
#' option-coding scheme throughout (e.g., every item numbered 1-5, or every
#' item lettered A-E); mixing schemes across items within a single call is
#' not well-supported, since the shared `opt` cannot represent two different
#' option universes at once.
#'
#' The two modes handle missing responses differently because they represent
#' different kinds of data. In selected-response mode, an omitted response is
#' itself a meaningful, common outcome for a selected-response item (the
#' examinee saw the item and chose not to answer it), so omissions are
#' tabulated as their own category rather than causing the examinee to be
#' dropped, exactly as in [irtQ::score_resp()]. In scored-category mode,
#' `data` is assumed to already be a scored response matrix (as used by
#' [ctt()]), where a missing cell instead means the
#' item was not administered to, or not observed for, that examinee, so the
#' examinee's row is excluded listwise for that analysis. When `total` is
#' supplied together with scored-category mode's listwise deletion, `total`
#' is assumed to align with the *original* rows of `data`; it is subset using
#' the same completeness mask before use, so its length need only match the
#' original `nrow(data)`.
#'
#' Point-biserial correlations are computed as the Pearson correlation
#' between a 0/1 indicator (whether the examinee chose that particular
#' option, or scored in that particular category) and the total score. As
#' with the item-level discrimination statistics computed by [ctt()], both
#' the raw (uncorrected, item/option-included) and
#' corrected (item-excluded) versions are always computed and returned as
#' separate columns (`pb_raw`/`pb_corrected`); `correct` only selects which
#' of the two feeds the `crit.distractor` flagging criterion in
#' selected-response mode. A point-biserial correlation is undefined, and
#' reported as `NA`, whenever an option/category was chosen by zero or all
#' examinees (a constant indicator).
#'
#' @return A list with these elements:
#' \item{distr}{A data frame in long format with one row per item x option
#'   (selected-response mode) or item x category (scored-category mode),
#'   containing `item`, `option` (the option number, option label, or score
#'   category; numeric for a numeric-coded or scored-category test form,
#'   character for a letter- or general-label-coded test form),
#'   `is_key` (logical; whether this option is the item's correct answer, or
#'   `NA` in scored-category mode), `freq` (the raw number of examinees
#'   choosing this option/category), `pct` (the same count expressed as a
#'   percentage of examinees), `pb_raw` (raw, uncorrected point-biserial
#'   correlation with the total score), `pb_corrected` (corrected,
#'   item-excluded point-biserial correlation with the total score), and
#'   `flag` (a description of `crit.distractor` violations, based on
#'   whichever of `pb_raw`/`pb_corrected` is selected by `correct`, in
#'   selected-response mode, or `""` otherwise).}
#' \item{omit}{In selected-response mode, a data frame with one row per item
#'   giving the blank/omission and double-marked response percentages
#'   (`pct_blank`, `pct_double`). `NULL` in scored-category mode.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::score_resp()], [ctt()]
#'
#' @examples
#' # Selected-response mode: raw option responses + an answer key
#' raw <- data.frame(
#'   V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
#'   V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
#' )
#' key <- c(1, 4)
#' out <- ctt_distr(data = raw, key = key)
#' out$distr
#' out$omit
#'
#' # Scored-category mode: an already-scored polytomous matrix, no key
#' set.seed(1)
#' scored <- data.frame(matrix(sample(0:2, 8 * 3, replace = TRUE), nrow = 8))
#' out2 <- ctt_distr(data = scored)
#' out2$distr
#'
#' # Scored-category mode also handles a genuinely mixed dichotomous/
#' # polytomous test form in a single call, since `cats` is always a
#' # per-item vector (see Note): item 3 below contributes 3 category rows
#' # while items 1-2 contribute 2 each
#' mixed <- data.frame(
#'   I1 = c(1, 0, 1, 1, 0, 1, 0, 1),   # dichotomous item (cats = 2)
#'   I2 = c(0, 1, 0, 1, 1, 0, 0, 1),   # dichotomous item (cats = 2)
#'   I3 = c(2, 0, 1, 2, 1, 0, 2, 1)    # polytomous item (cats = 3)
#' )
#' out2b <- ctt_distr(data = mixed, cats = c(2, 2, 3))
#' out2b$distr
#'
#' # Selected-response mode with letter-coded (A-E) options: each item's
#' # coding scheme is inferred from its own key value, exactly as in
#' # irtQ::score_resp()
#' raw3 <- data.frame(
#'   V1 = c("A", "B", "A", "A", "C", "B", "A", "D"),
#'   V2 = c("D", "D", "B", "B,D", "D", "A", "D", "D")
#' )
#' key3 <- c("A", "D")
#' out3 <- ctt_distr(data = raw3, key = key3)
#' out3$distr
#'
#' @export
ctt_distr <- function(data, item.id = NULL, key = NULL, opt = NULL,
                       cats = NULL, total = NULL, correct = FALSE,
                       missing = NA, crit.distractor = 0) {

  # coerce to a plain data frame so column-wise access behaves consistently
  # for matrix/tibble/data.frame input alike
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  n_item <- ncol(data)
  if (n_item < 1L) {
    stop("`data` must contain at least one item.", call. = FALSE)
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

  mc_mode <- !is.null(key)   # TRUE = selected-response mode, FALSE = scored-category mode

  if (mc_mode) {

    # ---- Selected-response (raw option / distractor analysis) mode -------

    # convert every column to character up front, then recode a
    # user-specified missing sentinel to NA (mirrors score_resp()); blanks
    # remain in the sample here rather than triggering row deletion, since an
    # omission is itself a meaningful response category in this mode
    data[] <- lapply(data, as.character)
    if (!is.na(missing)) {
      data[data == as.character(missing)] <- NA
    }

    # score the raw responses against the key to get the 0/1 item scores
    # needed both for the total score (if not supplied) and for the
    # "corrected" reference total used in each option's point-biserial.
    # score_resp() also validates `key` here (errors on any blank/NA key
    # value), so by the time execution reaches the re-derivation below,
    # `key` is already known to be valid
    scored_out <- score_resp(data = data, key = key, missing = NA)
    scored <- scored_out$scored

    # resolve the total score: use the caller-supplied vector if given
    # (must align with the original rows of `data`), otherwise compute it
    # from the just-scored responses
    if (is.null(total)) {
      total <- rowSums(scored)
    } else if (length(total) != nrow(data)) {
      stop("length(total) must equal nrow(data).", call. = FALSE)
    }

    # resolve the correct-option key into a trimmed character vector, aligned
    # to item column order - mirrors irtQ::score_resp()'s own key resolution
    # exactly, so both functions agree on what a given `key` means
    if (is.data.frame(key)) {
      # order by the numeric item number so character or factor item columns
      # are matched by item number, as in score_resp()
      item_order <- order(as.integer(as.character(key$item)))
      key_vec <- trimws(as.character(key$key[item_order]))
    } else {
      key_vec <- trimws(as.character(key))
    }

    # classify each item's key value into the same three option-coding
    # schemes used by irtQ::score_resp(): numeric, Latin-letter
    # (case-insensitive), or general label (any other non-blank token, e.g.
    # Hangul syllables, Roman numerals in a non-Latin script, circled-number
    # symbols). This per-item flag drives how that item's column in `data`,
    # and the shared `opt` universe below, are compared further down
    is_key_numeric <- grepl("^[0-9]+(\\.[0-9]+)?$", key_vec)
    is_key_latin <- grepl("^[A-Za-z]+$", key_vec)

    # the full set of possible option values, applied to every item (see
    # @details for the single-shared-scheme assumption this implies). If the
    # caller supplies `opt` explicitly, coerce it to the type implied by the
    # items' coding scheme (numeric if every item is numeric-coded,
    # upper-cased character otherwise) without reordering it. Otherwise,
    # infer it as the sorted set of all single (non-blank, non-double-marked)
    # response tokens observed anywhere in `data`, using each item's own
    # scheme (from `is_key_numeric`/`is_key_latin` above) to decide what
    # counts as a single response - exactly mirroring how
    # irtQ::score_resp() itself decides validity per item
    if (is.null(opt)) {
      all_tokens <- character(0)
      for (j in seq_len(n_item)) {
        chr_j <- trimws(data[[j]])
        is_blank_j <- is.na(chr_j) | chr_j == ""
        is_double_j <- !is_blank_j & grepl(",", chr_j, fixed = TRUE)
        if (is_key_numeric[j]) {
          is_single_j <- !is_blank_j & !is_double_j &
            !is.na(suppressWarnings(as.numeric(chr_j)))
        } else if (is_key_latin[j]) {
          is_single_j <- !is_blank_j & !is_double_j &
            grepl("^[A-Za-z]+$", chr_j)
        } else {
          is_single_j <- !is_blank_j & !is_double_j
        }
        # upper-case before pooling so "a"/"A" collapse to one option; a
        # no-op for digits and non-Latin scripts such as Hangul
        all_tokens <- c(all_tokens, toupper(chr_j[is_single_j]))
      }
      uniq_tokens <- unique(all_tokens)
      if (length(uniq_tokens) > 0L && all(is_key_numeric)) {
        # every item is numeric-coded: sort as numbers and keep `opt`
        # numeric
        opt <- sort(as.numeric(uniq_tokens))
      } else {
        # at least one item is letter- or general-label-coded: sort as
        # character strings and keep `opt` character
        opt <- sort(uniq_tokens)
      }
    } else if (all(is_key_numeric)) {
      opt <- as.numeric(opt)          # keep caller's values/order, coerce type
    } else {
      opt <- toupper(as.character(opt))   # keep caller's order, canonicalize case
    }

    # accumulate one long-format row per item x option, plus one row per
    # item for the omission/double-marking percentages
    distr_rows <- vector("list", n_item * length(opt))
    omit_rows <- vector("list", n_item)
    row_i <- 1L

    for (j in seq_len(n_item)) {

      resp_chr <- trimws(data[[j]])                          # this item's responses
      is_blank <- is.na(resp_chr) | resp_chr == ""            # omitted responses
      is_double <- !is_blank & grepl(",", resp_chr, fixed = TRUE)  # double-marked

      # this item's own option-coding scheme (numeric / Latin-letter /
      # general label), matching irtQ::score_resp()'s per-item detection
      # exactly; `resp_val` is the response in the canonical form used for
      # comparison below (a parsed number for a numeric-coded item, an
      # upper-cased string otherwise), and `is_single` marks which responses
      # count as a single valid option for this item's scheme
      if (is_key_numeric[j]) {
        resp_val <- suppressWarnings(as.numeric(resp_chr))
        is_single <- !is_blank & !is_double & !is.na(resp_val)
      } else if (is_key_latin[j]) {
        resp_val <- toupper(resp_chr)
        is_single <- !is_blank & !is_double & grepl("^[A-Za-z]+$", resp_chr)
      } else {
        resp_val <- toupper(resp_chr)
        is_single <- !is_blank & !is_double
      }

      # this item's key value, in that same canonical form
      key_val_j <- if (is_key_numeric[j]) {
        as.numeric(key_vec[j])
      } else {
        toupper(key_vec[j])
      }

      item_score <- scored[[j]]                       # this item's 0/1 score
      ref_total_raw <- total                            # raw (item-included) reference total
      ref_total_corrected <- total - item_score         # corrected (item-excluded) reference total

      n_examinee <- length(resp_chr)

      for (o in opt) {
        # this particular option value, in the same canonical form as
        # `resp_val`/`key_val_j` above, so the comparisons below are
        # apples-to-apples regardless of the item's coding scheme
        o_val <- if (is_key_numeric[j]) as.numeric(o) else toupper(as.character(o))

        chosen <- is_single & resp_val == o_val
        freq <- as.integer(sum(chosen))    # raw count of examinees choosing option o
        pct <- 100 * freq / n_examinee     # same count, expressed as a percentage
        chosen_num <- as.numeric(chosen)
        sd_chosen_ok <- stats::sd(chosen) > 0

        # point-biserial correlations of "chose option o" with the raw and
        # corrected reference totals; NA when the indicator is constant
        # (chosen by none or all), regardless of `correct`
        pb_raw <- if (sd_chosen_ok && stats::sd(ref_total_raw) > 0) {
          stats::cor(chosen_num, ref_total_raw)
        } else {
          NA_real_
        }
        pb_corrected <- if (sd_chosen_ok && stats::sd(ref_total_corrected) > 0) {
          stats::cor(chosen_num, ref_total_corrected)
        } else {
          NA_real_
        }
        pb_for_flag <- if (correct) pb_corrected else pb_raw  # drives crit.distractor below

        is_key <- isTRUE(o_val == key_val_j)   # is this option the correct answer?

        # flag a distractor (non-key option) whose point-biserial correlation
        # with the total score exceeds crit.distractor (attracts high scorers)
        flag_txt <- if (!is_key && !is.na(pb_for_flag) && pb_for_flag > crit.distractor) {
          "distractor more attractive to high scorers"
        } else {
          ""
        }

        distr_rows[[row_i]] <- data.frame(
          item = item_names[j], option = o, is_key = is_key,
          freq = freq, pct = round(pct, 2), pb_raw = round(pb_raw, 3),
          pb_corrected = round(pb_corrected, 3), flag = flag_txt,
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1L
      }

      # per-item omission/double-marking percentages
      omit_rows[[j]] <- data.frame(
        item = item_names[j],
        pct_blank = round(100 * sum(is_blank) / n_examinee, 2),
        pct_double = round(100 * sum(is_double) / n_examinee, 2),
        stringsAsFactors = FALSE
      )
    }

    distr <- do.call(rbind, distr_rows)
    omit <- do.call(rbind, omit_rows)

  } else {

    # ---- Scored-category (already-scored polytomous) mode ----------------

    # recode a user-specified missing sentinel to NA, then listwise-delete
    # any examinee with a remaining missing response, as in ctt_item()
    if (!is.na(missing)) {
      data[data == missing] <- NA
    }
    complete_rows <- stats::complete.cases(data)
    n_dropped <- sum(!complete_rows)
    if (n_dropped > 0L) {
      warning(n_dropped, " examinee(s) with missing item responses were ",
              "excluded listwise from ctt_distr().", call. = FALSE)
    }
    data <- data[complete_rows, , drop = FALSE]

    # stop when a score or a `cats` value is not valid
    check_ctt_scores(data, cats, item_names)

    # resolve the total score: use the caller-supplied vector (subset to the
    # same complete rows, so it stays aligned with `data` after deletion) if
    # given, otherwise compute it as the row sums of the scored data
    if (is.null(total)) {
      total <- rowSums(data)
    } else {
      if (length(total) != length(complete_rows)) {
        stop("length(total) must equal the original nrow(data).",
             call. = FALSE)
      }
      total <- total[complete_rows]
    }

    # infer the number of score categories per item, when not supplied
    # explicitly (cats[j] = observed max score + 1, with a minimum of two
    # categories), as in ctt_item()
    if (is.null(cats)) {
      cats <- pmax(vapply(data, function(x) max(x) + 1, numeric(1)), 2)
    }

    n_examinee <- nrow(data)
    distr_rows <- vector("list", sum(cats))
    row_i <- 1L

    for (j in seq_len(n_item)) {

      item_score <- data[[j]]                              # this item's score
      ref_total_raw <- total                                 # raw (item-included) reference total
      ref_total_corrected <- total - item_score              # corrected (item-excluded) reference total
      categories <- seq(0, cats[j] - 1)                      # 0, 1, ..., cats[j]-1

      for (o in categories) {
        chosen <- item_score == o
        freq <- as.integer(sum(chosen))    # raw count of examinees in category o
        pct <- 100 * freq / n_examinee     # same count, expressed as a percentage
        chosen_num <- as.numeric(chosen)
        sd_chosen_ok <- stats::sd(chosen) > 0

        # point-biserial correlations against both the raw and corrected
        # reference totals; both are always computed (see @details)
        pb_raw <- if (sd_chosen_ok && stats::sd(ref_total_raw) > 0) {
          stats::cor(chosen_num, ref_total_raw)
        } else {
          NA_real_
        }
        pb_corrected <- if (sd_chosen_ok && stats::sd(ref_total_corrected) > 0) {
          stats::cor(chosen_num, ref_total_corrected)
        } else {
          NA_real_
        }

        distr_rows[[row_i]] <- data.frame(
          item = item_names[j], option = o, is_key = NA,
          freq = freq, pct = round(pct, 2), pb_raw = round(pb_raw, 3),
          pb_corrected = round(pb_corrected, 3), flag = "",
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1L
      }
    }

    distr <- do.call(rbind, distr_rows)
    omit <- NULL   # no omission/double-marking concept in scored-category mode
  }

  list(distr = distr, omit = omit)
}
