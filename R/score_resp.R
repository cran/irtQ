#' Score Selected-Response Item Data Against an Answer Key
#'
#' This function converts raw selected-response (e.g., multiple-choice) item
#' data - coded as the selected option (an option number, e.g., 1-5 for a
#' five-option item; a single-letter option label, e.g., `"A"`-`"E"`; or any
#' other option label in a different script or notation, such as Korean
#' syllable labels or Roman numerals), a missing-response indicator for an
#' omitted response, or a comma-separated string (e.g., `"1,5"` or `"A,C"`)
#' for a double-marked (multiple-option) response - into a dichotomously
#' scored item-response matrix (0 = incorrect, 1 = correct) using a supplied
#' answer key. The option coding scheme is detected independently for each
#' item from its own `key` value, so a single test form may freely mix
#' numerically coded, letter-coded, and other label-coded items. Letter
#' (and other case-bearing) option labels are matched case-insensitively
#' (e.g., a response of `"a"` matches a key of `"A"`). Omitted and
#' double-marked responses are both scored as incorrect (0), and their
#' frequencies are tabulated separately for reporting.
#'
#' @note This function applies to **dichotomous, single-key, selected-response
#'   items only** (i.e., items with exactly one correct option, scored into
#'   exactly two categories: 0 = incorrect and 1 = correct). It does not
#'   generalize to items with more than two score categories (e.g.,
#'   partial-credit polytomous items scored 0, 1, 2, ...), since a
#'   selected-response item cannot yield a partial-credit score by
#'   construction. If polytomous item response data are available, supply an
#'   already-scored response matrix (with categories 0, 1, 2, ..., cats - 1)
#'   directly to downstream item-analysis functions instead of using
#'   `score_resp()`. See **Details** for further discussion.
#'
#' @param data A data frame or matrix of raw item responses, with examinees in
#'   rows and items in columns (in the same left-to-right order as `key`).
#'   Do not include non-item columns (e.g., group or examinee identifiers);
#'   subset those out before calling `score_resp()`. Each cell should contain
#'   either a single selected option - an option number (e.g., `2`), a
#'   single-letter option label (e.g., `"C"`), or any other option label
#'   (e.g., a Korean syllable label, or a Roman numeral such as `"II"`) -
#'   matching whichever coding scheme is used for that item's `key`
#'   value, a missing-response indicator (see `missing`), or a character
#'   string of comma-separated options (e.g., `"1,5"` or `"A,C"`) for a
#'   double-marked response. To keep double-marked responses as character
#'   strings rather than having them coerced to `NA` on import, read the
#'   source file with all item columns imported as text (e.g.,
#'   `readxl::read_excel(..., col_types = "text")`) before passing `data` to
#'   this function.
#' @param key Either (a) a vector of correct options in item order (length
#'   must equal `ncol(data)`), where each element is a numeric/integer option
#'   number (e.g., `4`), a single-letter option label (e.g., `"D"`), or any
#'   other non-blank option label (e.g., a Korean syllable label, or a
#'   Roman numeral such as `"II"`), or (b) a data frame
#'   with columns `item` and `key` giving the item number and its correct
#'   option, respectively (the `key` column may likewise mix option formats
#'   across rows). Each item's coding scheme is inferred independently from
#'   its own key value, so items using different option formats may be
#'   freely mixed within the same vector or data frame; letter (and other
#'   case-bearing) values are matched case-insensitively. When a data frame
#'   is supplied, it is internally sorted by `item` before use, so its row
#'   order does not need to match the column order of `data`.
#' @param missing A value indicating missing (omitted) responses in `data`,
#'   analogous to the `missing` argument in [irtQ::est_irt()] and
#'   [irtQ::est_score()]. Any cell equal to `missing` is recoded to `NA`
#'   before scoring. Default is `NA`, meaning `data` is assumed to already use
#'   `NA` (or an empty string) for omitted responses, and no recoding is
#'   performed. Set this to whatever sentinel value a particular data set uses
#'   for a missing response (e.g., `-9`, `"9"`) when it differs from `NA`.
#'
#' @details
#' For each item, a response is scored as correct (1) only when it is a
#' single, non-missing option that matches the corresponding value in `key`.
#' All other cases - an omitted response, a double-marked response, or a
#' single but incorrect option - are scored as incorrect (0). Missing and
#' double-marked responses are tallied separately in `resp_summary` so that
#' omission and double-marking rates can be reported independently, even
#' though both are scored as 0.
#'
#' Each item's option coding scheme is determined independently from its own
#' `key` value, using one of three rules:
#' \itemize{
#'   \item If the key value parses as a number (e.g., `4`), responses to
#'     that item are compared numerically.
#'   \item If the key value consists of one or more Latin letters (e.g.,
#'     `"D"`), responses are compared as letters, case-insensitively (`"d"`
#'     and `"D"` are treated as the same option), and a response must itself
#'     consist only of Latin letters to count as a single valid option (see
#'     below for what happens otherwise).
#'   \item Otherwise (e.g., a Korean syllable label, a Roman numeral written
#'     with non-Latin numeral characters, a circled-number symbol, or any
#'     other non-numeric, non-Latin-letter label), the key value is treated
#'     as a general option label: any
#'     non-blank, non-double-marked response to that item is accepted as a
#'     single valid option, and is compared to the key with a
#'     case-insensitive (`toupper()`-based, which only affects any embedded
#'     Latin characters) exact string match. There is no universal,
#'     script-independent rule for "a well-formed option label" analogous to
#'     the numeric or Latin-letter checks above, so this scheme cannot
#'     distinguish a genuinely incorrect option from a garbled response
#'     token; both are simply scored 0, and `n_invalid` is always 0 for
#'     items using this scheme (see below).
#' }
#' This per-item detection means a single test form may freely mix items
#' using any of the three schemes - only the value supplied in `key` for a
#' given item determines how that item's column in `data` is interpreted. A
#' `key` value that is blank or `NA` causes `score_resp()` to stop with an
#' error, since no correct option was actually supplied for that item.
#'
#' A response cell is classified as double-marked when, after coercion to
#' character, it contains a comma (e.g., `"1,5"`, `"B,D"`). A response cell is
#' classified as blank/omitted when it equals `missing` (recoded to `NA`
#' beforehand), or is otherwise `NA` or an empty/whitespace-only string. For
#' numeric- and Latin-letter-coded items, any response that is none of blank,
#' double-marked, or a single valid option in the coding scheme used for that
#' item (e.g., a numeric token where a letter was expected for a
#' letter-coded item, or stray text from a data-entry artifact) is
#' classified as invalid; it is still scored 0, but tallied separately in
#' `n_invalid` (with a `warning()` raised) so it is visible rather than
#' silently merged into the "wrong" count. As noted above, general
#' label-coded items never populate `n_invalid`, since any non-blank,
#' non-double-marked token is accepted as a single valid (though possibly
#' incorrect) option for that scheme. Note that `n_correct + n_wrong` always
#' equals `n`, but `n_blank`, `n_double`, and `n_invalid` are reported as
#' separate diagnostic tallies, not as a four-way disjoint partition of
#' `n_wrong`.
#'
#' This function performs dichotomous (single-key, selected-response) scoring
#' only, since a selected-response item has exactly one correct option by
#' design and therefore only two possible score categories (correct/incorrect).
#' It is not generalized to items with more than two score categories (e.g.,
#' partial-credit polytomous items). If such data become available, they
#' should be supplied directly, already scored, to downstream item-analysis
#' functions that support polytomous data, bypassing `score_resp()`.
#'
#' @return A list with two elements:
#' \item{scored}{A data frame with the same dimensions as `data`, containing
#'   the dichotomously scored (0/1) item responses. Column names follow
#'   `colnames(data)` (or `V1, V2, ...` when `data` has no column names).}
#' \item{resp_summary}{A data frame with one row per item, reporting the
#'   number and percentage of blank and double-marked responses, the number
#'   of unrecognized/invalid response tokens (`n_invalid`), and the number
#'   of examinees scored correct/incorrect, for each item.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [ctt_distr()], [irtQ::est_irt()], [irtQ::est_score()]
#'
#' @examples
#' # A small 5-examinee, 3-item toy example with a blank and a double-mark
#' raw <- data.frame(
#'   V1 = c("1", "2", NA, "1", "3"),
#'   V2 = c("4", "4", "4", "2,4", "4"),
#'   V3 = c("2", "1", "3", "3", NA)
#' )
#' key <- c(1, 4, 3)
#' out <- score_resp(data = raw, key = key)
#' out$scored
#' out$resp_summary
#'
#' # Same data, but omitted responses are coded as "9" instead of NA
#' raw2 <- data.frame(
#'   V1 = c("1", "2", "9", "1", "3"),
#'   V2 = c("4", "4", "4", "2,4", "4"),
#'   V3 = c("2", "1", "3", "3", "9")
#' )
#' out2 <- score_resp(data = raw2, key = key, missing = "9")
#' identical(out$scored, out2$scored)
#'
#' # A letter-coded (A-E) example, including a double-marked response and a
#' # lowercase response ("a"), which is matched case-insensitively
#' raw3 <- data.frame(
#'   V1 = c("A", "B", NA, "a", "C"),
#'   V2 = c("D", "D", "D", "B,D", "D"),
#'   V3 = c("B", "A", "C", "C", NA)
#' )
#' key3 <- c("A", "D", "C")
#' out3 <- score_resp(data = raw3, key = key3)
#' out3$scored
#'
#' # A mixed-format test: item 1 is numerically coded, items 2-3 are
#' # letter-coded; the coding scheme is inferred separately for each item
#' # from its own key value
#' raw4 <- data.frame(
#'   V1 = c("1", "2", "1", "1", "3"),
#'   V2 = c("D", "D", "D", "B,D", "D"),
#'   V3 = c("B", "A", "C", "C", NA)
#' )
#' key4 <- c(1, "D", "C")
#' out4 <- score_resp(data = raw4, key = key4)
#' out4$scored
#'
#' # A full five-option (A-E) multiple-choice item, scored against key "C"
#' raw5 <- data.frame(
#'   V1 = c("A", "B", "C", "D", "E", "c")
#' )
#' key5 <- c("C")
#' out5 <- score_resp(data = raw5, key = key5)
#' out5$scored
#'
#' # A general-label-coded example using non-Latin option labels: the Korean
#' # syllables romanized "ga"/"na"/"da"/"ra"/"ma", written here with \uXXXX
#' # escapes for ASCII portability. Because these labels are neither numbers
#' # nor Latin letters, they are handled by the general label-coding scheme
#' # described in Details
#' raw6 <- data.frame(
#'   V1 = c("\uAC00", "\uB098", "\uB2E4", "\uB2E4", NA) # ga, na, da, da, NA
#' )
#' key6 <- c("\uB2E4") # da
#' out6 <- score_resp(data = raw6, key = key6)
#' out6$scored
#'
#' @export
score_resp <- function(data, key, missing = NA) {

  # coerce the input responses to a plain data frame so that column-wise
  # extraction (data[[j]]) behaves consistently for matrix/tibble/data.frame
  # input alike
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  # number of examinees (rows) and items (columns) in the raw response data
  n_examinee <- nrow(data)
  n_item <- ncol(data)

  # basic guard: need at least one examinee row and one item column to score
  if (n_examinee == 0L || n_item == 0L) {
    stop("`data` must have at least one row and one column.", call. = FALSE)
  }

  # convert every column to character up front so that the `missing` sentinel
  # comparison below, and the per-item blank/double-mark classification
  # further down, behave consistently regardless of whether a column was
  # originally numeric or character (a double-marked cell such as "1,5" must
  # survive as text either way)
  data[] <- lapply(data, as.character)

  # recode a user-specified missing-value sentinel (e.g., -9, "9") to NA
  # before classifying responses; this mirrors the `missing` argument
  # convention used by irtQ::est_irt()/irtQ::est_score(). The default
  # `missing = NA` means the input already uses R's native NA (or an empty
  # string) for omitted responses, so nothing is recoded here
  if (!is.na(missing)) {
    data[data == as.character(missing)] <- NA
  }

  # resolve the `key` argument into a character vector of correct options
  # (trimmed strings), ordered to match the column order of `data`. the key
  # is deliberately kept as character here (rather than coerced with
  # as.numeric()) so that each item's correct option may be either a numeric
  # option number (e.g., 4) or a single-letter option label (e.g., "D"); the
  # coding scheme actually used for scoring is inferred separately for each
  # item further below, from that item's own key value
  if (is.data.frame(key)) {

    # a key data frame must supply both an item index column and a key column
    if (!all(c("item", "key") %in% names(key))) {
      stop("When `key` is a data frame, it must contain columns ",
           "'item' and 'key' (one row per item, giving the item number and ",
           "its correct option).", call. = FALSE)
    }

    # sort by item number so the row order in the file need not match the
    # column order of `data`
    key <- key[order(key$item), ]
    key_vec <- trimws(as.character(key$key))

    # the key's item numbers must be exactly 1, 2, ..., n_item with no
    # duplicates and no gaps; otherwise key values could silently misalign
    # with data columns after sorting (e.g., a duplicated or skipped item
    # number in a mistyped key would go undetected without this check).
    # compare as integer (not identical() on raw `key$item`) so that numeric
    # vs. integer storage mode never causes a false mismatch; only the
    # actual item numbers matter here
    if (!identical(sort(as.integer(key$item)), seq_len(n_item))) {
      stop("`key$item` must contain exactly the integers 1:ncol(data), ",
           "with no duplicates or gaps. Check the answer key for a ",
           "mistyped or missing item number.", call. = FALSE)
    }

  } else {

    # `key` was already supplied as a plain atomic vector of correct options;
    # keep it as trimmed character so numeric option numbers and letter
    # option labels can be freely mixed across items
    key_vec <- trimws(as.character(key))
  }

  # the number of key values must match the number of item columns in `data`
  if (length(key_vec) != n_item) {
    stop("length(key) must equal ncol(data): one key value per item column.",
         call. = FALSE)
  }

  # a key value must not be blank/NA - that is the only requirement on `key`
  # itself, since every non-blank value (a number, a Latin letter, a Korean
  # syllable, a Roman numeral, or any other label) is handled as a valid
  # option identifier by the per-item scheme detection below
  if (any(is.na(key_vec) | key_vec == "")) {
    bad_pos <- which(is.na(key_vec) | key_vec == "")
    stop("`key` contains missing or blank value(s) at position(s): ",
         paste(bad_pos, collapse = ", "), ". Every item must have a ",
         "non-blank correct-option value.", call. = FALSE)
  }

  # classify each item's key value into one of three option-coding schemes,
  # used to select the comparison rule for that item's column in the
  # scoring loop below; this per-item flagging is what allows a single test
  # form to freely mix items using different option formats:
  #   - numeric-coded (e.g., "4"): compared as numbers
  #   - Latin-letter-coded (e.g., "D"): compared as letters, case-insensitive
  #   - anything else (e.g., a Korean syllable label, a Roman numeral
  #     written in a non-Latin numeral script, a circled-number symbol):
  #     compared as a general label (see Details for how this scheme
  #     differs from the other two)
  is_key_numeric <- grepl("^[0-9]+(\\.[0-9]+)?$", key_vec)
  is_key_latin <- grepl("^[A-Za-z]+$", key_vec)

  # column names to carry over to the scored output; fall back to V1, V2, ...
  # when the input has no column names (mirrors irtQ's default item labeling)
  item_names <- colnames(data)
  if (is.null(item_names)) item_names <- paste0("V", seq_len(n_item))

  # pre-allocate the scored (0/1) matrix and a list to collect per-item
  # summary rows (filled in the loop below)
  scored <- matrix(0L, nrow = n_examinee, ncol = n_item,
                    dimnames = list(NULL, item_names))
  resp_summary <- vector("list", n_item)

  # loop over items; each item is scored independently against its own key
  # value, since option coding and correctness are item-specific
  for (j in seq_len(n_item)) {

    # this item's raw responses, already character (and missing-recoded)
    # thanks to the up-front conversion above; trim incidental whitespace
    resp_chr <- trimws(data[[j]])

    # a response is blank/omitted when it is NA or an empty string once
    # surrounding whitespace is removed
    is_blank <- is.na(resp_chr) | resp_chr == ""

    # a response is double-marked when it contains a comma (e.g., "1,5" or
    # "A,C"), since single valid responses never contain one
    is_double <- !is_blank & grepl(",", resp_chr, fixed = TRUE)

    # a single valid response is defined by this item's own coding scheme
    # (determined above from its key value): see the three-way rule
    # described in Details
    item_score <- integer(n_examinee)
    if (is_key_numeric[j]) {

      # numeric-coded item: parse each response as a number and compare
      # numerically to the key
      resp_num <- suppressWarnings(as.numeric(resp_chr))
      is_single <- !is_blank & !is_double & !is.na(resp_num)
      item_score[is_single & resp_num == as.numeric(key_vec[j])] <- 1L

    } else if (is_key_latin[j]) {

      # Latin-letter-coded item: a single valid response is one or more
      # Latin letters; compare to the key case-insensitively (e.g., "a"
      # matches a key of "D" only if they are literally the same letter once
      # both are upper-cased)
      is_single <- !is_blank & !is_double & grepl("^[A-Za-z]+$", resp_chr)
      item_score[is_single & toupper(resp_chr) == toupper(key_vec[j])] <- 1L

    } else {

      # general label-coded item (e.g., Korean syllable labels, Roman
      # numerals in a non-Latin numeral script, circled-number symbols, or
      # any other non-numeric, non-Latin-letter label): there is no
      # universal, script-independent format check for "a well-formed
      # option label" the way there is for numbers or Latin letters, so
      # every non-blank, non-double-marked token is accepted as a single
      # valid response for this item; toupper() is applied for consistency
      # with the other two schemes, but is a no-op for non-Latin scripts
      # such as Hangul, so this effectively reduces to an exact string match
      is_single <- !is_blank & !is_double
      item_score[is_single & toupper(resp_chr) == toupper(key_vec[j])] <- 1L
    }

    # any response that is none of blank/double/single-valid is an
    # unexpected token (e.g., a numeric answer on a letter-coded item, or
    # stray text from a data-entry or OCR artifact); it is still scored 0
    # above, but is tallied separately so it is visible rather than silently
    # folded into "wrong" with no trace
    is_invalid <- !is_blank & !is_double & !is_single

    scored[, j] <- item_score

    # tabulate this item's blank/double-mark/invalid/correctness counts
    n_blank <- sum(is_blank)               # count of omitted responses
    n_double <- sum(is_double)             # count of double-marked responses
    n_invalid <- sum(is_invalid)           # count of unrecognized tokens
    n_correct <- sum(item_score)           # count scored correct (1)

    # warn (once per item) if an unrecognized token shows up, since valid
    # data is expected to contain only option numbers, missing indicators, or
    # comma strings; anything else likely signals a data-quality issue
    if (n_invalid > 0L) {
      warning("Item '", item_names[j], "' has ", n_invalid,
              " response(s) that are neither a valid option number, blank, ",
              "nor double-marked; these are scored 0 but flagged as ",
              "n_invalid in resp_summary.", call. = FALSE)
    }

    resp_summary[[j]] <- data.frame(
      item = item_names[j],                            # item label
      key = key_vec[j],                                 # correct option
      n = n_examinee,                                    # examinees scored
      n_correct = n_correct,                             # number correct
      n_wrong = n_examinee - n_correct,                  # number incorrect
      n_blank = n_blank,                                 # number blank
      n_double = n_double,                               # number double-marked
      n_invalid = n_invalid,                             # number unrecognized
      pct_blank = round(100 * n_blank / n_examinee, 2),  # blank rate (%)
      pct_double = round(100 * n_double / n_examinee, 2),# double-mark rate (%)
      stringsAsFactors = FALSE
    )
  }

  # row-bind the per-item summary rows into a single data frame
  resp_summary <- do.call(rbind, resp_summary)

  # convert the scored matrix to a data frame, matching irtQ's convention of
  # returning data frames (e.g., as consumed by est_irt()/est_score())
  scored <- as.data.frame(scored)

  # return both the scored response matrix and the per-item response summary
  list(scored = scored, resp_summary = resp_summary)
}
