#' Frequency Distribution Table for Total Scores
#'
#' Computes a frequency distribution table for a vector of total (raw)
#' scores - frequency, percentage, and cumulative percentage for each score
#' value - commonly reported alongside classical test theory (CTT) item
#' analysis results. This is the same function [ctt()] calls internally to
#' build the total-score frequency distribution included in its output, but
#' it is also exported and fully usable on its own for any vector of integer
#' total scores.
#'
#' @param score A numeric vector of total (raw) scores, one value per
#'   examinee.
#' @param missing A value indicating missing scores in `score`, analogous to
#'   the `missing` argument in [irtQ::est_irt()] and [irtQ::score_resp()]. Any
#'   element equal to `missing` is recoded to `NA` before tabulation. Default
#'   is `NA`. Examinees with a (remaining) missing score are excluded from
#'   the table, with a warning reporting how many were dropped.
#'
#' @details
#' The table includes every integer score value spanning the observed
#' minimum to maximum score (not just the values that were actually
#' observed), so that a score with zero examinees still appears in the table
#' with a frequency of 0, matching how a raw-score frequency table is
#' conventionally reported. Percentages are computed relative to the number
#' of non-missing scores and rounded to two decimal places; cumulative
#' percentages are the running sum of the unrounded percentages, rounded to
#' two decimal places only in the final output, so that rounding error does
#' not accumulate across rows and the final cumulative percentage totals
#' almost exactly 100 (subject only to ordinary rounding).
#'
#' This function assumes scores already lie on an integer (or otherwise
#' evenly spaced discrete) scale, as is standard for a raw total score; it
#' does not bin or group continuous values, and raises an error if any
#' non-integer score is supplied rather than silently dropping it.
#'
#' The output table always spans `min(score)` to `max(score)`, so its size
#' scales with the observed score *range*, not the sample size; a single
#' unusually large or small outlier score will produce a correspondingly
#' large table.
#'
#' @return A data frame with one row per score value from `min(score)` to
#'   `max(score)`, containing:
#' \item{score}{the score value.}
#' \item{freq}{the number of examinees with that score.}
#' \item{pct}{the percentage of examinees with that score.}
#' \item{cum_pct}{the cumulative percentage up to and including that score.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [ctt()]
#'
#' @examples
#' set.seed(1)
#' score <- rbinom(500, size = 20, prob = 0.6)
#' freq_score(score)
#'
#' @export
freq_score <- function(score, missing = NA) {

  # coerce to a plain numeric vector so factor/character input does not
  # silently break the arithmetic below
  score <- as.numeric(score)

  # recode a user-specified missing-value sentinel to NA before tabulation,
  # mirroring the `missing` argument convention used elsewhere in irtQ
  if (!is.na(missing)) {
    score[score == missing] <- NA
  }

  # drop any missing scores, warning how many were excluded
  n_dropped <- sum(is.na(score))
  if (n_dropped > 0L) {
    warning(n_dropped, " missing score(s) were excluded from freq_score().",
            call. = FALSE)
    score <- score[!is.na(score)]
  }

  # need at least one valid score to build a table
  if (length(score) == 0L) {
    stop("`score` has no non-missing values to tabulate.", call. = FALSE)
  }

  # this function only supports integer-valued scores (see @details); a
  # non-integer score would silently fall outside the seq.int()-based bins
  # below and be dropped without any other warning, so check explicitly here
  if (!isTRUE(all.equal(score, round(score)))) {
    stop("`score` must contain only integer-valued scores; freq_score() ",
         "does not bin or group non-integer values.", call. = FALSE)
  }

  # full range of integer score values (including any with zero observed
  # frequency), from the observed minimum to the observed maximum
  score_range <- seq.int(min(score), max(score))

  # count of examinees at each value in the full range; matching against
  # score_range (rather than using table() alone) ensures unobserved values
  # in the middle of the range still appear with frequency 0
  freq <- vapply(score_range, function(s) sum(score == s), integer(1))

  # percentage of (non-missing) examinees at each score value
  n_total <- length(score)
  pct <- 100 * freq / n_total

  # cumulative percentage: running sum of the unrounded percentages, with
  # rounding applied only once at the very end (see @details)
  cum_pct <- cumsum(pct)

  # assemble the frequency table, rounding percentages for readable reporting
  data.frame(
    score = score_range,
    freq = freq,
    pct = round(pct, 2),
    cum_pct = round(cum_pct, 2)
  )
}
