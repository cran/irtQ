#' Classical Test Theory (CTT) Item and Test Analysis
#'
#' Computes classical test theory (CTT) statistics from a scored item
#' response data set (dichotomous or polytomous). For each item, computes
#' its difficulty, an item-total correlation reflecting how well the item
#' discriminates between high- and low-scoring examinees, and Cronbach's
#' alpha recomputed with that item removed, with optional flagging of items
#' whose difficulty or discrimination falls outside conventional quality
#' thresholds. At the test level, computes Cronbach's alpha as a reliability
#' estimate (in both raw and standardized forms), the associated standard
#' error of measurement (SEM), and the average item difficulty and
#' discrimination across the test. The distribution of total (raw) scores
#' across examinees is also tabulated. The result is returned as an object of
#' class `"ctt"`, with [print.ctt()] and [summary.ctt()] methods that display
#' a condensed and a full report, respectively, analogous to how
#' [irtQ::est_irt()] pairs with `print.est_irt()`/`summary.est_irt()`.
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
#'   where an item is correlated with the total score that includes its own
#'   contribution - and the corrected item-total correlation - excluding its
#'   own contribution - are always computed and reported as separate columns.
#'   This argument only controls which of the two is used for flagging: if
#'   `TRUE`, flagging is based on the corrected value; if `FALSE` (default),
#'   on the raw value.
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
#' Difficulty for item j is the mean observed item score divided by the
#' item's maximum possible score,
#' \deqn{p_j = \frac{\bar{X}_j}{m_j},}
#' where \eqn{m_j} is `cats[j] - 1`. For a dichotomous item (`cats[j] = 2`)
#' this reduces to the familiar proportion-correct difficulty index; for a
#' polytomous item it expresses the average score as a proportion of the
#' maximum attainable score, keeping difficulty on a common 0-1 scale
#' regardless of the number of score categories.
#'
#' Discrimination for item j is the Pearson correlation between the item
#' score and the total score, always computed and reported in two forms: the
#' raw (uncorrected) item-total correlation, which correlates the item with
#' the total score that includes the item's own contribution, and the
#' corrected (item-excluded) item-total correlation, which excludes it. For a
#' dichotomous item, the raw item-total correlation is mathematically
#' equivalent to the point-biserial correlation. `correct` selects which of
#' the two feeds the discrimination flagging criterion (`crit.dis`); both are
#' always reported as separate columns regardless of `correct` (this mirrors
#' how, e.g., `psych::alpha()` reports the raw item-total correlation as
#' `raw.r` and the corrected version as `r.drop`).
#'
#' Alpha-with-item-removed for item j is Cronbach's alpha (see below)
#' recomputed using only the remaining items, so that a value above the
#' overall alpha identifies an item whose removal would increase the
#' reliability of the test.
#'
#' At the test level, two forms of Cronbach's alpha are always computed and
#' reported. Raw alpha uses the standard variance-based formula
#' \deqn{\alpha = \frac{k}{k - 1}\left(1 - \frac{\sum_{j} \sigma_j^2}{\sigma_X^2}\right),}
#' where \eqn{k} is the number of items, \eqn{\sigma_j^2} is the variance of
#' item j, and \eqn{\sigma_X^2} is the variance of the total score X. Raw
#' alpha reflects the reliability of the actual (unweighted) total score
#' obtained by simply summing the item scores - the score most tests
#' actually use for reporting and decisions. Standardized alpha instead first
#' standardizes every item to unit variance before combining them,
#' \deqn{\alpha_{std} = \frac{k \bar{r}}{1 + (k - 1)\bar{r}},}
#' where \eqn{\bar{r}} is the average pairwise correlation among all items.
#' Standardized alpha differs meaningfully from raw alpha mainly when items
#' vary substantially in scale or format (e.g., a mix of dichotomous and
#' polytomous items with very different score ranges); when all items share
#' the same scale and format, the two are typically close, and raw alpha
#' remains the more directly interpretable of the two since it matches the
#' reliability of the score actually used in practice. See Cronbach (1951)
#' for the original derivation of coefficient alpha, and Osburn (2000) for a
#' discussion contrasting the raw (covariance-based) and standardized
#' (correlation-based) forms.
#'
#' The standard error of measurement (SEM) is
#' \eqn{SEM = SD(X)\sqrt{1 - \alpha}} (using raw alpha), the standard CTT
#' relationship between test reliability and measurement precision.
#'
#' The total-score frequency distribution is computed internally by
#' [freq_score()], using the row sums of `data` after applying the same
#' `missing` recoding and listwise deletion of incomplete rows used
#' throughout the rest of this function, so that the frequency distribution
#' reflects exactly the same set of examinees and the same total-score
#' definition used elsewhere in the analysis. [freq_score()] is also
#' exported separately, for callers who only need a frequency distribution
#' for an arbitrary vector of integer total scores.
#'
#' @return An object of class `"ctt"`, a list with elements:
#' \item{item}{A data frame with one row per item, containing the item
#'   label, number of score categories, difficulty, the raw (uncorrected)
#'   item-total correlation (`discrimination_raw`), the corrected
#'   (item-excluded) item-total correlation (`discrimination_corrected`),
#'   alpha-with-item-removed, and, if `flag = TRUE`, a `flag` column.}
#' \item{crit}{A list echoing the `crit.p` and `crit.dis` thresholds used for
#'   flagging.}
#' \item{alpha}{A one-row test-level summary data frame containing
#'   `n_examinee`, `n_item`, `alpha`, `alpha_std`, `sem`, `mean_difficulty`,
#'   `mean_discrimination_raw`, and `mean_discrimination_corrected`. Both the
#'   raw (`alpha`) and standardized (`alpha_std`) forms of Cronbach's alpha
#'   are always included; standardized alpha standardizes every item to unit
#'   variance before combining them, which differs meaningfully from raw
#'   alpha mainly when items vary widely in scale or format.}
#' \item{freq}{The total-score frequency distribution table returned by
#'   [freq_score()] (`score`, `freq`, `pct`, `cum_pct`).}
#' \item{call}{The matched call, as used by the `print()`/`summary()`
#'   methods for `"ctt"` objects.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [freq_score()], [print.ctt()], [summary.ctt()]
#'
#' @references
#'   Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
#'   tests. *Psychometrika, 16*(3), 297-334. \doi{10.1007/BF02310555}.
#'
#'   Osburn, H. G. (2000). Coefficient alpha and related internal consistency
#'   reliability coefficients. *Psychological Methods, 5*(3), 343-355.
#'   \doi{10.1037/1082-989X.5.3.343}.
#'
#' @examples
#' # A dichotomous example: simulate the responses of 300 examinees to 15
#' # 3PLM items, then run ctt() on the simulated data
#' set.seed(1)
#' x_bin <- shape_df(
#'   par.drm = list(a = rep(1.5, 15), b = seq(-1.5, 1.5, length.out = 15),
#'                  g = rep(0.2, 15)),
#'   cats = 2, model = "3PLM"
#' )
#' dat <- simdat(x = x_bin, theta = rnorm(300), D = 1)
#' out <- ctt(data = dat)
#' out
#' summary(out)
#'
#' # A more realistic, mixed-format example: simulate response data for a
#' # 55-item test (50 dichotomous 3PLM items + 5 polytomous GRM items) from
#' # item parameters imported from a flexMIRT output file bundled with irtQ,
#' # then run ctt() on the simulated data
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#' set.seed(2)
#' theta <- rnorm(1000)
#' dat_mixed <- simdat(x = x, theta = theta, D = 1)
#' out_mixed <- ctt(data = dat_mixed, item.id = x$id, cats = x$cats)
#' summary(out_mixed)
#'
#' @export
ctt <- function(data, item.id = NULL, cats = NULL, correct = FALSE,
                 missing = NA, flag = TRUE, crit.p = c(0.10, 0.95),
                 crit.dis = 0.20) {

  # capture the matched call for display in print()/summary() methods, as
  # done for irtQ's own est_irt()/est_item()/est_mg() objects
  call <- match.call()

  # coerce to a plain data frame so column-wise access behaves consistently
  # for matrix/tibble/data.frame input alike
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  # recode a user-specified missing-value sentinel to NA before analysis
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # listwise-delete any examinee with a remaining missing response, warning
  # once, so every statistic below uses the same set of examinees
  complete_rows <- stats::complete.cases(data)
  n_dropped <- sum(!complete_rows)
  if (n_dropped > 0L) {
    warning(n_dropped, " examinee(s) with missing item responses were ",
            "excluded listwise from ctt().", call. = FALSE)
  }
  data <- data[complete_rows, , drop = FALSE]

  # at least one examinee with complete responses is required
  if (nrow(data) == 0L) {
    stop("No examinee has complete item responses.", call. = FALSE)
  }

  # item-level statistics: difficulty, raw/corrected discrimination,
  # alpha-with-item-removed, and flagging; `data` has no missing values left
  item_out <- ctt_item(data = data, item.id = item.id, cats = cats,
                        correct = correct, missing = NA, flag = flag,
                        crit.p = crit.p, crit.dis = crit.dis)

  # test-level reliability summary (alpha, SEM, mean difficulty/discrimination)
  alpha_out <- ctt_alpha(data = data, item.id = item.id, cats = cats,
                          correct = correct, missing = NA)

  # total score of the retained examinees, which freq_score() tabulates
  total <- rowSums(data)

  # total-score frequency distribution (no missing values remain)
  freq_out <- freq_score(score = total, missing = NA)

  # bundle all three module results together under a single S3-classed object
  out <- list(item = item_out$item, crit = item_out$crit, alpha = alpha_out,
              freq = freq_out, call = call)
  class(out) <- "ctt"
  out
}
