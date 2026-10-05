# ctt_distr() computes option/category response distributions and
# point-biserial correlations, in both selected-response (with a key) and
# scored-category (already-scored, no key) modes.

test_that("ctt_distr() selected-response mode reproduces a hand-worked distribution", {
  raw <- data.frame(
    V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
    V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
  )
  key <- c(1, 4)
  out <- ctt_distr(data = raw, key = key)

  # V1: option 1 chosen by rows 1,3,4,7 -> freq 4 of 8 -> 50%
  v1_opt1 <- out$distr[out$distr$item == "V1" & out$distr$option == 1, ]
  expect_equal(v1_opt1$freq, 4L)
  expect_equal(v1_opt1$pct, 50)
  expect_true(v1_opt1$is_key)

  # V2: option 4 is the key; row 4 ("2,4") is double-marked and must not count
  # as choosing either option 2 or option 4
  v2_opt4 <- out$distr[out$distr$item == "V2" & out$distr$option == 4, ]
  expect_equal(v2_opt4$freq, 5L)  # rows 1,2,5,7,8
  expect_false(v2_opt4$freq == 6L)

  # omission/double-marking: V2 has exactly one double-marked response
  expect_equal(out$omit$pct_double[out$omit$item == "V2"], round(100 * 1 / 8, 2))
  expect_equal(out$omit$pct_blank[out$omit$item == "V1"], 0)
})

test_that("ctt_distr() flags a distractor whose point-biserial correlation is positive", {
  # 10 examinees ranked by an externally supplied total score; the wrong
  # option "2" is chosen only by the 5 highest-scoring examinees, and the
  # keyed option "1" only by the 5 lowest-scoring examinees, so option 2's
  # point-biserial correlation with total must be strongly positive.
  raw <- data.frame(target = c("2", "2", "2", "2", "2", "1", "1", "1", "1", "1"))
  total <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
  out <- ctt_distr(data = raw, key = c(1), total = total)
  distractor_row <- out$distr[out$distr$option == 2, ]
  expect_gt(distractor_row$pb_raw, 0)
  expect_match(distractor_row$flag, "more attractive")
})

test_that("ctt_distr() scored-category mode reports category frequencies with no key concept", {
  scored <- data.frame(
    I1 = c(0, 1, 2, 1, 0, 2, 1, 1),
    I2 = c(1, 1, 0, 2, 2, 1, 0, 1)
  )
  out <- ctt_distr(data = scored)
  expect_null(out$omit)
  expect_true(all(is.na(out$distr$is_key)))

  i1_cat0 <- out$distr[out$distr$item == "V1" & out$distr$option == 0, ]
  expect_equal(i1_cat0$freq, sum(scored$I1 == 0))
})

test_that("ctt_distr() errors when a supplied total score has the wrong length", {
  raw <- data.frame(V1 = c("1", "2", "1"))
  expect_error(ctt_distr(data = raw, key = c(1), total = c(1, 2)),
               "length\\(total\\) must equal nrow\\(data\\)")
})

test_that("ctt_distr() pb_raw/pb_corrected are NA for a constant (all-key) option", {
  raw <- data.frame(V1 = c("1", "1", "1", "1"))
  out <- ctt_distr(data = raw, key = c(1))
  row <- out$distr[out$distr$item == "V1" & out$distr$option == 1, ]
  expect_true(is.na(row$pb_raw))
  expect_true(is.na(row$pb_corrected))
})

test_that("ctt_distr() scores a letter-coded (A-E) test form identically to its numeric analog", {
  # same responses/keys as the hand-worked numeric example above, but
  # expressed with letter option labels (1->A, 2->B, 3->C, 4->D); every
  # numeric output column should match row-for-row once options are mapped
  raw_num <- data.frame(
    V1 = c("1", "2", "1", "1", "3", "2", "1", "4"),
    V2 = c("4", "4", "2", "2,4", "4", "1", "4", "4")
  )
  out_num <- ctt_distr(data = raw_num, key = c(1, 4))

  raw_let <- data.frame(
    V1 = c("A", "B", "A", "A", "C", "B", "A", "D"),
    V2 = c("D", "D", "B", "B,D", "D", "A", "D", "D")
  )
  out_let <- ctt_distr(data = raw_let, key = c("A", "D"))

  option_map <- c("A" = 1, "B" = 2, "C" = 3, "D" = 4)
  out_let$distr$option_num <- option_map[out_let$distr$option]
  merged <- merge(out_num$distr, out_let$distr,
                   by.x = c("item", "option"), by.y = c("item", "option_num"),
                   suffixes = c("_num", "_let"))

  expect_equal(nrow(merged), nrow(out_num$distr))
  expect_equal(merged$freq_num, merged$freq_let)
  expect_equal(merged$pct_num, merged$pct_let)
  expect_equal(merged$pb_raw_num, merged$pb_raw_let)
  expect_equal(merged$pb_corrected_num, merged$pb_corrected_let)
  expect_equal(merged$is_key_num, merged$is_key_let)
  expect_equal(merged$flag_num, merged$flag_let)
  expect_equal(out_num$omit$pct_blank, out_let$omit$pct_blank)
  expect_equal(out_num$omit$pct_double, out_let$omit$pct_double)
})

test_that("ctt_distr() matches letter-coded options case-insensitively", {
  raw <- data.frame(V1 = c("a", "A", "b", "B", "a", "B", "A", "a"))
  out <- ctt_distr(data = raw, key = c("A"))
  row_a <- out$distr[out$distr$option == "A", ]
  expect_equal(row_a$freq, 5L)
  expect_true(row_a$is_key)
})

test_that("ctt_distr() scores a general-label (non-Latin) test form correctly", {
  # option labels are Korean syllables (romanized ga/na/da), written with
  # \uXXXX escapes for ASCII source portability, as in score_resp()
  raw <- data.frame(V1 = c("\uAC00", "\uB098", "\uB2E4", "\uAC00",
                            "\uB098", "\uB2E4", "\uAC00", "\uAC00"))
  out <- ctt_distr(data = raw, key = c("\uAC00"))
  row_key <- out$distr[out$distr$option == "\uAC00", ]
  expect_equal(row_key$freq, 4L)
  expect_true(row_key$is_key)
})

test_that("ctt_distr() respects an explicit character opt= for a letter-coded item", {
  raw <- data.frame(
    V1 = c("A", "B", "A", "A", "C", "B", "A", "D"),
    V2 = c("D", "D", "B", "B,D", "D", "A", "D", "D")
  )
  out <- ctt_distr(data = raw, key = c("A", "D"), opt = c("A", "B", "C", "D"))
  expect_equal(sort(unique(out$distr$option[out$distr$item == "V1"])),
               c("A", "B", "C", "D"))
})
