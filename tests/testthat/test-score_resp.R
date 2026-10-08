# score_resp() dichotomously scores raw selected-response data against an
# answer key. These tests check correct scoring, missing/double-mark
# handling, key validation, and the resp_summary tallies.

test_that("score_resp() scores a small hand-worked example correctly", {
  raw <- data.frame(
    V1 = c("1", "2", NA, "1", "3"),
    V2 = c("4", "4", "4", "2,4", "4"),
    V3 = c("2", "1", "3", "3", NA)
  )
  key <- c(1, 4, 3)
  out <- score_resp(data = raw, key = key)

  # hand-worked expected scores:
  # V1: 1=correct, 2=wrong, NA=blank->0, 1=correct, 3=wrong  -> 1,0,0,1,0
  # V2: 4=correct, 4=correct, 4=correct, "2,4"=double->0, 4=correct -> 1,1,1,0,1
  # V3: 2=wrong, 1=wrong, 3=correct, 3=correct, NA=blank->0 -> 0,0,1,1,0
  expect_equal(out$scored$V1, c(1L, 0L, 0L, 1L, 0L))
  expect_equal(out$scored$V2, c(1L, 1L, 1L, 0L, 1L))
  expect_equal(out$scored$V3, c(0L, 0L, 1L, 1L, 0L))

  # resp_summary tallies
  expect_equal(out$resp_summary$n_blank, c(1L, 0L, 1L))
  expect_equal(out$resp_summary$n_double, c(0L, 1L, 0L))
  expect_equal(out$resp_summary$n_correct, c(2L, 4L, 2L))
  expect_equal(out$resp_summary$n_wrong, out$resp_summary$n - out$resp_summary$n_correct)
})

test_that("score_resp() recodes a custom missing sentinel identically to NA", {
  raw <- data.frame(
    V1 = c("1", "2", NA, "1", "3"),
    V2 = c("4", "4", "4", "2,4", "4"),
    V3 = c("2", "1", "3", "3", NA)
  )
  raw2 <- data.frame(
    V1 = c("1", "2", "9", "1", "3"),
    V2 = c("4", "4", "4", "2,4", "4"),
    V3 = c("2", "1", "3", "3", "9")
  )
  key <- c(1, 4, 3)
  out1 <- score_resp(data = raw, key = key)
  out2 <- score_resp(data = raw2, key = key, missing = "9")
  expect_identical(out1$scored, out2$scored)
})

test_that("score_resp() accepts key as a data frame, sorted by item", {
  raw <- data.frame(
    V1 = c("1", "2", "1"),
    V2 = c("4", "4", "2"),
    V3 = c("2", "1", "3")
  )
  # rows deliberately out of column order
  key_df <- data.frame(item = c(3, 1, 2), key = c(3, 1, 4))
  out_vec <- score_resp(data = raw, key = c(1, 4, 3))
  out_df <- score_resp(data = raw, key = key_df)
  expect_identical(out_vec$scored, out_df$scored)
})

test_that("score_resp() errors on key/data length mismatch", {
  raw <- data.frame(V1 = c("1", "2"), V2 = c("2", "1"))
  expect_error(score_resp(data = raw, key = c(1, 2, 3)),
               "length\\(key\\) must equal ncol\\(data\\)")
})

test_that("score_resp() errors when key data frame is missing required columns", {
  raw <- data.frame(V1 = c("1", "2"), V2 = c("2", "1"))
  bad_key <- data.frame(item = c(1, 2), answer = c(1, 2))
  expect_error(score_resp(data = raw, key = bad_key),
               "must contain columns")
})

test_that("score_resp() errors when key$item has duplicates or gaps", {
  raw <- data.frame(V1 = c("1", "2"), V2 = c("2", "1"))
  bad_key <- data.frame(item = c(1, 1), key = c(1, 2))
  expect_error(score_resp(data = raw, key = bad_key),
               "no duplicates or gaps")
})

test_that("score_resp() warns on and still scores invalid response tokens", {
  raw <- data.frame(V1 = c("1", "x", "2"))
  key <- c(1)
  expect_warning(out <- score_resp(data = raw, key = key), "n_invalid")
  expect_equal(out$scored$V1, c(1L, 0L, 0L))
  expect_equal(out$resp_summary$n_invalid, 1L)
})

test_that("score_resp() errors on empty data", {
  expect_error(score_resp(data = data.frame(), key = numeric(0)),
               "at least one row and one column")
})

test_that("score_resp() falls back to V1, V2, ... when data has no column names", {
  raw <- matrix(c("1", "2", "1", "4"), nrow = 2, dimnames = NULL)
  colnames(raw) <- NULL
  out <- score_resp(data = as.data.frame(raw), key = c(1, 4))
  expect_true(all(colnames(out$scored) %in% c("V1", "V2")))
})

test_that("score_resp() scores a letter-coded item identically to its numeric analog", {
  # same correctness pattern as the numeric hand-worked example above, but
  # expressed with letter option labels and a lowercase response ("a")
  raw3 <- data.frame(
    V1 = c("A", "B", NA, "a", "C"),
    V2 = c("D", "D", "D", "B,D", "D"),
    V3 = c("B", "A", "C", "C", NA)
  )
  key3 <- c("A", "D", "C")
  out3 <- score_resp(data = raw3, key = key3)

  expect_equal(out3$scored$V1, c(1L, 0L, 0L, 1L, 0L))
  expect_equal(out3$scored$V2, c(1L, 1L, 1L, 0L, 1L))
  expect_equal(out3$scored$V3, c(0L, 0L, 1L, 1L, 0L))
  expect_equal(out3$resp_summary$n_blank, c(1L, 0L, 1L))
  expect_equal(out3$resp_summary$n_double, c(0L, 1L, 0L))
})

test_that("score_resp() matches letter responses case-insensitively", {
  raw <- data.frame(V1 = c("a", "A", "b", "B"))
  out_lower_key <- score_resp(data = raw, key = c("a"))
  out_upper_key <- score_resp(data = raw, key = c("A"))
  expect_identical(out_lower_key$scored, out_upper_key$scored)
  expect_equal(out_upper_key$scored$V1, c(1L, 1L, 0L, 0L))
})

test_that("score_resp() supports a mixed test form (numeric- and letter-coded items)", {
  raw4 <- data.frame(
    V1 = c("1", "2", "1", "1", "3"),
    V2 = c("D", "D", "D", "B,D", "D"),
    V3 = c("B", "A", "C", "C", NA)
  )
  key4 <- c(1, "D", "C")
  out4 <- score_resp(data = raw4, key = key4)

  expect_equal(out4$scored$V1, c(1L, 0L, 1L, 1L, 0L))
  expect_equal(out4$scored$V2, c(1L, 1L, 1L, 0L, 1L))
  expect_equal(out4$scored$V3, c(0L, 0L, 1L, 1L, 0L))
})

test_that("score_resp() flags a numeric response on a letter-coded item as invalid", {
  raw <- data.frame(V1 = c("A", "1", "B"))
  key <- c("A")
  expect_warning(out <- score_resp(data = raw, key = key), "n_invalid")
  expect_equal(out$scored$V1, c(1L, 0L, 0L))
  expect_equal(out$resp_summary$n_invalid, 1L)
})

test_that("score_resp() treats a non-Latin, non-numeric key as a general label", {
  # a mixed alphanumeric key value like "1A" is accepted as a general option
  # label (matched by exact, case-insensitive string comparison), since
  # score_resp() does not restrict option labels to plain numbers or
  # single-script Latin letters
  raw <- data.frame(V1 = c("1", "2"), V2 = c("A", "B"))
  out <- score_resp(data = raw, key = c("1A", "B"))
  # V1's key "1A" is a general label ("1" and "2" both fail to match it
  # exactly); V2's key "B" is Latin-letter-coded ("A" != "B",
  # "B" == "B")
  expect_equal(out$scored$V1, c(0L, 0L))
  expect_equal(out$scored$V2, c(0L, 1L))
})

test_that("score_resp() errors only when a key value is blank or NA", {
  raw <- data.frame(V1 = c("1", "2"), V2 = c("A", "B"))
  expect_error(score_resp(data = raw, key = c(NA, "B")),
               "missing or blank value")
  expect_error(score_resp(data = raw, key = c("", "B")),
               "missing or blank value")
})

test_that("score_resp() scores a Korean-labeled (general label scheme) item correctly", {
  raw6 <- data.frame(V1 = c("\uAC00", "\uB098", "\uB2E4", "\uB2E4", NA)) # ga/na/da/da/NA
  key6 <- c("\uB2E4") # da
  out6 <- score_resp(data = raw6, key = key6)
  expect_equal(out6$scored$V1, c(0L, 0L, 1L, 1L, 0L))
  expect_equal(out6$resp_summary$n_blank, 1L)
  expect_equal(out6$resp_summary$n_invalid, 0L)
})

test_that("score_resp() never flags n_invalid for a general label-coded item", {
  # under the general label scheme, any non-blank, non-double-marked token
  # (even one that looks nothing like the key, or that would have been
  # "invalid" under the numeric/Latin-letter schemes) is scored 0 as simply
  # wrong, not tallied as n_invalid, since there is no universal format
  # check for an arbitrary label alphabet
  raw <- data.frame(V1 = c("\uAC00", "1", "xyz", "!!")) # ga, 1, xyz, !!
  out <- score_resp(data = raw, key = c("\uAC00")) # ga
  expect_equal(out$scored$V1, c(1L, 0L, 0L, 0L))
  expect_equal(out$resp_summary$n_invalid, 0L)
})

test_that("score_resp() scores a full five-option (A-E) letter item, case-insensitively", {
  raw5 <- data.frame(V1 = c("A", "B", "C", "D", "E", "c"))
  out5 <- score_resp(data = raw5, key = c("C"))
  expect_equal(out5$scored$V1, c(0L, 0L, 1L, 0L, 0L, 1L))
})

test_that("score_resp() matches a key data frame by item number for any item column type", {
  # twelve items so that string order differs from numeric order
  set.seed(11)
  n_item <- 12
  raw <- as.data.frame(
    matrix(sample(as.character(1:4), 30 * n_item, replace = TRUE), 30, n_item)
  )
  key_vec <- rep(c(1, 4, 3, 2), length.out = n_item)
  out_vec <- score_resp(data = raw, key = key_vec)

  # the same key rows in a shuffled order, with a numeric, character, and
  # factor item column
  shuffled <- c(5, 12, 1, 9, 3, 11, 7, 2, 10, 6, 8, 4)
  key_num <- data.frame(item = shuffled, key = key_vec[shuffled])
  key_chr <- data.frame(item = as.character(shuffled), key = key_vec[shuffled])
  key_fct <- data.frame(item = factor(shuffled), key = key_vec[shuffled])
  expect_identical(score_resp(data = raw, key = key_num)$scored, out_vec$scored)
  expect_identical(score_resp(data = raw, key = key_chr)$scored, out_vec$scored)
  expect_identical(score_resp(data = raw, key = key_fct)$scored, out_vec$scored)
})

test_that("score_resp() invalid-response warning does not refer to option numbers", {
  # a letter-coded item with an unrecognized token
  raw <- data.frame(V1 = c("A", "B", "7", "A"))
  expect_warning(out <- score_resp(data = raw, key = "A"),
                 "neither a valid option for this item")
  expect_equal(out$resp_summary$n_invalid, 1L)
})
