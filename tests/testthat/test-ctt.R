# ctt_item() and ctt_alpha() are internal (not exported); access via :::,
# following the same convention as test-confirm_df.R. ctt(), print.ctt(), and
# summary.ctt() are exported and tested directly.
ctt_item <- irtQ:::ctt_item
ctt_alpha <- irtQ:::ctt_alpha

# ── ctt_item(): dichotomous ─────────────────────────────────────────────────

test_that("ctt_item() computes difficulty and item-total correlations matching a from-scratch calculation", {
  set.seed(1)
  dat <- data.frame(matrix(rbinom(50 * 6, 1, 0.6), nrow = 50))
  out <- ctt_item(data = dat, flag = FALSE)$item

  total <- rowSums(dat)
  for (j in seq_len(ncol(dat))) {
    expect_equal(out$difficulty[j], round(mean(dat[[j]]), 3))
    expect_equal(out$discrimination_raw[j], round(stats::cor(dat[[j]], total), 3))
    expect_equal(out$discrimination_corrected[j],
                 round(stats::cor(dat[[j]], total - dat[[j]]), 3))
  }
})

test_that("ctt_item() generalizes correctly to polytomous items", {
  set.seed(2)
  dat <- data.frame(matrix(sample(0:3, 40 * 4, replace = TRUE), nrow = 40))
  out <- ctt_item(data = dat, cats = rep(4, 4), flag = FALSE)$item

  expect_equal(out$difficulty, unname(round(vapply(dat, mean, numeric(1)) / 3, 3)))
  expect_true(all(out$cats == 4))
})

test_that("ctt_item() flags extreme difficulty as documented, at the correct boundary", {
  dat <- data.frame(
    easy = c(rep(1, 19), 0),   # difficulty = 19/20 = 0.95 -> boundary, NOT > 0.95
    hard = c(rep(0, 19), 1),   # difficulty = 1/20  = 0.05 -> < 0.10, should flag
    ok   = rep(c(1, 0), 10)    # difficulty = 0.50 -> within bounds
  )
  out <- ctt_item(data = dat, item.id = c("easy", "hard", "ok"),
                   crit.p = c(0.10, 0.95), crit.dis = 0.20)$item
  expect_match(out$flag[out$item == "hard"], "difficulty too low")
  expect_false(grepl("difficulty", out$flag[out$item == "easy"]))
  expect_false(grepl("difficulty", out$flag[out$item == "ok"]))
})

test_that("ctt_item() uses custom item.id and does not fall back to colnames", {
  dat <- data.frame(foo = c(1, 0, 1, 0), bar = c(1, 1, 0, 0))
  out <- ctt_item(data = dat, item.id = c("Q1", "Q2"), flag = FALSE)$item
  expect_equal(out$item, c("Q1", "Q2"))

  out_default <- ctt_item(data = dat, flag = FALSE)$item
  expect_equal(out_default$item, c("V1", "V2"))  # not "foo"/"bar"
})

test_that("ctt_item() errors when item.id length does not match ncol(data)", {
  dat <- data.frame(a = c(1, 0), b = c(0, 1))
  expect_error(ctt_item(data = dat, item.id = c("Q1", "Q2", "Q3")),
               "length\\(item.id\\) must equal ncol\\(data\\)")
})

test_that("ctt_item() excludes examinees with missing responses listwise, with a warning", {
  dat <- data.frame(I1 = c(1, 0, NA, 1), I2 = c(1, 1, 0, 0), I3 = c(0, 1, 1, 1))
  expect_warning(out <- ctt_item(data = dat, flag = FALSE), "excluded listwise")
  expect_equal(nrow(out$item), 3)  # still 3 items
})

test_that("ctt_item() recodes a custom missing sentinel before listwise deletion", {
  dat <- data.frame(I1 = c(1, 0, -9, 1), I2 = c(1, 1, 0, 0), I3 = c(0, 1, 1, 1))
  expect_warning(out <- ctt_item(data = dat, missing = -9, flag = FALSE),
                 "excluded listwise")
  expect_equal(out$item$difficulty, ctt_item(
    data = data.frame(I1 = c(1, 0, 1), I2 = c(1, 1, 0), I3 = c(0, 1, 1)),
    flag = FALSE
  )$item$difficulty)
})

test_that("ctt_item() errors with fewer than two items", {
  dat <- data.frame(I1 = c(1, 0, 1))
  expect_error(ctt_item(data = dat), "at least two items")
})

# ── ctt_alpha() ──────────────────────────────────────────────────────────────

test_that("ctt_alpha() raw alpha matches the standard variance-based formula", {
  set.seed(4)
  dat <- data.frame(matrix(rbinom(100 * 5, 1, 0.55), nrow = 100))
  out <- ctt_alpha(data = dat)

  k <- ncol(dat)
  item_var <- vapply(dat, stats::var, numeric(1))
  total_var <- stats::var(rowSums(dat))
  alpha_ref <- (k / (k - 1)) * (1 - sum(item_var) / total_var)

  expect_equal(out$alpha, round(alpha_ref, 3))
  expect_equal(out$sem, round(stats::sd(rowSums(dat)) * sqrt(1 - alpha_ref), 3))
})

test_that("ctt_alpha() standardized alpha matches the mean-inter-item-correlation formula", {
  set.seed(5)
  dat <- data.frame(matrix(rbinom(100 * 5, 1, 0.5), nrow = 100))
  out <- ctt_alpha(data = dat)

  k <- ncol(dat)
  item_cor <- stats::cor(dat)
  r_bar <- mean(item_cor[upper.tri(item_cor)])
  alpha_std_ref <- (k * r_bar) / (1 + (k - 1) * r_bar)

  expect_equal(out$alpha_std, round(alpha_std_ref, 3))
})

test_that("ctt_alpha() mean discrimination columns match ctt_item()'s own output", {
  set.seed(6)
  dat <- data.frame(matrix(rbinom(80 * 6, 1, 0.6), nrow = 80))
  item_out <- ctt_item(data = dat, flag = FALSE)$item
  alpha_out <- ctt_alpha(data = dat)

  expect_equal(alpha_out$mean_discrimination_raw,
               round(mean(item_out$discrimination_raw, na.rm = TRUE), 3))
  expect_equal(alpha_out$mean_discrimination_corrected,
               round(mean(item_out$discrimination_corrected, na.rm = TRUE), 3))
})

test_that("ctt_alpha() returns NA alpha when total score has zero variance", {
  dat <- data.frame(I1 = c(1, 1, 1), I2 = c(0, 0, 0), I3 = c(1, 1, 1))
  # total score is constant (2,2,2) for all examinees
  out <- ctt_alpha(data = dat)
  expect_true(is.na(out$alpha))
  expect_true(is.na(out$sem))
})

# ── ctt(): combined wrapper + print/summary S3 methods ──────────────────────

test_that("ctt() bundles the same item/alpha content as calling the pieces separately", {
  set.seed(7)
  dat <- data.frame(matrix(rbinom(150 * 8, 1, 0.6), nrow = 150))
  out <- ctt(data = dat)

  expect_s3_class(out, "ctt")
  expect_identical(out$item, ctt_item(data = dat)$item)
  expect_identical(out$alpha, ctt_alpha(data = dat))
})

test_that("ctt()'s internally derived total score matches freq_score()'s tabulated scores", {
  set.seed(8)
  dat <- data.frame(matrix(rbinom(60 * 5, 1, 0.5), nrow = 60))
  out <- ctt(data = dat)
  ref_freq <- freq_score(rowSums(dat))
  expect_equal(out$freq, ref_freq)
})

test_that("print.ctt() and summary.ctt()/print.summary.ctt() run without error and return invisibly", {
  set.seed(9)
  dat <- data.frame(matrix(rbinom(40 * 5, 1, 0.6), nrow = 40))
  out <- ctt(data = dat)

  expect_output(print(out), "Classical Test Theory")
  expect_true(inherits(print(out), "ctt"))

  smry <- summary(out)
  expect_s3_class(smry, "summary.ctt")
  expect_output(print(smry), "Item-Level Statistics")
  expect_output(print(smry), "Test-Level Reliability Summary")
  expect_output(print(smry), "Total-Score Frequency Distribution")
})

test_that("ctt() works on the bundled LSAT6 dataset without error", {
  skip_if_not(exists("LSAT6"), "LSAT6 dataset not available")
  out <- ctt(data = as.data.frame(LSAT6))
  expect_s3_class(out, "ctt")
  expect_true(out$alpha$alpha >= 0 && out$alpha$alpha <= 1)
})
