# freq_score() tabulates a frequency distribution table (freq, pct, cum_pct)
# for a vector of integer total scores.

test_that("freq_score() produces a table spanning min to max with correct totals", {
  set.seed(1)
  score <- rbinom(200, size = 10, prob = 0.5)
  out <- freq_score(score)

  expect_equal(out$score, seq.int(min(score), max(score)))
  expect_equal(sum(out$freq), length(score))
  # cumulative percentage must reach (approximately) 100 at the final row
  expect_equal(out$cum_pct[nrow(out)], 100, tolerance = 0.01)
  # pct must independently reconstruct to freq / n * 100
  expect_equal(out$pct, round(100 * out$freq / length(score), 2))
})

test_that("freq_score() includes zero-frequency scores within the observed range", {
  # scores 0, 1, 3 observed but 2 is never observed -> must still appear with freq 0
  score <- c(0, 0, 1, 3, 3, 3)
  out <- freq_score(score)
  expect_equal(out$score, 0:3)
  expect_equal(out$freq, c(2L, 1L, 0L, 3L))
})

test_that("freq_score() recodes a custom missing sentinel and drops it with a warning", {
  score <- c(1, 2, -9, 3, 2)
  expect_warning(out <- freq_score(score, missing = -9), "missing score")
  expect_equal(sum(out$freq), 4)
})

test_that("freq_score() warns and drops NA scores under the default missing = NA", {
  score <- c(1, 2, NA, 3)
  expect_warning(out <- freq_score(score), "missing score")
  expect_equal(sum(out$freq), 3)
})

test_that("freq_score() errors when every score is missing", {
  expect_error(suppressWarnings(freq_score(c(NA, NA))),
               "no non-missing values")
})

test_that("freq_score() errors on non-integer scores", {
  expect_error(freq_score(c(1, 2.5, 3)),
               "integer-valued scores")
})
