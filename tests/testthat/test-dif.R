# ── Shared fixtures ───────────────────────────────────────────────────────────

# Import 10 × 3PLM operational item parameters (flexMIRT sample)
prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_dif    <- bring.flexmirt(file = prm_file, "par")$Group1$full_df[1:10, ]

# Simulate 600 examinees (300 reference, 300 focal)
set.seed(55)
n_each    <- 300L
theta_ref <- rnorm(n_each, mean =  0.0, sd = 1)
theta_foc <- rnorm(n_each, mean = -0.3, sd = 1)   # mild impact
theta_all <- c(theta_ref, theta_foc)
group_vec <- c(rep(0L, n_each), rep(1L, n_each))   # 0 = reference, 1 = focal

resp_all  <- simdat(x = x_dif, theta = theta_all, D = 1)

# Pre-computed ML ability estimates (used as input to avoid re-scoring inside each test)
score_all <- est_score(
  x = x_dif, data = resp_all, D = 1,
  method = "ML", range = c(-6, 6)
)

# ── rdif() ────────────────────────────────────────────────────────────────────

test_that("rdif() returns a list with no_purify component", {
  res <- rdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1,
    alpha      = 0.05
  )
  expect_type(res, "list")
  expect_true("no_purify" %in% names(res))
})

test_that("rdif() dif_stat has one row per item", {
  res <- rdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  expect_equal(nrow(res$no_purify$dif_stat), 10L)
})

test_that("rdif() dif_stat contains RDIF_R, RDIF_S, RDIF_RS columns", {
  res <- rdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  cols <- colnames(res$no_purify$dif_stat)
  expect_true(any(grepl("rdifr",  tolower(cols))))
  expect_true(any(grepl("rdifs",  tolower(cols))))
  expect_true(any(grepl("rdifrs", tolower(cols))))
})

test_that("rdif() purify=FALSE sets purify slot to FALSE", {
  res <- rdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1,
    purify     = FALSE
  )
  expect_false(res$purify)
})

test_that("rdif() internal scoring (score=NULL) works without error", {
  expect_no_error(
    rdif(
      x          = x_dif,
      data       = resp_all,
      group      = group_vec,
      focal.name = 1L,
      D          = 1,
      method     = "ML",
      range      = c(-6, 6),
      verbose    = FALSE
    )
  )
})

test_that("rdif() p-values are in [0, 1]", {
  res <- rdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  pvals <- res$no_purify$dif_stat[, grep("^p\\.|pvalue|p_", colnames(res$no_purify$dif_stat),
                                          ignore.case = TRUE)]
  if (length(pvals) > 0) {
    expect_true(all(unlist(pvals) >= 0 & unlist(pvals) <= 1, na.rm = TRUE))
  }
})


# ── crdif() ───────────────────────────────────────────────────────────────────

test_that("crdif() returns a list with no_purify component", {
  res <- crdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  expect_type(res, "list")
  expect_true("no_purify" %in% names(res))
})

test_that("crdif() dif_stat has one row per item", {
  res <- crdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  expect_equal(nrow(res$no_purify$dif_stat), 10L)
})

test_that("crdif() dif_stat contains CRDIF columns", {
  res <- crdif(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  cols <- tolower(colnames(res$no_purify$dif_stat))
  expect_true(any(grepl("crdif", cols)))
})

test_that("crdif() and rdif() return identical results for dichotomous items", {
  res_crdif <- crdif(
    x = x_dif, data = resp_all, score = score_all$est.theta,
    group = group_vec, focal.name = 1L, D = 1
  )
  res_rdif <- rdif(
    x = x_dif, data = resp_all, score = score_all$est.theta,
    group = group_vec, focal.name = 1L, D = 1
  )
  # RDIF_RS and CRDIF_RS should be numerically equal for binary items
  rdifrs_col  <- grep("rdifrs",  tolower(colnames(res_rdif$no_purify$dif_stat)),  value = TRUE)[1]
  crdifrs_col <- grep("crdifrs", tolower(colnames(res_crdif$no_purify$dif_stat)), value = TRUE)[1]
  if (!is.na(rdifrs_col) && !is.na(crdifrs_col)) {
    expect_equal(
      res_rdif$no_purify$dif_stat[[rdifrs_col]],
      res_crdif$no_purify$dif_stat[[crdifrs_col]],
      tolerance = 1e-6
    )
  }
})

test_that("crdif() internal scoring (score=NULL) works without error", {
  expect_no_error(
    crdif(
      x          = x_dif,
      data       = resp_all,
      group      = group_vec,
      focal.name = 1L,
      D          = 1,
      method     = "ML",
      range      = c(-6, 6),
      verbose    = FALSE
    )
  )
})


# ── catsib() ─────────────────────────────────────────────────────────────────

test_that("catsib() returns a list with no_purify component", {
  res <- catsib(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  expect_type(res, "list")
  expect_true("no_purify" %in% names(res))
})

test_that("catsib() dif_stat has one row per item", {
  res <- catsib(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  expect_equal(nrow(res$no_purify$dif_stat), 10L)
})

test_that("catsib() dif_stat contains beta and p-value columns", {
  res <- catsib(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1
  )
  cols <- tolower(colnames(res$no_purify$dif_stat))
  expect_true(any(grepl("beta", cols)))
  expect_true(any(grepl("^p\\b|pvalue|p_value", cols)))
})

test_that("catsib() alpha slot matches supplied alpha", {
  res <- catsib(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1,
    alpha      = 0.01
  )
  expect_equal(res$alpha, 0.01)
})

test_that("catsib() purify=FALSE slot is FALSE", {
  res <- catsib(
    x          = x_dif,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1,
    purify     = FALSE
  )
  expect_false(res$purify)
})

test_that("catsib() internal scoring (score=NULL) works without error", {
  expect_no_error(
    catsib(
      x          = x_dif,
      data       = resp_all,
      group      = group_vec,
      focal.name = 1L,
      D          = 1,
      method     = "ML",
      range      = c(-6, 6)
    )
  )
})

test_that("catsib() retained bins satisfy min.binsize threshold for both groups", {
  res3 <- suppressWarnings(catsib(
    x          = NULL,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1,
    min.binsize = 3
  ))
  res5 <- suppressWarnings(catsib(
    x          = NULL,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1,
    min.binsize = 5
  ))
  # Every retained bin in the contingency table must satisfy the min.binsize
  # threshold for BOTH groups.  Note: n.total is NOT monotone in min.binsize
  # because the bin-count selection loop also changes num.bin, which can yield
  # fewer but larger bins that include more examinees in total.
  check_bins <- function(contingency, thresh) {
    vapply(contingency, function(ct) {
      bins <- head(ct, -1)   # drop adorn_totals summary row
      if (nrow(bins) == 0L) return(TRUE)
      all(bins$n.ref >= thresh & bins$n.foc >= thresh)
    }, logical(1L))
  }
  expect_true(all(check_bins(res3$no_purify$contingency, 3L)))
  expect_true(all(check_bins(res5$no_purify$contingency, 5L)))
})

test_that("catsib() min.binsize is enforced in the final bin filter, not just the bin-count loop", {
  res5 <- catsib(
    x          = NULL,
    data       = resp_all,
    score      = score_all$est.theta,
    se         = score_all$se.theta,
    group      = group_vec,
    focal.name = 1L,
    D          = 1,
    min.binsize = 5
  )
  # Every retained per-bin row across all contingency tables must satisfy
  # n.ref >= 5 AND n.foc >= 5.  The last row of each contingency table is the
  # adorn_totals summary row and is excluded via head(-1).
  all_ok <- vapply(res5$no_purify$contingency, function(ct) {
    bins <- head(ct, -1)
    if (nrow(bins) == 0L) return(TRUE)
    all(bins$n.ref >= 5L & bins$n.foc >= 5L)
  }, logical(1L))
  expect_true(all(all_ok))
})
