# ── Common arguments ─────────────────────────────────────────────────────────
# est_item() is a one-shot FAPC routine (no EM iteration), so the
# control set is small.  use.gprior = FALSE keeps the 3PLM tests
# from drifting on the small samples used here, and verbose = FALSE
# silences the "Starting..." / "Parsing input..." messages.
ITEM_ARGS <- list(D = 1, use.aprior = FALSE, use.gprior = FALSE,
                  verbose = FALSE)

# ── Fixtures ─────────────────────────────────────────────────────────────────
prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_full   <- bring.flexmirt(file = prm_file, "par")$Group1$full_df  # 55 items

# 2PLM-only bank for single-model coverage
x_2plm5 <- shape_df(
  par.drm = list(a = c(1.0, 0.9, 1.2, 1.1, 0.95),
                 b = c(-1.5, -0.5, 0.0, 0.5, 1.5),
                 g = rep(0, 5)),
  cats = rep(2L, 5L), model = "2PLM"
)

# All-1PLM bank for the constrained-a code path (fix.a.1pl = FALSE)
x_1plm5 <- shape_df(
  par.drm = list(a = rep(1.0, 5),
                 b = c(-2.0, -1.0, 0.0, 1.0, 2.0),
                 g = rep(0, 5)),
  cats = rep(2L, 5L), model = "1PLM"
)

# Mixed 1PLM-constrained + 2PLM + 3PLM bank.  This is THE scenario
# that exercises the bug the recent fix targeted: estimation order
# is [loc_1p_const items first, then loc_else items], so c(1,3,5,
# 2,4,6) -- different from natural order -- and the (id, params)
# pairing has to be permuted back via order(c(loc_1p_const,
# loc_else)) before the downstream cbind to x[, 1:3].
x_mixed_drm6 <- shape_df(
  par.drm = list(
    a = c(1.0, 1.0, 1.0, 1.2, 1.0, 0.9),    # 1PLM share a; 2PLM/3PLM own
    b = c(-1.5, -0.5, 0.0, 0.5, 1.0, 1.5),
    g = c(0,    0,    0,    0,   0,   0.18)
  ),
  cats = rep(2L, 6L),
  model = c("1PLM", "2PLM", "1PLM", "2PLM", "1PLM", "3PLM")
)

# Mixed-format bank (DRM + GRM) for the polytomous code path.
# Both GRM items use cats = 4 to exercise the same threshold-count
# branch (the small-cats variants of hess_item_prm have separate
# numerical fragilities that are out of scope for this test file).
x_mixed_grm <- shape_df(
  par.drm = list(a = c(1.0, 1.1), b = c(-0.4, 0.4), g = c(0, 0)),
  par.prm = list(a = c(1.0, 1.1),
                 d = list(c(-0.6, 0.2, 0.7), c(-0.5, 0.3, 0.8))),
  item.id = c("D1", "D2", "P1", "P2"),
  cats  = c(2L, 2L, 4L, 4L),
  model = c("2PLM", "2PLM", "GRM", "GRM")
)

# Helper to simulate response data + true ability scores from a bank.
# Returning the same theta as score mimics the FAPC use case where
# examinee abilities are already known (e.g. from a calibrated form).
sim_with_score <- function(x, n = 1500L, seed = 1L) {
  set.seed(seed)
  theta <- rnorm(n)
  list(
    data  = simdat(x = x, theta = theta, D = 1),
    score = theta
  )
}

run_item <- function(x, sim, ...) {
  do.call(est_item,
          c(list(x = x, data = sim$data, score = sim$score, ...), ITEM_ARGS))
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Output class and structure
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_item() returns class 'est_item'", {
  sim <- sim_with_score(x_2plm5, seed = 11)
  fit <- run_item(x_2plm5, sim, fix.a.1pl = TRUE)
  expect_s3_class(fit, "est_item")
})

test_that("est_item() top-level slots are present", {
  sim <- sim_with_score(x_2plm5, seed = 12)
  fit <- run_item(x_2plm5, sim, fix.a.1pl = TRUE)
  expected <- c("estimates", "par.est", "se.est", "pos.par", "covariance",
                "loglikelihood", "group.par", "data", "score", "scale.D",
                "convergence", "nitem", "deleted.item", "npar.est",
                "n.response", "TotalTime", "call")
  expect_true(all(expected %in% names(fit)))
})

test_that("est_item() par.est / se.est are data.frames with correct nrow", {
  sim <- sim_with_score(x_mixed_drm6, seed = 13)
  fit <- run_item(x_mixed_drm6, sim, fix.a.1pl = FALSE)
  expect_s3_class(fit$par.est, "data.frame")
  expect_s3_class(fit$se.est,  "data.frame")
  expect_equal(nrow(fit$par.est), nrow(x_mixed_drm6))
  expect_equal(nrow(fit$se.est),  nrow(x_mixed_drm6))
})

test_that("est_item() loglikelihood is a finite scalar", {
  sim <- sim_with_score(x_2plm5, seed = 14)
  fit <- run_item(x_2plm5, sim, fix.a.1pl = TRUE)
  expect_true(is.finite(fit$loglikelihood))
})

test_that("est_item() covariance is a square npar x npar matrix", {
  sim <- sim_with_score(x_mixed_drm6, seed = 15)
  fit <- run_item(x_mixed_drm6, sim, fix.a.1pl = FALSE)
  expect_true(is.matrix(fit$covariance))
  expect_equal(nrow(fit$covariance), ncol(fit$covariance))
  expect_equal(nrow(fit$covariance), fit$npar.est)
})


# ══════════════════════════════════════════════════════════════════════════════
# 2. Bug regression: (id, parameter) row alignment for mixed banks
# ══════════════════════════════════════════════════════════════════════════════
# Pre-fix, the per-item estimation loop in est_item.R appended results
# in [loc_1p_const items first, then loc_else items] order, but the
# downstream pipeline used dplyr::arrange("loc") -- which sorts by the
# literal string "loc" (a no-op) -- so par_df / se_df rows stayed in
# estimation order while the cbind to x[, 1:3] used natural order.
# Result: the row labelled with id "M2" carried params for some other
# item.  The shared-a check below is the discriminating signal,
# because the 1PLM-constrained code path forces all 1PLM items to
# share a single estimated `a` -- if mis-pairing brings non-1PLM
# values into 1PLM rows, the rows no longer agree.

test_that("est_item() mixed 1PLM-const + others: par.est$id stays in input order", {
  sim <- sim_with_score(x_mixed_drm6, seed = 21)
  fit <- run_item(x_mixed_drm6, sim, fix.a.1pl = FALSE)
  expect_identical(as.character(fit$par.est$id), as.character(x_mixed_drm6$id))
  expect_identical(as.character(fit$par.est$model), as.character(x_mixed_drm6$model))
})

test_that("est_item() mixed 1PLM-const + others: 1PLM rows share a single `a`", {
  sim <- sim_with_score(x_mixed_drm6, seed = 22)
  fit <- run_item(x_mixed_drm6, sim, fix.a.1pl = FALSE)
  is_1plm <- as.character(fit$par.est$model) == "1PLM"
  expect_true(any(is_1plm))
  # constrained-a model: every 1PLM row must carry the SAME estimated a;
  # if rows are mis-paired with non-1PLM items, the values diverge
  a_vals <- fit$par.est$par.1[is_1plm]
  expect_false(any(is.na(a_vals)))
  expect_lt(max(a_vals) - min(a_vals), 1e-10)
})

test_that("est_item() mixed: se.est$id matches x order (parallel to par.est)", {
  sim <- sim_with_score(x_mixed_drm6, seed = 23)
  fit <- run_item(x_mixed_drm6, sim, fix.a.1pl = FALSE)
  expect_identical(as.character(fit$se.est$id), as.character(x_mixed_drm6$id))
})


# ══════════════════════════════════════════════════════════════════════════════
# 3. Single-model paths
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_item() all-2PLM: par.est rows align with input order", {
  sim <- sim_with_score(x_2plm5, seed = 31)
  fit <- run_item(x_2plm5, sim, fix.a.1pl = TRUE)
  expect_identical(as.character(fit$par.est$id), as.character(x_2plm5$id))
  expect_true(all(as.character(fit$par.est$model) == "2PLM"))
})

test_that("est_item() all-1PLM with fix.a.1pl = FALSE: shared `a` across rows", {
  sim <- sim_with_score(x_1plm5, seed = 32)
  fit <- run_item(x_1plm5, sim, fix.a.1pl = FALSE)
  expect_identical(as.character(fit$par.est$id), as.character(x_1plm5$id))
  a_vals <- fit$par.est$par.1
  expect_lt(max(a_vals) - min(a_vals), 1e-10)
})

test_that("est_item() all-1PLM with fix.a.1pl = TRUE: par.1 = a.val.1pl per row", {
  sim <- sim_with_score(x_1plm5, seed = 33)
  fit <- run_item(x_1plm5, sim, fix.a.1pl = TRUE, a.val.1pl = 1)
  # when slope is fixed, par.1 takes the user-supplied a.val.1pl
  expect_equal(unname(fit$par.est$par.1), rep(1, nrow(x_1plm5)),
               tolerance = 1e-10)
})


# ══════════════════════════════════════════════════════════════════════════════
# 4. Mixed format with polytomous items
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_item() mixed DRM + GRM: par.est preserves model column per item", {
  sim <- sim_with_score(x_mixed_grm, seed = 41)
  fit <- run_item(x_mixed_grm, sim, fix.a.1pl = TRUE)
  expect_identical(as.character(fit$par.est$id),    as.character(x_mixed_grm$id))
  expect_identical(as.character(fit$par.est$model), as.character(x_mixed_grm$model))
  # GRM items must carry one extra threshold column compared to 2PLM
  grm_rows <- which(fit$par.est$model == "GRM")
  expect_true(all(!is.na(fit$par.est$par.3[grm_rows])))
})

test_that("est_item() mixed DRM + GRM: covariance is finite + symmetric", {
  sim <- sim_with_score(x_mixed_grm, seed = 42)
  fit <- run_item(x_mixed_grm, sim, fix.a.1pl = TRUE)
  expect_true(all(is.finite(fit$covariance)))
  # symmetry within numerical tolerance of the bdiag + reorder pipeline
  expect_lt(max(abs(fit$covariance - t(fit$covariance))), 1e-10)
})

test_that("est_item() handles cats = 3 GRM without diag<- crash", {
  # cats == 3 GRM was a long-standing bug in hess_item_prm_inner():
  # the off-diagonal b-block fill via diag(hess[2:m, 3:(m+1)]) <-
  # hess_b1b2 silently dropped the 1x1 sub-matrix to a scalar when
  # m == 2 (cats == 3), after which diag<-() raised
  #   "only matrix diagonals can be replaced".
  # Fixed by switching to 2-column index-matrix assignment, which
  # treats the m == 2 single-cell case the same as m >= 3 multi-
  # cell cases.  This test guards against regression.
  set.seed(701)
  x_grm3 <- shape_df(
    par.prm = list(a = c(1.0, 1.1),
                   d = list(c(-0.6, 0.5), c(-0.5, 0.4))),
    item.id = c("P1", "P2"),
    cats  = c(3L, 3L),
    model = c("GRM", "GRM")
  )
  sim <- sim_with_score(x_grm3, n = 1500L, seed = 701)
  fit <- run_item(x_grm3, sim, fix.a.1pl = TRUE)
  expect_s3_class(fit, "est_item")
  # the SE pipeline only succeeds if hess_item_prm() returned a
  # full hessian for both cats=3 GRM rows; non-finite entries
  # would propagate into se.est$par.* values
  expect_true(all(is.finite(fit$covariance)))
  expect_lt(max(abs(fit$covariance - t(fit$covariance))), 1e-10)
  # cats=3 GRM has 1 slope + 2 thresholds = 3 params per item;
  # 2 items -> 6 params total in the per-item covariance block-
  # diagonal, then reordered by reloc.par
  expect_equal(nrow(fit$covariance), 6L)
})


# ══════════════════════════════════════════════════════════════════════════════
# 5. Parameter recovery (loose tolerance)
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_item() recovers 2PLM b parameters within loose tolerance", {
  set.seed(51)
  true_b <- c(-1.5, -0.5, 0.0, 0.5, 1.5)
  x_pop  <- shape_df(par.drm = list(a = rep(1.0, 5), b = true_b,
                                     g = rep(0, 5)),
                     cats = rep(2L, 5L), model = "2PLM")
  theta <- rnorm(3000)
  data  <- simdat(x = x_pop, theta = theta, D = 1)
  fit   <- do.call(est_item,
                   c(list(x = x_pop, data = data, score = theta,
                          fix.a.1pl = TRUE), ITEM_ARGS))
  expect_equal(unname(fit$par.est$par.2), true_b, tolerance = 0.3)
})


# ══════════════════════════════════════════════════════════════════════════════
# 6. summary() and print() do not error
# ══════════════════════════════════════════════════════════════════════════════

test_that("summary.est_item() runs without error", {
  sim <- sim_with_score(x_2plm5, seed = 61)
  fit <- run_item(x_2plm5, sim, fix.a.1pl = TRUE)
  expect_no_error(capture.output(summary(fit)))
})

test_that("print.est_item() runs without error", {
  sim <- sim_with_score(x_2plm5, seed = 62)
  fit <- run_item(x_2plm5, sim, fix.a.1pl = TRUE)
  expect_no_error(capture.output(print(fit)))
})
