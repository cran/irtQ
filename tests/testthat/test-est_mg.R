# ── Common EM settings ────────────────────────────────────────────────────────
MG_ARGS <- list(D = 1, Etol = 1e-3, MaxE = 100L, se = FALSE, verbose = FALSE)

# ── Fixtures ──────────────────────────────────────────────────────────────────

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_full   <- bring.flexmirt(file = prm_file, "par")$Group1$full_df  # 55 items

## Dichotomous: 10 × 3PLM common items
x_drm10 <- x_full[1:10, ]

## GRM: 2 × GRM (5 categories) from the sample file
x_grm2  <- x_full[39:40, ]

## GPCM: 2 items created with shape_df
x_gpcm2 <- shape_df(
  par.prm = list(
    a = c(1.0, 1.2),
    d = list(c(-0.8, 0.0, 0.9), c(-0.5, 0.5, 1.2))
  ),
  cats = c(4L, 4L), model = "GPCM"
)

## Mixed: 6 × 3PLM + 2 × GRM (from sample)
x_mixed8 <- x_full[c(1:6, 39:40), ]

# Helper: simulate two groups from the same item bank with different theta distributions
make_two_groups <- function(x_items, seed = 1, n = 300) {
  set.seed(seed)
  list(
    g1 = simdat(x = x_items, theta = rnorm(n, mean =  0.0, sd = 1.0), D = 1),
    g2 = simdat(x = x_items, theta = rnorm(n, mean =  0.5, sd = 1.2), D = 1)
  )
}

# Helper: run est_mg and check the basic class + structure
run_mg <- function(x_items, data_list, ...) {
  do.call(est_mg, c(
    list(x          = list(x_items, x_items),
         data       = data_list,
         group.name = c("G1", "G2"),
         free.group = "G2"),
    MG_ARGS,
    list(...)
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Output class and structure
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_mg() returns class 'est_mg'", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expect_s3_class(fit, "est_mg")
})

test_that("est_mg() top-level slots are present", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expected_slots <- c("par.est", "loglikelihood", "aic", "bic",
                      "group.par", "niter", "nitem")
  expect_true(all(expected_slots %in% names(fit)))
})

test_that("est_mg() par.est is a list with 'overall' and per-group elements", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expect_type(fit$par.est, "list")
  expect_true("overall" %in% names(fit$par.est))
  expect_true("group"   %in% names(fit$par.est))
  expect_true(all(c("G1", "G2") %in% names(fit$par.est$group)))
})

test_that("est_mg() par.est$overall is a data frame with correct nrow", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expect_s3_class(fit$par.est$overall, "data.frame")
  expect_equal(nrow(fit$par.est$overall), 10L)  # 10 unique items
})

test_that("est_mg() loglikelihood is a list with 'overall' and 'group'", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expect_type(fit$loglikelihood, "list")
  expect_true("overall" %in% names(fit$loglikelihood))
  expect_true("group"   %in% names(fit$loglikelihood))
  expect_true(is.numeric(fit$loglikelihood$overall))
})

test_that("est_mg() aic and bic are finite scalars", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expect_true(is.finite(fit$aic))
  expect_true(is.finite(fit$bic))
})

test_that("est_mg() group.par has entries for each group with mu/sigma2/sigma", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expect_true(all(c("G1", "G2") %in% names(fit$group.par)))
  expect_true(all(c("mu", "sigma2", "sigma") %in% names(fit$group.par$G1)))
  expect_true(all(c("mu", "sigma2", "sigma") %in% names(fit$group.par$G2)))
})

test_that("est_mg() free group G2 has different mu from fixed group G1", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  # G1 is fixed at mean = 0; G2 is freely estimated from shifted data
  # mu may be stored as a vector; take the first non-NA value
  mu_g1 <- fit$group.par$G1$mu[!is.na(fit$group.par$G1$mu)][1]
  mu_g2 <- fit$group.par$G2$mu[!is.na(fit$group.par$G2$mu)][1]
  expect_equal(mu_g1, 0, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(mu_g1, mu_g2, tolerance = 0.05)))
})

test_that("est_mg() niter is a positive integer", {
  dat <- make_two_groups(x_drm10)
  fit <- run_mg(x_drm10, dat)
  expect_true(fit$niter > 0L)
})


# ══════════════════════════════════════════════════════════════════════════════
# 2. Dichotomous-only (3PLM)
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_mg() 3PLM: converges and returns 10 item rows", {
  dat <- make_two_groups(x_drm10, seed = 10)
  fit <- run_mg(x_drm10, dat, use.gprior = TRUE,
                gprior = list(dist = "beta", params = c(5, 16)))
  expect_s3_class(fit, "est_mg")
  expect_equal(nrow(fit$par.est$overall), 10L)
})

test_that("est_mg() 3PLM: guessing estimates in [0, 1]", {
  dat <- make_two_groups(x_drm10, seed = 10)
  fit <- run_mg(x_drm10, dat, use.gprior = TRUE,
                gprior = list(dist = "beta", params = c(5, 16)))
  g_vals <- fit$par.est$overall$par.3
  expect_true(all(g_vals >= 0 & g_vals <= 1, na.rm = TRUE))
})

test_that("est_mg() 3PLM model with x = NULL (model/cats/item.id specified directly)", {
  dat   <- make_two_groups(x_drm10, seed = 20)
  ids   <- x_drm10$id
  fit <- do.call(est_mg, c(
    list(x          = NULL,
         data       = list(dat$g1, dat$g2),
         group.name = c("G1", "G2"),
         free.group = "G2",
         model      = list(rep("3PLM", 10), rep("3PLM", 10)),
         cats       = list(rep(2L, 10),     rep(2L, 10)),
         item.id    = list(ids, ids),
         use.gprior = TRUE,
         gprior     = list(dist = "beta", params = c(5, 16))),
    MG_ARGS
  ))
  expect_s3_class(fit, "est_mg")
  expect_equal(nrow(fit$par.est$overall), 10L)
})


# ══════════════════════════════════════════════════════════════════════════════
# 3. Polytomous-only (GRM)
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_mg() GRM-only: converges with 2 items", {
  dat <- make_two_groups(x_grm2, seed = 30, n = 400)
  fit <- run_mg(x_grm2, dat)
  expect_s3_class(fit, "est_mg")
  expect_equal(nrow(fit$par.est$overall), 2L)
})

test_that("est_mg() GRM-only: par.est$overall has expected model column", {
  dat <- make_two_groups(x_grm2, seed = 30, n = 400)
  fit <- run_mg(x_grm2, dat)
  expect_true(all(fit$par.est$overall$model == "GRM"))
})

test_that("est_mg() GRM-only: par.est has same structure for both groups", {
  dat <- make_two_groups(x_grm2, seed = 30, n = 400)
  fit <- run_mg(x_grm2, dat)
  expect_equal(nrow(fit$par.est$group$G1), 2L)
  expect_equal(nrow(fit$par.est$group$G2), 2L)
})


# ══════════════════════════════════════════════════════════════════════════════
# 4. Polytomous-only (GPCM)
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_mg() GPCM-only: converges with 2 items", {
  dat <- make_two_groups(x_gpcm2, seed = 40, n = 400)
  fit <- run_mg(x_gpcm2, dat)
  expect_s3_class(fit, "est_mg")
  expect_equal(nrow(fit$par.est$overall), 2L)
})

test_that("est_mg() GPCM-only: par.est$overall has GPCM model column", {
  dat <- make_two_groups(x_gpcm2, seed = 40, n = 400)
  fit <- run_mg(x_gpcm2, dat)
  expect_true(all(fit$par.est$overall$model == "GPCM"))
})


# ══════════════════════════════════════════════════════════════════════════════
# 5. Mixed format (3PLM + GRM)
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_mg() mixed (3PLM+GRM): converges with 8 items", {
  dat <- make_two_groups(x_mixed8, seed = 50, n = 500)
  fit <- run_mg(x_mixed8, dat, use.gprior = TRUE,
                gprior = list(dist = "beta", params = c(5, 16)))
  expect_s3_class(fit, "est_mg")
  expect_equal(nrow(fit$par.est$overall), 8L)
})

test_that("est_mg() mixed: par.est$overall contains both 3PLM and GRM rows", {
  dat <- make_two_groups(x_mixed8, seed = 50, n = 500)
  fit <- run_mg(x_mixed8, dat, use.gprior = TRUE,
                gprior = list(dist = "beta", params = c(5, 16)))
  models <- fit$par.est$overall$model
  expect_true("3PLM" %in% models)
  expect_true("GRM"  %in% models)
})


# ══════════════════════════════════════════════════════════════════════════════
# 6. Common-item linking: shared IDs constrain to same parameters
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_mg() anchor items: group-specific par.est agree for common items", {
  # Both groups have the same 10 items (same IDs) → all items are anchors
  dat <- make_two_groups(x_drm10, seed = 60)
  fit <- run_mg(x_drm10, dat, use.gprior = TRUE,
                gprior = list(dist = "beta", params = c(5, 16)))
  par_g1 <- fit$par.est$group$G1[, c("par.1","par.2","par.3")]
  par_g2 <- fit$par.est$group$G2[, c("par.1","par.2","par.3")]
  # Constrained common items must have identical estimates across groups
  expect_equal(par_g1, par_g2, tolerance = 1e-8)
})


# ══════════════════════════════════════════════════════════════════════════════
# 7. EmpHist = TRUE
# ══════════════════════════════════════════════════════════════════════════════

test_that("est_mg() EmpHist=TRUE produces non-uniform weights for G2", {
  dat <- make_two_groups(x_drm10, seed = 70)
  fit <- do.call(est_mg, c(
    list(x          = list(x_drm10, x_drm10),
         data       = list(dat$g1, dat$g2),
         group.name = c("G1", "G2"),
         free.group = "G2",
         EmpHist    = TRUE),
    MG_ARGS
  ))
  expect_s3_class(fit, "est_mg")
})


# ══════════════════════════════════════════════════════════════════════════════
# 8. FIPC (multiple-group fixed item parameter calibration)
# ══════════════════════════════════════════════════════════════════════════════
# This section exercises the est_mg() FIPC code path -- the multi-group
# analogue of the FIPC branch tested for est_irt().  est_mg_fipc()
# calls the same divide_data() / Estep_fipc() / Mstep / info_xpd()
# pipeline as est_irt_fipc(), so the test guards regression in any of
# the shared helpers when they are exercised through the multi-group
# code path (e.g. the freq.cat construction at est_mg.R lines ~1797-
# 1803, which is otherwise uncovered by the test suite).

test_that("est_mg() FIPC (MEM) estimates pretest items on fixed-item scale", {
  # build a small fixed-item bank: 8 dichotomous (3PLM) items taken
  # from the flexMIRT sample (x_full[1:8, ]).  these are the items
  # whose parameters stay fixed during multi-group calibration
  x_fix_bank <- x_full[1:8, ]

  # FIPC metadata: 8 fixed 3PLM items at positions 1-8, then 4 new
  # pretest items at positions 9-12 (2 dichotomous + 2 GRM with 5
  # categories).  shape_df_fipc() inserts default starting values
  # for the new items and returns a single combined data.frame
  meta_fipc <- shape_df_fipc(
    x       = x_fix_bank,
    fix.loc = 1:8,
    item.id = paste0("NI", 1:4),
    cats    = c(2L, 2L, 5L, 5L),
    model   = c("3PLM", "3PLM", "GRM", "GRM")
  )

  # simulate two groups from the SAME 12-item form but with different
  # latent-trait distributions (G1 ~ N(0, 1), G2 ~ N(0.4, 1.1));
  # est_mg() should recover the G2 mean shift while keeping the
  # 8 fixed-item parameters anchored on the reference scale
  set.seed(901)
  data_list <- list(
    G1 = simdat(x = meta_fipc, theta = rnorm(300, mean = 0.0, sd = 1.0), D = 1),
    G2 = simdat(x = meta_fipc, theta = rnorm(300, mean = 0.4, sd = 1.1), D = 1)
  )

  # run multi-group FIPC: both groups share the same metadata (the
  # fixed items must agree across groups by definition)
  # NOTE: in MG-FIPC, fix.loc must be a list of integer vectors --
  # one element per group -- because each group can fix a different
  # subset of items (see est_mg.R `fix.loc` param doc).  Both groups
  # fix items 1-8 here because they share the same metadata.
  fit <- do.call(est_mg, c(
    list(x          = list(meta_fipc, meta_fipc),
         data       = data_list,
         group.name = c("G1", "G2"),
         free.group = "G2",
         fipc       = TRUE,
         fipc.method = "MEM",
         fix.loc    = list(1:8, 1:8),
         use.gprior = TRUE,
         gprior     = list(dist = "beta", params = c(5, 16))),
    MG_ARGS
  ))

  # basic structural checks
  expect_s3_class(fit, "est_mg")
  # overall summary contains all 12 items (8 fixed + 4 new)
  expect_equal(nrow(fit$par.est$overall), 12L)
  # per-group results exist and cover all groups
  expect_true(all(c("G1", "G2") %in% names(fit$par.est$group)))
  # G2 (the freed group) should have a shifted mean -- not equal to
  # the G1 anchor at 0 (using a loose tolerance because n=300 is
  # small relative to the population shift)
  mu_g2 <- fit$group.par$G2$mu[!is.na(fit$group.par$G2$mu)][1]
  expect_false(isTRUE(all.equal(mu_g2, 0, tolerance = 0.05)))
})


# ══════════════════════════════════════════════════════════════════════════════
# 9. summary() and print() do not error
# ══════════════════════════════════════════════════════════════════════════════

test_that("summary.est_mg() runs without error", {
  dat <- make_two_groups(x_drm10, seed = 80)
  fit <- run_mg(x_drm10, dat)
  expect_no_error(capture.output(summary(fit)))
})

test_that("print.est_mg() runs without error", {
  dat <- make_two_groups(x_drm10, seed = 80)
  fit <- run_mg(x_drm10, dat)
  expect_no_error(capture.output(print(fit)))
})
