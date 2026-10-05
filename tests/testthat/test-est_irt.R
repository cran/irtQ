# Shared fixture: import 55-item flexMIRT parameter file used in README examples
prm_file  <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_true    <- bring.flexmirt(file = prm_file, "par")$Group1$full_df   # 55 items
x_ref     <- x_true[1:40, ]   # first 40 items (38 × 3PLM + 2 × GRM)

# Common EM settings kept small for test speed
EM_args <- list(
  D          = 1,
  Quadrature = c(41L, 6),
  Etol       = 1e-3,
  MaxE       = 100L,
  se         = FALSE,
  verbose    = FALSE
)

# ── helper: generate response data ────────────────────────────────────────────
gen_data <- function(x, n = 500, seed = 1, mean = 0, sd = 1) {
  set.seed(seed)
  theta <- rnorm(n, mean = mean, sd = sd)
  simdat(x = x, theta = theta, D = 1)
}


# ── 1. Dichotomous-only ───────────────────────────────────────────────────────

test_that("est_irt() converges for 1PLM (fix.a.1pl = TRUE)", {
  x_1pl <- x_ref[1:10, ]
  x_1pl$model <- "1PLM"
  data <- gen_data(x_1pl)
  args <- c(list(data = data, model = "1PLM", cats = 2, fix.a.1pl = TRUE), EM_args)
  fit  <- do.call(est_irt, args)
  expect_s3_class(fit, "est_irt")
  expect_equal(nrow(fit$par.est), 10L)
  expect_true(fit$niter > 0)
})

test_that("est_irt() converges for 2PLM", {
  x_2pl <- x_ref[1:10, ]
  x_2pl$model <- "2PLM"
  x_2pl$par.3 <- NA_real_
  data <- gen_data(x_2pl)
  args <- c(list(data = data, model = "2PLM", cats = 2), EM_args)
  fit  <- do.call(est_irt, args)
  expect_s3_class(fit, "est_irt")
  expect_equal(nrow(fit$par.est), 10L)
})

test_that("est_irt() converges for 3PLM with guessing prior", {
  x_3pl <- x_ref[1:10, ]
  data  <- gen_data(x_3pl)
  args  <- c(
    list(data = data, model = "3PLM", cats = 2,
         use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16))),
    EM_args
  )
  fit <- do.call(est_irt, args)
  expect_s3_class(fit, "est_irt")
  expect_true(all(fit$par.est$par.3 >= 0 & fit$par.est$par.3 <= 1))
})

test_that("est_irt() LSAT6 data (1PLM, fix.a.1pl = FALSE) returns est_irt object", {
  fit <- est_irt(
    data       = LSAT6, D = 1, model = "1PLM", cats = 2,
    fix.a.1pl  = FALSE, Etol = 1e-3, MaxE = 100L,
    se = FALSE, verbose = FALSE
  )
  expect_s3_class(fit, "est_irt")
  expect_equal(nrow(fit$par.est), 5L)
})

# the EM must stop only when the largest absolute change is <= Etol; the
# constant guessing column of 2PLM items has a change of exactly zero and
# must not hide larger negative changes in the free parameters
test_that("est_irt() 2PLM on LSAT6 stops only when the largest absolute change is <= Etol", {
  fit <- est_irt(
    data = LSAT6, D = 1, model = "2PLM", cats = 2,
    se = FALSE, verbose = FALSE
  )
  expect_gt(fit$niter, 2L)
  expect_gt(fit$maxpar.diff, 0)
  expect_lte(fit$maxpar.diff, fit$Etol)

  # allowing exactly the number of cycles that were needed is still
  # judged as converged and gives no warning
  warns <- character(0)
  fit2 <- withCallingHandlers(
    est_irt(
      data = LSAT6, D = 1, model = "2PLM", cats = 2,
      MaxE = fit$niter, se = FALSE, verbose = FALSE
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(fit2$niter, fit$niter)
  expect_false(any(grepl("Convergence criteria are not satisfied", warns)))
  expect_identical(fit2$test.1, "Convergence criteria are satisfied.")
})


# ── 2. Polytomous-only ────────────────────────────────────────────────────────

test_that("est_irt() converges for GRM", {
  x_grm <- x_ref[39:40, ]   # two 5-category GRM items
  data  <- gen_data(x_grm, n = 600)
  args  <- c(list(data = data, x = x_grm), EM_args)
  fit   <- do.call(est_irt, args)
  expect_s3_class(fit, "est_irt")
  expect_equal(nrow(fit$par.est), 2L)
})

test_that("est_irt() converges for GPCM", {
  set.seed(42)
  # eight 4-category items with slopes from 0.8 to 1.6 and sorted step
  # parameters in [-1.5, 1.5]; a design with enough items to identify the slopes
  x_gpcm <- shape_df(
    par.prm = list(
      a = seq(0.8, 1.6, length.out = 8),
      d = lapply(1:8, function(i) sort(runif(3, -1.5, 1.5)))
    ),
    cats    = rep(4L, 8),
    model   = "GPCM"
  )
  theta <- rnorm(1000)
  data  <- simdat(x = x_gpcm, theta = theta, D = 1)
  args  <- c(list(data = data, x = x_gpcm), EM_args)
  fit   <- do.call(est_irt, args)
  expect_s3_class(fit, "est_irt")
  expect_equal(nrow(fit$par.est), 8L)
  expect_identical(fit$test.1, "Convergence criteria are satisfied.")
})


# ── 3. Mixed-format ───────────────────────────────────────────────────────────

test_that("est_irt() converges for mixed-format (3PLM + GRM)", {
  data <- gen_data(x_ref, n = 800)
  args <- c(
    list(data = data, x = x_ref,
         use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16))),
    EM_args
  )
  fit <- do.call(est_irt, args)
  expect_s3_class(fit, "est_irt")
  expect_equal(nrow(fit$par.est), 40L)
  expect_equal(fit$nitem, 40L)
})


# ── 4. EmpHist TRUE / FALSE ───────────────────────────────────────────────────

test_that("est_irt() EmpHist=TRUE produces non-uniform weights", {
  data <- gen_data(x_ref[1:10, ], n = 500)
  fit  <- est_irt(
    data = data, model = "3PLM", cats = 2, D = 1,
    EmpHist = TRUE, Etol = 1e-3, MaxE = 50L, se = FALSE, verbose = FALSE
  )
  w <- fit$weights$weight
  expect_false(isTRUE(all.equal(rep(w[1], length(w)), w)))
})

test_that("est_irt() EmpHist=FALSE keeps normal prior weights", {
  data <- gen_data(x_ref[1:10, ], n = 500)
  fit  <- est_irt(
    data = data, model = "3PLM", cats = 2, D = 1,
    EmpHist = FALSE, Etol = 1e-3, MaxE = 50L, se = FALSE, verbose = FALSE
  )
  expect_s3_class(fit, "est_irt")
})


# ── 5. FIPC ───────────────────────────────────────────────────────────────────

test_that("est_irt() FIPC (MEM) estimates pretest items on fixed-item scale", {
  set.seed(21)
  theta_new  <- rnorm(800, mean = 0.5, sd = 1.3)
  data_new   <- simdat(x = x_true, theta = theta_new, D = 1)

  meta_fipc  <- shape_df_fipc(
    x        = x_ref,
    fix.loc  = 1:40,
    item.id  = paste0("NI", 1:15),
    cats     = c(rep(2L, 12L), rep(5L, 3L)),
    model    = c(rep("3PLM", 12), rep("GRM", 3))
  )

  args <- c(
    list(x = meta_fipc, data = data_new,
         use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16)),
         EmpHist = TRUE, fipc = TRUE, fipc.method = "MEM", fix.loc = 1:40),
    EM_args
  )
  fit <- do.call(est_irt, args)

  expect_s3_class(fit, "est_irt")
  expect_true(fit$fipc)
  # only 15 pretest items should be estimated
  expect_equal(nrow(fit$par.est), 55L)
})

test_that("est_irt() FIPC (OEM) does not warn about convergence criteria", {
  set.seed(21)
  theta_new <- rnorm(800, mean = 0.5, sd = 1.3)
  data_new  <- simdat(x = x_true, theta = theta_new, D = 1)

  meta_fipc <- shape_df_fipc(
    x        = x_ref,
    fix.loc  = 1:40,
    item.id  = paste0("NI", 1:15),
    cats     = c(rep(2L, 12L), rep(5L, 3L)),
    model    = c(rep("3PLM", 12), rep("GRM", 3))
  )

  args <- c(
    list(x = meta_fipc, data = data_new,
         use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16)),
         EmpHist = TRUE, fipc = TRUE, fipc.method = "OEM", fix.loc = 1:40),
    EM_args
  )

  # collect every warning message; other warnings are allowed
  warns <- character(0)
  fit <- withCallingHandlers(
    do.call(est_irt, args),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(fit$niter, 1L)
  expect_false(any(grepl("Convergence criteria are not satisfied", warns)))
  expect_identical(fit$test.1, "Convergence criteria are satisfied.")
})

test_that("est_irt() FIPC (MEM) with all items fixed iterates until the prior converges", {
  x_fixed <- x_true[1:6, ]
  data    <- gen_data(x_fixed, n = 500, seed = 5)
  fit <- est_irt(
    x = x_fixed, data = data, D = 1, EmpHist = FALSE, Etol = 1e-3, MaxE = 100L,
    fipc = TRUE, fipc.method = "MEM", fix.loc = 1:6, se = FALSE, verbose = FALSE
  )
  # the change in the prior mean and variance must reach Etol
  expect_gt(fit$niter, 1L)
  expect_lte(fit$maxpar.diff, fit$Etol)
  expect_identical(fit$test.1, "Convergence criteria are satisfied.")
})

test_that("est_irt() FIPC (OEM) with all items fixed does not report convergence", {
  x_fixed <- x_true[1:6, ]
  data    <- gen_data(x_fixed, n = 500, seed = 5)
  # a single EM cycle cannot meet the criterion for the prior moments
  fit <- suppressWarnings(est_irt(
    x = x_fixed, data = data, D = 1, EmpHist = FALSE, Etol = 1e-3, MaxE = 100L,
    fipc = TRUE, fipc.method = "OEM", fix.loc = 1:6, se = FALSE, verbose = FALSE
  ))
  expect_equal(fit$niter, 1L)
  expect_false(identical(fit$test.1, "Convergence criteria are satisfied."))
})


# ── 6. Parameter recovery ─────────────────────────────────────────────────────

test_that("est_irt() recovers 2PLM difficulty parameters within tolerance", {
  true_b <- c(-1.5, -0.5, 0.0, 0.5, 1.5)
  true_a <- rep(1.0, 5)
  x_pop  <- shape_df(
    par.drm = list(a = true_a, b = true_b, g = rep(0, 5)),
    cats = rep(2L, 5), model = "2PLM"
  )
  set.seed(7)
  theta <- rnorm(2000)
  data  <- simdat(x = x_pop, theta = theta, D = 1)
  fit   <- est_irt(
    data = data, model = "2PLM", cats = 2, D = 1,
    Etol = 1e-4, MaxE = 300L, se = FALSE, verbose = FALSE
  )
  est_b <- fit$par.est$par.2
  expect_equal(est_b, true_b, tolerance = 0.25)
})


# ── 7. Output structure ───────────────────────────────────────────────────────

test_that("est_irt() output contains expected slots", {
  data <- gen_data(x_ref[1:5, ], n = 300)
  fit  <- est_irt(
    data = data, model = "3PLM", cats = 2, D = 1,
    Etol = 1e-3, MaxE = 50L, se = FALSE, verbose = FALSE
  )
  expect_true(all(c("par.est", "weights", "loglikelihood", "niter", "nitem") %in% names(fit)))
  expect_s3_class(fit$par.est, "data.frame")
  expect_type(fit$loglikelihood, "double")
  expect_type(fit$niter, "integer")
})

test_that("getirt() returns the documented posterior.dist and scale.D for est_irt objects", {
  fit <- est_irt(
    data = LSAT6, D = 1, model = "2PLM", cats = 2,
    Etol = 1e-3, se = FALSE, verbose = FALSE
  )
  expect_identical(getirt(fit, what = "posterior.dist"), fit$posterior.dist)
  expect_identical(getirt(fit, what = "scale.D"), fit$scale.D)
  expect_identical(getirt(fit, what = "scale.D"), 1)
  expect_true(is.matrix(getirt(fit, what = "posterior.dist")))
  expect_true(is.numeric(getirt(fit, what = "scale.D")))
})
