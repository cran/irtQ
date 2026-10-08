# fixtures: the simMST 1-3-3 panel (dichotomous 3PLM items on the D = 1.702 scale)

x_mst     <- simMST$item_bank
mod_mst   <- simMST$module
map_mst   <- simMST$route_map
cut_mst   <- simMST$cut_score

# responses of 30 examinees to all 56 items of the bank
set.seed(2027)
theta_mst <- rnorm(30)
resp_mst  <- simdat(x = x_mst, theta = theta_mst, D = 1.702)

# fixtures: a small mixed-format 1-2-2 panel (three 3PLM items and one GRM item
# per module), so that every module has the same maximum sum score of 5
set.seed(2028)
n_mod_mix <- 5L
x_mix <- shape_df(
  par.drm = list(a = runif(3 * n_mod_mix, 0.8, 1.6),
                 b = rnorm(3 * n_mod_mix, 0, 0.8),
                 g = rep(0.15, 3 * n_mod_mix)),
  par.prm = list(a = runif(n_mod_mix, 0.8, 1.4),
                 d = replicate(n_mod_mix, sort(rnorm(2, 0, 0.7)), simplify = FALSE)),
  cats    = c(rep(2L, 3 * n_mod_mix), rep(3L, n_mod_mix)),
  model   = c(rep("3PLM", 3 * n_mod_mix), rep("GRM", n_mod_mix))
)
# items 1-3 of module 1, 4-6 of module 2, ...; the GRM items are 16-20
mod_mix <- matrix(0L, nrow = nrow(x_mix), ncol = n_mod_mix)
for (m in seq_len(n_mod_mix)) {
  mod_mix[c(((m - 1L) * 3L + 1L):(m * 3L), 15L + m), m] <- 1L
}
map_mix <- matrix(0L, n_mod_mix, n_mod_mix)
map_mix[1, 2:3] <- 1L
map_mix[2:3, 4:5] <- 1L
set.seed(2029)
theta_mix <- rnorm(30)
resp_mix  <- simdat(x = x_mix, theta = theta_mix, D = 1)

# helper: item rows of the modules administered up to stage s for examinee i
items_upto <- function(path_row, module, s) {
  unlist(lapply(path_row[seq_len(s)], function(m) which(module[, m] == 1)))
}

# helper: max abs difference between the routing estimates of run_mst() and
# est_score() applied to the cumulative responses of stages 1..s
max_diff_cum <- function(fit, x, module, resp, D, method, stages = 1:2) {
  diffs <- vapply(seq_len(nrow(fit$path)), function(i) {
    max(vapply(stages, function(s) {
      it  <- items_upto(fit$path[i, ], module, s)
      ref <- est_score(x = x[it, ], data = resp[i, it, drop = FALSE],
                       D = D, method = method)$est.theta
      abs(ref - fit$theta.route[i, s])
    }, numeric(1L)))
  }, numeric(1L))
  max(diffs)
}


# 1. structure of the result

test_that("run_mst() returns the documented structure for every routing method", {
  for (meth in c("ML", "WL", "MLF", "MAP", "EAP", "EAP.SUM", "INV.TCC")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst[1:10], response = resp_mst[1:10, ], D = 1.702,
      route_method = "bmat",
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_s3_class(fit, "run_mst")
    expect_equal(dim(fit$theta.route), c(10L, 3L))
    expect_equal(dim(fit$path), c(10L, 3L))
    # routing estimates of stages 1 and 2 are finite
    expect_true(all(is.finite(fit$theta.route[, 1:2])))
    # the last column of theta.route is the final estimate
    expect_equal(unname(fit$theta.route[, 3]), fit$est.theta)
  }
})


# 2. routing estimates use the responses to all modules administered so far

test_that("ML and EAP routing estimates equal est_score() on the cumulative responses", {
  for (meth in c("ML", "EAP")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst, response = resp_mst, D = 1.702,
      route_method = "bmat",
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_lt(max_diff_cum(fit, x_mst, mod_mst, resp_mst, 1.702, meth), 1e-6)
  }
})

test_that("WL, MAP, and MLF routing estimates equal est_score() on the cumulative responses", {
  for (meth in c("WL", "MAP", "MLF")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst, response = resp_mst, D = 1.702,
      route_method = NULL, cut_score = cut_mst,
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_lt(max_diff_cum(fit, x_mst, mod_mst, resp_mst, 1.702, meth), 1e-6)
  }
})

test_that("the stage 2 routing estimate differs from the estimate based on the stage 2 module alone", {
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_mst, D = 1.702,
    route_method = "bmat",
    route_score = list(method = "ML"),
    final_score = list(method = "ML"), verbose = FALSE
  )
  # estimate from the responses to the stage 2 module only
  stage_only <- vapply(seq_len(nrow(fit$path)), function(i) {
    it <- which(mod_mst[, fit$path[i, 2]] == 1)
    est_score(x = x_mst[it, ], data = resp_mst[i, it, drop = FALSE],
              D = 1.702, method = "ML")$est.theta
  }, numeric(1L))
  expect_gt(max(abs(stage_only - fit$theta.route[, 2])), 1e-3)
})

test_that("routing estimates are cumulative with mixed-format modules", {
  for (meth in c("ML", "EAP")) {
    fit <- run_mst(
      x = x_mix, route_map = map_mix, module = mod_mix,
      theta = theta_mix, response = resp_mix, D = 1,
      route_method = "bmat",
      route_score = list(method = meth),
      final_score = list(method = "ML"), verbose = FALSE
    )
    expect_lt(max_diff_cum(fit, x_mix, mod_mix, resp_mix, 1, meth, stages = 1:2), 1e-6)
  }
})

test_that("routing estimates use only the observed responses when some are missing", {
  resp_na <- resp_mst
  # set about 10 percent of the responses to missing at fixed positions
  set.seed(2030)
  na_pos <- sample(length(resp_na), size = round(0.10 * length(resp_na)))
  resp_na[na_pos] <- NA
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_na, D = 1.702,
    route_method = "bmat",
    route_score = list(method = "ML"),
    final_score = list(method = "ML"), verbose = FALSE
  )
  # est_score() drops the missing responses and scores the remaining items of
  # the administered modules
  expect_lt(max_diff_cum(fit, x_mst, mod_mst, resp_na, 1.702, "ML"), 1e-6)
})

test_that("final estimates use the parameters of the observed items when some responses are missing", {
  resp_na <- resp_mst
  # set about 10 percent of the responses to missing at fixed positions
  set.seed(2033)
  resp_na[sample(length(resp_na), size = round(0.10 * length(resp_na)))] <- NA
  for (meth in c("ML", "EAP", "MLF")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst, response = resp_na, D = 1.702,
      route_method = "bmat",
      route_score = list(method = "EAP"),
      final_score = list(method = meth), verbose = FALSE
    )
    # est_score() on the observed items of the administered modules
    ref <- vapply(seq_len(nrow(fit$path)), function(i) {
      it <- items_upto(fit$path[i, ], mod_mst, 3L)
      est_score(x = x_mst[it, ], data = resp_na[i, it, drop = FALSE],
                D = 1.702, method = meth)$est.theta
    }, numeric(1L))
    expect_equal(fit$est.theta, ref, tolerance = 1e-6)
  }
})

test_that("final sum-score estimates count missing responses as zero, as est_score() does", {
  resp_na <- resp_mst
  set.seed(2034)
  resp_na[sample(length(resp_na), size = round(0.10 * length(resp_na)))] <- NA
  for (meth in c("EAP.SUM", "INV.TCC")) {
    fit <- run_mst(
      x = x_mst, route_map = map_mst, module = mod_mst,
      theta = theta_mst, response = resp_na, D = 1.702,
      route_method = "bmat",
      route_score = list(method = "EAP"),
      final_score = list(method = meth), verbose = FALSE
    )
    # est_score() replaces missing responses with zeros and warns
    ref <- suppressWarnings(vapply(seq_len(nrow(fit$path)), function(i) {
      it <- items_upto(fit$path[i, ], mod_mst, 3L)
      est_score(x = x_mst[it, ], data = resp_na[i, it, drop = FALSE],
                D = 1.702, method = meth)$est.par$est.theta
    }, numeric(1L)))
    expect_equal(fit$est.theta, ref, tolerance = 1e-6)
  }
})

test_that("INV.TCC routing estimates equal the inverse TCC lookup of the cumulative sum score", {
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_mst, D = 1.702,
    route_method = NULL, cut_score = cut_mst,
    route_score = list(method = "INV.TCC"),
    final_score = list(method = "INV.TCC"), verbose = FALSE
  )
  for (i in seq_len(nrow(fit$path))) {
    for (s in 1:2) {
      it  <- items_upto(fit$path[i, ], mod_mst, s)
      tbl <- est_score(x = x_mst[it, ], data = resp_mst[i, it, drop = FALSE],
                       D = 1.702, method = "INV.TCC")$score.table
      ref <- tbl$est.theta[tbl$sum.score == sum(resp_mst[i, it])]
      expect_equal(unname(fit$theta.route[i, s]), ref, tolerance = 1e-8)
    }
  }
})

test_that("cut-score routing assigns the module from the cumulative estimate", {
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_mst, response = resp_mst, D = 1.702,
    route_method = NULL, cut_score = cut_mst,
    route_score = list(method = "ML"),
    final_score = list(method = "ML"), verbose = FALSE
  )
  # stage 2 module follows the stage 1 estimate
  rank2 <- findInterval(fit$theta.route[, 1], cut_mst[[1]]) + 1L
  expect_equal(unname(fit$path[, 2]), c(2L, 3L, 4L)[rank2])

  # stage 3 module follows the estimate that pools the stage 1 and stage 2
  # responses; modules 2 and 4 reach only two of the three stage 3 modules, so
  # only the cut score that separates those two modules applies
  cut3 <- cut_mst[[2]]
  est2 <- unname(fit$theta.route[, 2])
  expected3 <- vapply(seq_len(nrow(fit$path)), function(i) {
    switch(as.character(fit$path[i, 2]),
      "2" = if (est2[i] <= cut3[1]) 5L else 6L,
      "3" = if (est2[i] <= cut3[1]) 5L else if (est2[i] <= cut3[2]) 6L else 7L,
      "4" = if (est2[i] <= cut3[2]) 6L else 7L)
  }, integer(1L))
  expect_equal(unname(fit$path[, 3]), expected3)
})

test_that("an examinee in module 4 with a middle estimate is routed to module 6", {
  # inverse TCC routing keeps the simulation fast for many examinees
  set.seed(2032)
  fit <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = rnorm(6000), D = 1.702,
    route_method = NULL, cut_score = cut_mst,
    route_score = list(method = "INV.TCC"),
    final_score = list(method = "INV.TCC"), verbose = FALSE
  )
  cut3    <- cut_mst[[2]]
  est2    <- fit$theta.route[, 2]
  in_mod4 <- fit$path[, 2] == 4L
  middle  <- in_mod4 & est2 > cut3[1] & est2 <= cut3[2]
  # the simulation has examinees in module 4 whose estimate lies between the cut scores
  expect_gt(sum(middle), 0L)
  # they go to module 6, the easier of the two modules reachable from module 4
  expect_true(all(fit$path[middle, 3] == 6L))
  # module 4 examinees above the second cut score go to module 7
  expect_true(all(fit$path[in_mod4 & est2 > cut3[2], 3] == 7L))
})


# 3. agreement with the recursion-based evaluation

test_that("run_mst() with cut scores and inverse TCC scoring approaches reval_mst()", {
  skip_on_cran()

  # simulation with 1000 examinees at each of five ability levels
  grid <- seq(-2, 2, 1)
  n_rep <- 1000L
  set.seed(2031)
  theta_rep <- rep(grid, each = n_rep)
  mc <- run_mst(
    x = x_mst, route_map = map_mst, module = mod_mst,
    theta = theta_rep, D = 1.702,
    route_method = NULL, cut_score = cut_mst,
    route_score = list(method = "INV.TCC", range.tcc = c(-7, 7)),
    final_score = list(method = "INV.TCC", range.tcc = c(-7, 7)),
    verbose = FALSE
  )
  rv <- reval_mst(
    x = x_mst, D = 1.702, route_map = map_mst, module = mod_mst,
    cut_score = cut_mst, theta = grid, range.tcc = c(-7, 7)
  )$eval.tb

  grp       <- factor(theta_rep)
  bias_mc   <- as.numeric(tapply(mc$est.theta - theta_rep, grp, mean))
  csem_mc   <- as.numeric(tapply(mc$est.theta, grp, stats::sd))

  # bias differs from the analytical value within four Monte Carlo standard errors
  expect_true(all(abs(bias_mc - rv$bias) < 4 * rv$csem / sqrt(n_rep)))
  # CSEM differs from the analytical value by less than 15 percent
  expect_true(all(abs(csem_mc - rv$csem) / rv$csem < 0.15))
})


# 4. input validation

test_that("run_mst() stops on invalid input", {
  expect_error(
    run_mst(x = x_mst, route_map = map_mst, module = mod_mst, verbose = FALSE),
    "At least one of"
  )
  expect_error(
    run_mst(x = x_mst, route_map = map_mst, module = mod_mst,
            theta = theta_mst, route_method = NULL, verbose = FALSE),
    "cut_score"
  )
  expect_error(
    run_mst(x = x_mst, route_map = map_mst, module = mod_mst,
            theta = theta_mst, route_score = list(method = "XYZ"),
            verbose = FALSE),
    "route_score"
  )
})


# 5. routing at a cut score

test_that("give_path() sends a score equal to a cut score to the higher module", {
  # scores at, below, and above the cut scores, including infinite scores
  out <- irtQ:::give_path(score = c(-0.5, 0, 0.5, -Inf, Inf),
                          cut_sc = c(-0.5, 0.5))
  expect_equal(out$path, c(2, 2, 3, 1, 3))

  # a single cut score: the cut score itself goes to the higher module
  out1 <- irtQ:::give_path(score = c(-1, 0, 1), cut_sc = 0)
  expect_equal(out1$path, c(1, 2, 2))

  # no cut score: every score, including an infinite one, goes to the only module
  out0 <- irtQ:::give_path(score = c(-Inf, 0, Inf), cut_sc = numeric(0))
  expect_equal(out0$path, c(1, 1, 1))

  # a missing score stays missing
  expect_true(is.na(irtQ:::give_path(score = NA_real_, cut_sc = 0)$path))
})
