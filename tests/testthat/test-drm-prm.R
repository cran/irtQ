test_that("drm() returns a matrix with correct dimensions", {
  P <- drm(theta = c(-1, 0, 1), a = c(1, 1.5), b = c(0, 0.5), g = c(0, 0.2), D = 1)
  expect_true(is.matrix(P))
  expect_equal(dim(P), c(3L, 2L))
})

test_that("drm() probabilities are bounded between g and 1", {
  g <- c(0.1, 0.2)
  P <- drm(theta = seq(-4, 4, by = 0.5), a = c(1, 1.5), b = c(0, 0.5), g = g, D = 1)
  expect_true(all(P > 0))
  expect_true(all(P < 1))
  # each column's minimum must be >= its guessing parameter
  expect_true(all(P[, 1] >= g[1]))
  expect_true(all(P[, 2] >= g[2]))
})

test_that("drm() 2PL: g=0 gives monotone increasing probabilities", {
  theta <- c(-2, -1, 0, 1, 2)
  P <- drm(theta = theta, a = 1.2, b = 0, D = 1)
  expect_true(all(diff(P[, 1]) > 0))
})

test_that("drm() 1PL: P(theta=b) is 0.5 when g=0", {
  P <- drm(theta = 0.5, a = 1, b = 0.5, g = 0, D = 1)
  expect_equal(as.numeric(P), 0.5, tolerance = 1e-6)
})

test_that("drm() 3PL: P(theta=b) equals (1+g)/2", {
  g_val <- 0.2
  b_val <- 1.0
  P <- drm(theta = b_val, a = 1, b = b_val, g = g_val, D = 1)
  expect_equal(as.numeric(P), (1 + g_val) / 2, tolerance = 1e-6)
})

test_that("drm() with NULL g defaults to 0 (2PL behaviour)", {
  P_null_g <- drm(theta = 0, a = 1, b = 0, D = 1)
  P_zero_g <- drm(theta = 0, a = 1, b = 0, g = 0, D = 1)
  expect_equal(as.numeric(P_null_g), as.numeric(P_zero_g), tolerance = 1e-10)
})

test_that("drm() single theta scalar input works", {
  P <- drm(theta = 0, a = 1, b = 0, g = 0.2, D = 1)
  expect_equal(dim(P), c(1L, 1L))
})

# ── prm() ──────────────────────────────────────────────────────────────────────

test_that("prm() GRM: category probabilities sum to 1", {
  P <- prm(theta = c(-1, 0, 1), a = 1.2, d = c(-1, 0, 1), D = 1, pr.model = "GRM")
  expect_equal(dim(P), c(3L, 4L))
  row_sums <- rowSums(P)
  expect_equal(row_sums, rep(1, 3), tolerance = 1e-8)
})

test_that("prm() GPCM: category probabilities sum to 1", {
  P <- prm(theta = c(-2, 0, 2), a = 1.4, d = c(-0.2, 0, 0.5), D = 1, pr.model = "GPCM")
  expect_equal(dim(P), c(3L, 4L))
  row_sums <- rowSums(P)
  expect_equal(row_sums, rep(1, 3), tolerance = 1e-8)
})

test_that("prm() GRM: all probabilities strictly positive", {
  P <- prm(theta = seq(-3, 3, by = 1), a = 1, d = c(-1.5, 0, 1.5), D = 1, pr.model = "GRM")
  expect_true(all(P > 0))
})

test_that("prm() GPCM: all probabilities strictly positive", {
  P <- prm(theta = seq(-3, 3, by = 1), a = 1, d = c(-0.5, 0.5), D = 1, pr.model = "GPCM")
  expect_true(all(P > 0))
})

test_that("prm() GRM five-category item returns 5 columns", {
  P <- prm(theta = 0, a = 1.2, d = c(-1.5, -0.5, 0.5, 1.5), D = 1, pr.model = "GRM")
  expect_equal(ncol(P), 5L)
})

test_that("prm() GPCM and GRM give different results for same parameters", {
  theta <- c(-1, 0, 1)
  d     <- c(-0.5, 0.5)
  P_grm  <- prm(theta, a = 1, d = d, D = 1, pr.model = "GRM")
  P_gpcm <- prm(theta, a = 1, d = d, D = 1, pr.model = "GPCM")
  expect_false(isTRUE(all.equal(P_grm, P_gpcm)))
})

test_that("drm() known probability: 2PL D=1.702 matches hand calculation", {
  # P = 1 / (1 + exp(-1.702 * 1 * (0 - 0))) = 0.5
  P <- drm(theta = 0, a = 1, b = 0, g = 0, D = 1.702)
  expect_equal(as.numeric(P), 0.5, tolerance = 1e-6)
})
