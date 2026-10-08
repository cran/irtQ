# simdat() simulates item responses from item metadata or from item parameter vectors

test_that("simdat() treats NA guessing parameters as zeros when parameters are given directly", {
  th <- rnorm(200)
  set.seed(1)
  res_na <- simdat(
    theta = th, a.drm = c(1, 1, 1), b.drm = c(0, 0, 0),
    g.drm = c(0.2, NA, NA), cats = c(2, 2, 2), D = 1
  )
  set.seed(1)
  res_zero <- simdat(
    theta = th, a.drm = c(1, 1, 1), b.drm = c(0, 0, 0),
    g.drm = c(0.2, 0, 0), cats = c(2, 2, 2), D = 1
  )
  # no response is missing and NA behaves as a zero guessing parameter
  expect_false(anyNA(res_na))
  expect_identical(res_na, res_zero)
  expect_true(all(res_na %in% c(0, 1)))
})

test_that("simdat() with NA guessing parameters matches the item metadata input", {
  th <- rnorm(200)
  x <- shape_df(
    par.drm = list(a = c(1, 1.2, 0.8), b = c(0, 0.5, -0.5), g = c(0.2, NA, NA)),
    cats = c(2, 2, 2), model = c("3PLM", "2PLM", "1PLM")
  )
  set.seed(5)
  res_meta <- simdat(x = x, theta = th, D = 1)
  set.seed(5)
  res_par <- simdat(
    theta = th, a.drm = c(1, 1.2, 0.8), b.drm = c(0, 0.5, -0.5),
    g.drm = c(0.2, NA, NA), cats = c(2, 2, 2), D = 1
  )
  expect_identical(res_par, res_meta)
})
