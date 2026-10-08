# fixtures: the simMST 1-3-3 panel (dichotomous 3PLM items on the D = 1.702 scale)

x_rv     <- simMST$item_bank
module_rv <- simMST$module
map_rv   <- simMST$route_map
cut_rv   <- simMST$cut_score


test_that("reval_mst() evaluates the default theta grid and returns the documented objects", {
  rv <- reval_mst(x = x_rv, D = 1.702, route_map = map_rv, module = module_rv,
                  cut_score = cut_rv)
  expect_named(rv, c("panel.info", "item.by.mod", "item.by.path", "eq.theta",
                     "cdist.by.mod", "jdist.by.path", "eval.tb"))
  # the default theta grid is seq(-5, 5, 1)
  expect_equal(rv$eval.tb$theta, seq(-5, 5, 1))
  expect_true(all(c("theta", "mu", "sigma2", "bias", "csem") %in% names(rv$eval.tb)))
  expect_true(all(is.finite(rv$eval.tb$csem)))
})

test_that("reval_mst() stops when modules in a stage differ in maximum sum score", {
  # drop one item from module 2 so that stage 2 has modules of 7, 8, and 8 items
  module_bad <- module_rv
  module_bad[which(module_bad[, 2] == 1)[1], 2] <- 0L
  expect_error(
    reval_mst(x = x_rv, D = 1.702, route_map = map_rv, module = module_bad,
              cut_score = cut_rv, theta = seq(-1, 1, 1)),
    "same maximum sum score"
  )
})
