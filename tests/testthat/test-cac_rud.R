# cac_rud() computes classification accuracy and consistency from ability estimates and standard errors

set.seed(11)
theta_cac <- rnorm(50)
se_cac <- rep(0.3, 50)

test_that("cac_rud() labels the total row of the marginal table", {
  res <- cac_rud(cutscore = c(-0.5, 0.8), theta = theta_cac, se = se_cac)
  expect_equal(nrow(res$marginal), 4L)
  expect_identical(as.character(res$marginal$level), c("1", "2", "3", "marginal"))
  expect_false(anyNA(res$marginal$level))
})

test_that("cac_rud() marginal accuracy and consistency are the weighted sums of the conditional values", {
  res <- cac_rud(cutscore = c(-0.5, 0.8), theta = theta_cac, se = se_cac)
  total <- res$marginal[res$marginal$level == "marginal", ]
  expect_equal(total$accuracy, mean(res$conditional$accuracy))
  expect_equal(total$consistency, mean(res$conditional$consistency))
})

test_that("cac_rud() handles a performance level that no examinee reaches", {
  res <- cac_rud(cutscore = c(-0.5, 10), theta = theta_cac, se = se_cac)
  # the empty level has zero accuracy and consistency
  expect_equal(nrow(res$marginal), 4L)
  expect_equal(res$marginal$accuracy[3], 0)
  expect_equal(res$marginal$consistency[3], 0)
  # the confusion matrix keeps one row and one column per level
  expect_equal(dim(res$confusion), c(3L, 3L))
  expect_equal(unname(res$confusion[3, ]), c(0, 0, 0))
  # the marginal values still add up across the occupied levels
  expect_equal(res$marginal$accuracy[4], sum(res$marginal$accuracy[1:3]))
})

test_that("cac_rud() handles an empty level when quadrature weights are supplied", {
  wts <- gen.weight(n = 41, dist = "norm", mu = 0, sigma = 1)
  res <- cac_rud(cutscore = c(-0.5, 10), weights = wts, se = rep(0.3, 41))
  expect_equal(dim(res$confusion), c(3L, 3L))
  expect_equal(nrow(res$marginal), 4L)
  expect_lt(res$marginal$accuracy[3], 1e-10)
})
