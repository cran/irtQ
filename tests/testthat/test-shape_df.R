# shape_df() builds the item metadata, including the default parameter values

test_that("shape_df() with default.par repeats a single cats value for every item", {
  models <- c("1PLM", "2PLM", "2PLM", "3PLM", "3PLM")
  meta1 <- shape_df(cats = 2, model = models, default.par = TRUE)
  # the same input written out in full
  meta2 <- shape_df(cats = rep(2, 5), model = models, default.par = TRUE)
  expect_identical(meta1, meta2)
  # item ids follow the number of items
  expect_identical(meta1$id, paste0("V", 1:5))
  # guessing parameters are NA for 1PLM and 2PLM items and 0.2 for 3PLM items
  expect_identical(meta1$par.3, c(NA, NA, NA, 0.2, 0.2))
})

test_that("shape_df() with default.par repeats a single cats value for polytomous items", {
  meta1 <- shape_df(cats = 3, model = c("GRM", "GPCM"), default.par = TRUE)
  meta2 <- shape_df(cats = c(3, 3), model = c("GRM", "GPCM"), default.par = TRUE)
  expect_identical(meta1, meta2)
  expect_identical(meta1$id, c("V1", "V2"))
  expect_identical(meta1$cats, c(3, 3))
})

test_that("shape_df() with default.par takes the number of items from item.id", {
  meta1 <- shape_df(cats = 2, model = "3PLM", item.id = paste0("it", 1:3), default.par = TRUE)
  expect_identical(meta1$id, paste0("it", 1:3))
  expect_identical(meta1$par.3, rep(0.2, 3))
})

test_that("shape_df() with default.par repeats a single model for every item", {
  meta1 <- shape_df(cats = rep(2, 4), model = "3PLM", default.par = TRUE)
  meta2 <- shape_df(cats = rep(2, 4), model = rep("3PLM", 4), default.par = TRUE)
  expect_identical(meta1, meta2)
  expect_identical(meta1$par.3, rep(0.2, 4))
})

test_that("startval_df() repeats a single cats value like a full vector", {
  models <- c("1PLM", "3PLM", "DRM")
  expect_identical(
    irtQ:::startval_df(cats = 2, model = models),
    irtQ:::startval_df(cats = rep(2, 3), model = models)
  )
})
