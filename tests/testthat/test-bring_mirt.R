# bring.mirt() needs the suggested mirt package, so skip when it is missing

test_that("bring.mirt() returns the item metadata of a fitted mirt model", {
  skip_if_not_installed("mirt")

  # fit the 2PL model to the LSAT6 data with mirt
  fit <- suppressMessages(
    mirt::mirt(as.data.frame(LSAT6), 1, itemtype = "2PL", verbose = FALSE)
  )
  out <- bring.mirt(fit)
  meta <- out$full_df

  # one row per item, with the item metadata columns
  expect_s3_class(meta, "data.frame")
  expect_equal(nrow(meta), 5L)
  expect_true(all(c("id", "cats", "model", "par.1", "par.2") %in% names(meta)))
  expect_true(all(meta$cats == 2L))
  expect_true(all(meta$model == "DRM"))
})
