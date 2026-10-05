# shape_df_fipc() builds the metadata of fixed and new items for FIPC

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_fix <- bring.flexmirt(file = prm_file, "par")$Group1$full_df[1:10, ]

test_that("shape_df_fipc() repeats a single cats and model for every new item", {
  new_ids <- paste0("N", 1:4)
  # single values
  meta1 <- shape_df_fipc(x = x_fix, fix.loc = 1:10, item.id = new_ids, cats = 2, model = "3PLM")
  # the same input written out in full
  meta2 <- shape_df_fipc(
    x = x_fix, fix.loc = 1:10, item.id = new_ids,
    cats = rep(2, 4), model = rep("3PLM", 4)
  )
  expect_identical(meta1, meta2)
  expect_equal(nrow(meta1), 14L)
  expect_identical(as.character(meta1$id[11:14]), new_ids)
  expect_true(all(meta1$model[11:14] == "3PLM"))
})

test_that("shape_df_fipc() stops when cats or model do not match the number of new items", {
  new_ids <- paste0("N", 1:3)
  expect_error(
    shape_df_fipc(x = x_fix, fix.loc = 1:10, item.id = new_ids, cats = c(2, 2), model = "3PLM"),
    "must be 1 or equal to the number of new items"
  )
  expect_error(
    shape_df_fipc(x = x_fix, fix.loc = 1:10, item.id = new_ids, cats = 2, model = c("3PLM", "2PLM")),
    "must be 1 or equal to the number of new items"
  )
})
