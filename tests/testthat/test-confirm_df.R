# confirm_df() is not exported; access via :::
confirm_df <- irtQ:::confirm_df

# ── helpers ───────────────────────────────────────────────────────────────────

make_drm_df <- function(model = "3PLM") {
  data.frame(
    id    = c("I1", "I2"),
    cats  = c(2L, 2L),
    model = model,
    par.1 = c(1.0, 1.2),
    par.2 = c(-0.5, 0.5),
    par.3 = c(0.2, 0.15),
    stringsAsFactors = FALSE
  )
}

make_grm_df <- function() {
  data.frame(
    id    = c("I1", "I2"),
    cats  = c(4L, 4L),
    model = "GRM",
    par.1 = c(1.0, 1.2),
    par.2 = c(-1.0, -0.5),
    par.3 = c(0.0, 0.5),
    par.4 = c(1.0, 1.5),
    stringsAsFactors = FALSE
  )
}

# ── valid inputs ──────────────────────────────────────────────────────────────

test_that("confirm_df() accepts a valid 3PLM data frame", {
  x <- make_drm_df("3PLM")
  result <- confirm_df(x)
  expect_s3_class(result, "data.frame")
  expect_equal(colnames(result)[1:3], c("id", "cats", "model"))
})

test_that("confirm_df() accepts a valid 2PLM data frame and sets par.3 = 0", {
  x <- data.frame(
    id = "I1", cats = 2L, model = "2PLM",
    par.1 = 1.0, par.2 = 0.0, par.3 = NA_real_,
    stringsAsFactors = FALSE
  )
  result <- confirm_df(x)
  expect_equal(result$par.3, 0)
})

test_that("confirm_df() accepts a valid 1PLM data frame and sets par.3 = 0", {
  x <- data.frame(
    id = "I1", cats = 2L, model = "1PLM",
    par.1 = 1.0, par.2 = 0.5, par.3 = NA_real_,
    stringsAsFactors = FALSE
  )
  result <- confirm_df(x)
  expect_equal(result$par.3, 0)
})

test_that("confirm_df() accepts a valid GRM data frame", {
  x <- make_grm_df()
  result <- confirm_df(x)
  expect_equal(result$model, c("GRM", "GRM"))
})

test_that("confirm_df() accepts a valid GPCM data frame", {
  x <- data.frame(
    id = "I1", cats = 4L, model = "GPCM",
    par.1 = 1.2, par.2 = -0.5, par.3 = 0.0, par.4 = 0.8,
    stringsAsFactors = FALSE
  )
  expect_no_error(confirm_df(x))
})

test_that("confirm_df() converts DRM to 3PLM with a warning", {
  x <- make_drm_df("DRM")
  expect_warning(result <- confirm_df(x), "DRM")
  expect_true(all(result$model == "3PLM"))
})

test_that("confirm_df() normalises model names to upper case", {
  x <- make_drm_df("3plm")
  result <- confirm_df(x)
  expect_equal(result$model, c("3PLM", "3PLM"))
})

test_that("confirm_df() renames columns to id/cats/model/par.1/par.2/...", {
  x <- make_drm_df("2PLM")
  x$par.3 <- NA_real_
  colnames(x) <- c("item", "ncat", "irt", "a", "b", "g")
  result <- confirm_df(x)
  expect_equal(colnames(result)[1:3], c("id", "cats", "model"))
  expect_true("par.1" %in% colnames(result))
})

test_that("confirm_df() g2na=TRUE sets par.3 to NA for 1PLM/2PLM items", {
  x <- make_drm_df("2PLM")
  x$par.3 <- NA_real_
  result <- confirm_df(x, g2na = TRUE)
  expect_true(all(is.na(result$par.3)))
})

test_that("confirm_df() drops all-NA trailing parameter columns", {
  x <- make_drm_df("3PLM")
  x$par.4 <- NA_real_
  x$par.5 <- NA_real_
  result <- confirm_df(x)
  expect_false("par.5" %in% colnames(result))
})

test_that("confirm_df() accepts a factor model column and converts it", {
  x <- make_drm_df("3PLM")
  x$model <- factor(x$model)
  result <- confirm_df(x)
  expect_true(is.character(result$model))
})

# ── invalid inputs ────────────────────────────────────────────────────────────

test_that("confirm_df() errors on unknown model name", {
  x <- make_drm_df("4PLM")
  expect_error(confirm_df(x), "mis-specified")
})

test_that("confirm_df() errors when cats < 1 (cats = 0)", {
  x <- make_drm_df("3PLM")
  x$cats[1] <- 0L
  expect_error(confirm_df(x), "score category")
})

test_that("confirm_df() errors when par.3 column missing and model is not all-2PLM", {
  x <- data.frame(
    id = "I1", cats = 2L, model = "3PLM",
    par.1 = 1.0, par.2 = 0.0,
    stringsAsFactors = FALSE
  )
  expect_error(confirm_df(x), "par.3")
})
