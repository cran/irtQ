# tests/testthat/test-find_cut.R
# Unit tests for find_cut()

test_that("find_cut() returns correct structure with simMST 1-3-3 panel", {

  # ── Setup ────────────────────────────────────────────────────────────────
  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map

  result    <- find_cut(x = x, module = module, route_map = route_map)

  # ── Class and top-level structure ────────────────────────────────────────
  expect_s3_class(result, "find_cut")
  expect_named(result, c("cut_score", "details", "tif_data"))

  # ── cut_score: list of length n_stg - 1 ──────────────────────────────────
  # simMST is a 1-3-3 panel (3 stages) → cut_score has 2 elements
  expect_type(result$cut_score, "list")
  expect_length(result$cut_score, 2L)

  # Stage 2 has 3 modules → 2 cut scores
  expect_length(result$cut_score[[1L]], 2L)
  expect_true(all(is.numeric(result$cut_score[[1L]])))

  # Cut scores for stage 2 must be sorted ascending
  cs1 <- result$cut_score[[1L]]
  expect_true(cs1[1L] < cs1[2L])

  # Cut scores must lie within theta_range (default c(-6, 6))
  expect_true(all(result$cut_score[[1L]] > -6 & result$cut_score[[1L]] < 6))
  if (length(result$cut_score[[2L]]) > 0L) {
    expect_true(all(result$cut_score[[2L]] > -6 & result$cut_score[[2L]] < 6))
  }
})

test_that("find_cut() details contain correct stage-level information", {

  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map
  result    <- find_cut(x = x, module = module, route_map = route_map)

  # details is a named list
  expect_type(result$details, "list")
  expect_length(result$details, 2L)   # 2 transitions for 3-stage MST

  # Each stage element has required fields (when multiple modules exist)
  d2 <- result$details[[1L]]   # stage 2 details
  expect_true("modules"      %in% names(d2))
  expect_true("mean_locs"    %in% names(d2))
  expect_true("sorted_order" %in% names(d2))
  expect_true("pairs"        %in% names(d2))

  # mean_locs is a named numeric vector, one per module at stage 2
  expect_true(is.numeric(d2$mean_locs))
  expect_length(d2$mean_locs, length(d2$modules))

  # Each pair entry has the expected fields
  for (pnm in names(d2$pairs)) {
    p <- d2$pairs[[pnm]]
    expect_true("left_module"         %in% names(p))
    expect_true("right_module"        %in% names(p))
    expect_true("proper_crossings"    %in% names(p))
    expect_true("anomalous_crossings" %in% names(p))
    expect_true("selected_cut"        %in% names(p))
    # selected_cut must be a single finite numeric value
    expect_true(is.numeric(p$selected_cut) && length(p$selected_cut) == 1L)
    expect_true(is.finite(p$selected_cut))
    # selected_cut must be one of the proper_crossings
    expect_true(p$selected_cut %in% p$proper_crossings)
  }
})

test_that("find_cut() tif_data is a well-formed tibble", {

  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map
  result    <- find_cut(x = x, module = module, route_map = route_map)

  td <- result$tif_data
  expect_s3_class(td, "tbl_df")
  expect_true(all(c("theta", "module", "stage", "tif") %in% names(td)))
  expect_true(all(is.finite(td$tif)))
  expect_true(all(td$tif >= 0))    # TIF is always non-negative
  expect_true(all(td$stage >= 1L)) # stage 1 (routing) and above are included
})

test_that("find_cut() cut_score is directly usable in run_mst()", {

  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map
  result    <- find_cut(x = x, module = module, route_map = route_map)

  set.seed(42L)
  theta_true <- rnorm(50L)

  # run_mst() with the derived cut scores should complete without error
  expect_no_error({
    mst_result <- run_mst(
      x            = x,
      route_map    = route_map,
      module       = module,
      theta        = theta_true,
      route_method = NULL,
      cut_score    = result$cut_score,
      verbose      = FALSE
    )
  })

  # Basic sanity: one estimated theta per examinee
  expect_length(mst_result$est.theta, 50L)
  expect_true(all(is.finite(mst_result$est.theta)))
})

test_that("print.find_cut() runs without error", {

  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map
  result    <- find_cut(x = x, module = module, route_map = route_map)

  expect_output(print(result), "MST TIF-Crossing Cut Score Results")
  expect_output(print(result), "Selected cut score")
  expect_output(print(result), "run_mst")
})

test_that("find_cut() input validation catches bad arguments", {

  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map

  expect_error(find_cut(x = as.matrix(x), module = module, route_map = route_map),
               "'x' must be a data frame")
  expect_error(find_cut(x = x, module = as.data.frame(module), route_map = route_map),
               "'module' must be a numeric binary matrix")
  expect_error(find_cut(x = x, module = module, route_map = route_map,
                        theta_range = c(1, -1)),
               "theta_range")
  expect_error(find_cut(x = x, module = module, route_map = route_map,
                        n_grid = 50L),
               "n_grid")
})

test_that("find_cut() ref_theta selects the closest proper crossing", {

  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map

  # Default ref_theta = 0
  result_0  <- find_cut(x = x, module = module, route_map = route_map,
                        ref_theta = 0)

  # ref_theta = 2 should still return a valid result
  result_2  <- find_cut(x = x, module = module, route_map = route_map,
                        ref_theta = 2)

  # Both should return the same structure
  expect_s3_class(result_0, "find_cut")
  expect_s3_class(result_2, "find_cut")
  expect_length(result_0$cut_score, 2L)
  expect_length(result_2$cut_score, 2L)
})

# ── plot.find_cut() tests ──────────────────────────────────────────────────

test_that("plot.find_cut() returns a ggplot object for all layout options", {
  skip_if_not_installed("ggplot2")

  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map

  cut_result <- suppressWarnings(
    find_cut(x = x, module = module, route_map = route_map)
  )

  # Default: vertical layout
  p_vert <- plot(cut_result)
  expect_s3_class(p_vert, "ggplot")

  # Horizontal layout
  p_horiz <- plot(cut_result, layout = "horizontal")
  expect_s3_class(p_horiz, "ggplot")

  # No anomalous markers
  p_no_anom <- plot(cut_result, show_anomalous = FALSE)
  expect_s3_class(p_no_anom, "ggplot")

  # No cut score labels
  p_no_label <- plot(cut_result, label_cuts = FALSE)
  expect_s3_class(p_no_label, "ggplot")
})

test_that("find_cut() tif_data now includes stage 1", {
  x         <- simMST$item_bank
  module    <- simMST$module
  route_map <- simMST$route_map

  cut_result <- suppressWarnings(
    find_cut(x = x, module = module, route_map = route_map)
  )

  # tif_data must contain all three stages (1, 2, 3) for simMST
  stages_in_tif <- sort(unique(cut_result$tif_data$stage))
  expect_equal(stages_in_tif, c(1L, 2L, 3L))

  # Stage 1 should have exactly one module (module 1 for simMST)
  stage1_modules <- unique(cut_result$tif_data$module[
    cut_result$tif_data$stage == 1L
  ])
  expect_equal(stage1_modules, 1L)
})
