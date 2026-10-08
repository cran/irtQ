# irtfit() computes observed proportions, residuals, and fit statistics for each item

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_mix <- bring.flexmirt(file = prm_file, "par")$Group1$full_df[c(1:15, 51:55), ]
set.seed(2026)
theta_fit <- rnorm(1000)
resp_fit <- simdat(x = x_mix, theta = theta_fit, D = 1)

test_that("irtfit() observed proportions are the category frequencies divided by the group size", {
  for (gm in c("equal.width", "equal.freq")) {
    fit <- irtfit(
      x = x_mix, score = theta_fit, data = resp_fit, group.method = gm,
      n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1
    )
    for (i in seq_along(fit$contingency.plot)) {
      tb <- fit$contingency.plot[[i]]
      freq <- tb[, grep("^obs\\.freq\\.[0-9]+$", names(tb)), drop = FALSE]
      prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb)), drop = FALSE]
      # proportions in each interval sum to one
      expect_equal(unname(rowSums(prop)), rep(1, nrow(tb)))
      # proportions equal the frequencies over the interval totals
      expect_equal(unname(as.matrix(prop)), unname(as.matrix(freq)) / tb$total)
    }
  }
})

test_that("irtfit() raw residuals equal observed proportions minus expected probabilities", {
  fit <- irtfit(
    x = x_mix, score = theta_fit, data = resp_fit, group.method = "equal.freq",
    n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1
  )
  tb <- fit$contingency.plot[[1]]
  prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb))]
  expp <- tb[, grep("^exp\\.prob\\.[0-9]+$", names(tb))]
  rsd <- tb[, grep("^raw\\.rsd\\.[0-9]+$", names(tb))]
  expect_equal(unname(as.matrix(rsd)), unname(as.matrix(prop - expp)))
})

test_that("irtfit() drops empty score groups and keeps proportions that sum to one", {
  # two extreme scores leave empty intervals in the equal-width grouping
  sc <- theta_fit
  sc[1] <- -6
  sc[2] <- 6
  fit <- irtfit(
    x = x_mix, score = sc, data = resp_fit, group.method = "equal.width",
    n.width = 10, loc.theta = "average", range.score = c(-7, 7), D = 1
  )
  n_group <- vapply(fit$contingency.plot, nrow, integer(1))
  expect_true(any(n_group < 10L))
  for (tb in fit$contingency.plot) {
    prop <- tb[, grep("^obs\\.prop\\.[0-9]+$", names(tb))]
    expect_equal(unname(rowSums(prop)), rep(1, nrow(tb)))
  }
})

test_that("plot.irtfit() Wald intervals use the two-sided critical value", {
  # the last plot is read back, which needs get_last_plot() in ggplot2
  skip_if_not("get_last_plot" %in% getNamespaceExports("ggplot2"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fit <- irtfit(
    x = x_mix, score = theta_fit, data = resp_fit, group.method = "equal.freq",
    n.width = 8, loc.theta = "average", range.score = c(-4, 4), D = 1, alpha = 0.05
  )
  plot(x = fit, item.loc = 1, type = "icc", ci.method = "wald", show.table = FALSE)
  built <- ggplot2::ggplot_build(ggplot2::get_last_plot())
  # the segment layer holds the interval (y = upper limit, yend = lower limit)
  seg <- Filter(function(d) all(c("y", "yend", "xend") %in% names(d)), built$data)[[1]]
  unclipped <- seg$yend > 0 & seg$y < 1
  half_width <- (seg$y - seg$yend)[unclipped] / 2
  tb <- fit$contingency.plot[[1]]
  se_all <- c(tb$se.0, tb$se.1)
  # each half width equals the 97.5th percentile of the normal times a standard error
  expect_true(length(half_width) > 0)
  expect_true(all(vapply(half_width, function(h) any(abs(h / stats::qnorm(0.975) - se_all) < 1e-6), logical(1))))
})
