#' Draw Raw and Standardized Residual Plots
#'
#' This method provides graphical displays of the residuals between
#' observed data and model-based predictions (Hambleton et al., 1991).
#' For each score category of an item, it generates two types of residual plots:
#' (a) the raw residual plot and (b) the standardized residual plot.
#' Note that for dichotomous items, residual plots are drawn only for score
#' category 1.
#'
#' @inheritParams plot.info
#' @param x x An object of class `irtfit` obtained from [irtQ::irtfit()].
#' @param item.loc An integer specifying the position of the item to be plotted
#'   (i.e., the *n*th item in the item set). See **Details** below.
##' @param type A character string indicating the type of residual plot to be
##'  displayed. Available options are:
##'  - `"icc"` for the raw residual plot
##'  - `"sr"` for the standardized residual plot
#'   - `"both"` for displaying both plots
#'
#'   Default is `"both"`.
#' @param ci.method A character string specifying the method used to compute
#'   confidence intervals for the raw residual plot. Available options are:
#'   - `"wald"` for the Wald method
#'   - `"wilson"` for the Wilson score interval
#'   - `"wilson.cr"` for the Wilson interval with continuity correction
#'
#'   Default is `"wald"`. See **Details** below.
#' @param show.table A logical value indicating whether to return the
#'   contingency table used for drawing the residual plots of the specified
#'   item. If `TRUE`, the function returns the same contingency table stored in the
#'   internal `contingency.plot` object of the `irtfit` object. Default is `TRUE`.
#' @param layout.col An integer specifying the number of columns in the panel layout
#'   when plotting residuals for a polytomous item. Default is 2.
#' @param xlab.text A character string specifying the title for the x-axis.
#'   If omitted, a default label is used.
#' @param ylab.text A character string specifying the title for the y-axis.
#'   If `type = "both"`, a character vector of length two can be provided,
#'   corresponding to the y-axis titles for the raw residual and standardized
#'   residual plots, respectively. If omitted, default labels are used.
#' @param main.text A character string specifying the main title for the plot.
#'   If `type = "both"`, a character vector of length two can be provided,
#'   corresponding to the main titles for the raw residual and standardized
#'   residual plots, respectively. If omitted, default titles are used.
#' @param point.size A numeric value specifying the size of points. Default is 2.5.
#' @param strip.size A numeric value specifying the size of facet label text. Default is 12.
#' @param ylim.icc A numeric vector of length two specifying the y-axis limits
#'   for the raw residual plot. Default is c(0, 1).
#' @param ylim.sr.adjust Logical. If `TRUE`, the y-axis range for the
#'   standardized residual plot is automatically adjusted based on the maximum
#'   residual value for each item. If `FALSE`, the range is fixed according to
#'   the values specified in the `ylim.sr` argument. Default is `FALSE`.
#' @param ylim.sr A numeric vector of length two specifying the y-axis limits
#'   for the standardized residual plot. Default is c(-4, 4).
#' @param ... Additional arguments passed to [ggplot2::ggplot()] from the
#'   \pkg{ggplot2} package.
#'
#' @details All plots are generated using the \pkg{ggplot2} package.
#'
#'   Once the IRT model fit analysis is completed using [irtQ::irtfit()], the
#'   resulting object of class `irtfit` can be used to draw raw and standardized
#'   residual plots.These plots are primarily based on the information stored in
#'   the internal object `contingency.plot`.
#'
#'   Because residual plots are generated for one item at a time, you must
#'   specify which item to evaluate by providing an integer value for the
#'   `item.loc` argument, indicating the item's position in the test form.
#'   For example, to draw residual plots for the third item, set `item.loc = 3`.
#'
#'   For the raw residual plot, the `ci.method` argument determines the method
#'   used to estimate confidence intervals. The available methods are:
#'   - `"wald"`: Wald interval based on the normal approximation (Laplace, 1812)
#'   - `"wilson"`: Wilson score interval (Wilson, 1927)
#'   - `"wilson.cr"`: Wilson score interval with continuity correction (Newcombe, 1998)
#'
#'   For more information, see
#'   <https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval>.
#'   Note that the width of the confidence interval is governed by the
#'   \eqn{\alpha}-level specified in the `alpha` argument of the
#'   [irtQ::irtfit()] function.
#'
#'   For the standardized residual plot, residuals exceeding the threshold
#'   specified in the `overSR` argument of the [irtQ::irtfit()] function
#'   are displayed as circles. Residuals that do not exceed the threshold
#'   are displayed as crosses.
#'
#' @return This method displays the IRT raw residual plot, standardized
#'   residual plot, or both for the specified item. When `show.table = TRUE`, a
#'   contingency table used to generate the residual plots is also returned.
#'   See [irtQ::irtfit()] for more details about the contingency table.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtfit()]
#'
#' @references Hambleton, R. K., Swaminathan, H., & Rogers, H. J.
#' (1991).*Fundamentals of item response theory*. Newbury Park, CA: Sage.
#'
#' Laplace, P. S. (1820).*Theorie analytique des probabilites* (in French).
#' Courcier.
#'
#' Newcombe, R. G. (1998). Two-sided confidence intervals for the single
#' proportion: comparison of seven methods.
#' *Statistics in medicine, 17*(8), 857-872.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and
#' statistical inference.
#' *Journal of the American Statistical Association, 22*(158), 209-212.
#'
#' @examples
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select the first two dichotomous items and the last polytomous item
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # Generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(1000, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # Compute fit statistics
#' fit <- irtfit(
#'   x = x, score = score, data = data, group.method = "equal.freq",
#'   n.width = 11, loc.theta = "average", range.score = c(-4, 4), D = 1,
#'   alpha = 0.05, overSR = 1.5
#' )
#'
#' # Residual plots for the first item (dichotomous item)
#' plot(x = fit, item.loc = 1, type = "both", ci.method = "wald",
#'      show.table = TRUE, ylim.sr.adjust = TRUE)
#'
#' # Residual plots for the third item (polytomous item)
#' plot(x = fit, item.loc = 3, type = "both", ci.method = "wald",
#'      show.table = FALSE, ylim.sr.adjust = TRUE)
#'
#' # Raw residual plot for the third item (polytomous item)
#' plot(x = fit, item.loc = 3, type = "icc", ci.method = "wald",
#'      show.table = TRUE, ylim.sr.adjust = TRUE)
#'
#' # Standardized residual plot for the third item (polytomous item)
#' plot(x = fit, item.loc = 3, type = "sr", ci.method = "wald",
#'      show.table = TRUE, ylim.sr.adjust = TRUE)
#' }
#'
#' @import ggplot2 dplyr
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @importFrom gridExtra grid.arrange
#'
#' @export
plot.irtfit <- function(x,
                        item.loc = NULL,
                        type = "both",
                        ci.method = c("wald", "wilson", "wilson.cr"),
                        show.table = TRUE,
                        layout.col = 2,
                        xlab.text,
                        ylab.text,
                        main.text,
                        lab.size = 15,
                        main.size = 15,
                        axis.size = 15,
                        line.size = 1.0,
                        point.size = 2.5,
                        strip.size = 12,
                        ylim.icc = c(0, 1),
                        ylim.sr.adjust = FALSE,
                        ylim.sr = c(-4, 4),
                        ...) {

  ## -------------------------------------------------------------------------
  if (is.null(item.loc)) {
    stop("The location of item being plotted should be specified in 'item.loc'.", call. = FALSE)
  }

  if (!is.null(item.loc)) {
    if (length(item.loc) > 1) {
      stop("A length of 'item.loc' must be 1.", call. = FALSE)
    }
  }

  # prepare data
  # extract an meta information of an item in which residual plot is drawn
  item_meta <- x$item_df[item.loc, ]

  # obtain the number of score categories of the item
  score.cats <- item_meta$cats

  # extract a score range for the plot
  range.score <- x$ancillary$range.score

  # extract a scaling factor
  D <- x$ancillary$scale.D

  # extract an alpha level
  alpha <- x$ancillary$alpha

  # a criterion of standardized residuals
  overSR <- x$ancillary$overSR

  # extract contingency table for the item
  ctg.tb <- x$contingency.plot[[item.loc]]

  # subset of the observed probability table for score categories
  obs.prop <-
    dplyr::select(ctg.tb, theta = "point", dplyr::starts_with("obs.prop")) %>%
    dplyr::rename_all(.funs = function(k) gsub(pattern = "obs.prop.", replacement = "", x = k))

  # subset of the frequency table for score categories
  obs.freq <- dplyr::select(ctg.tb, dplyr::starts_with("obs.freq"))

  # extract total frequencies for the intervals
  total <- ctg.tb$total

  # extract the standard errors
  se <- dplyr::select(ctg.tb, dplyr::starts_with("se"))

  # extract standardize the raw residuals
  std.rsd <- dplyr::select(ctg.tb, dplyr::starts_with("std.rsd"), theta = "point")

  # find a z-score corresponding to significance level
  zscore <- stats::qnorm(1 - alpha)

  # a data.frame including the standardized residuals and and information
  # to see if the standardized residuals are greater than a specified SR criterion.
  resid_df <-
    std.rsd %>%
    reshape2::melt(id.vars = "theta", variable.name = "score", value.name = "std.rsd") %>%
    dplyr::mutate(lessSR = factor(ifelse(abs(.data$std.rsd) > overSR, 1, 0), levels = c(0, 1))) %>%
    dplyr::mutate_at(.vars = "score", ~ {
      gsub(pattern = "std.rsd.", replacement = "", x = paste0("Score: ", .x))
    })

  # when the item is a DRM item
  if (score.cats == 2) {
    resid_df <-
      resid_df %>%
      dplyr::filter(.data$score == "Score: 1")
  }

  ## -------------------------------------------------------------------------
  # observed data information for ICC plot
  # reorganize the "obs.prop" data.frame
  tb1 <- reshape2::melt(data = obs.prop, variable.name = "score", id.vars = "theta", value.name = "obs.prop")

  # compute the confidence interval
  ci.method <- tolower(ci.method)
  ci.method <- match.arg(ci.method)
  ci <- suppressWarnings(
    switch(ci.method,
      wald =
        purrr::map2(.x = obs.prop[, -1], .y = zscore * se, .f = function(k, j) data.frame(lower = k - j, upper = k + j)) %>%
          do.call(what = "rbind"),
      wilson =
        purrr::map(obs.freq, .f = function(k) {
          purrr::map2(.x = k, .y = total, .f = function(i, j) {
            stats::prop.test(i, j, alternative = "two.sided", conf.level = 1 - alpha, correct = FALSE)$conf.int
          }) %>%
            bind.fill(type = "rbind") %>%
            data.frame() %>%
            dplyr::rename("lower" = "X1", "upper" = "X2")
        }) %>%
          do.call(what = "rbind"),
      wilson.cr =
        purrr::map(obs.freq, .f = function(k) {
          purrr::map2(.x = k, .y = total, .f = function(i, j) {
            stats::prop.test(i, j, alternative = "two.sided", conf.level = 1 - alpha, correct = TRUE)$conf.int
          }) %>%
            bind.fill(type = "rbind") %>%
            data.frame() %>%
            dplyr::rename("lower" = "X1", "upper" = "X2")
        }) %>%
          do.call(what = "rbind")
    )
  )

  # restrict the confidence intervals from 0 to 1 when Wald statistic is used
  if (ci.method == "wald") {
    ci[ci < 0] <- 0
    ci[ci > 1] <- 1
  }

  # a data.frame including the observed proportion and the confidence interval
  obs_df <-
    cbind(tb1, ci) %>%
    dplyr::mutate_at(.vars = "score", ~ {
      gsub(pattern = "^*", replacement = "Score: ", x = .x)
    })

  # define theta nodes for x-axis
  theta <- seq(range.score[1], range.score[2], 0.1)

  # a data.frame including the ICC across all score categories
  icc_df <-
    data.frame(theta, traceline(x = item_meta, theta, D = D)$prob.cats[[1]]) %>%
    dplyr::rename_all(.funs = function(k) gsub(pattern = "resp.", replacement = "", x = k)) %>%
    reshape2::melt(id.vars = "theta", variable.name = "score", value.name = "icc") %>%
    dplyr::mutate_at(.vars = "score", ~ {
      gsub(pattern = "^*", replacement = "Score: ", x = .x)
    })

  # when the item is dichotomous
  if (score.cats == 2) {
    obs_df <-
      obs_df %>%
      dplyr::filter(.data$score == "Score: 1")
    icc_df <-
      icc_df %>%
      dplyr::filter(.data$score == "Score: 1")
  }

  ## -------------------------------------------------------------------------
  # draw plots
  type <- tolower(type)
  # (1) draw ICC plots
  if (type == "icc") {
    if (missing(xlab.text)) xlab.text <- expression(theta)
    if (missing(ylab.text)) ylab.text <- "Probability"
    if (missing(main.text)) main.text <- paste0("Raw Residuals: ", item_meta$id)

    p <-
      ggplot2::ggplot(data = icc_df, mapping = ggplot2::aes(x = .data$theta, y = .data$icc), ...) +
      ggplot2::geom_line(color = "blue", linewidth = line.size) +
      ggplot2::geom_point(
        data = obs_df, mapping = ggplot2::aes(x = .data$theta, y = .data$obs.prop),
        shape = 4, size = point.size, color = "red"
      ) +
      ggplot2::geom_segment(
        data = obs_df,
        mapping = ggplot2::aes(x = .data$theta, y = .data$upper, xend = .data$theta, yend = .data$lower),
        lineend = "square", size = 0.7
      ) +
      ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
      ggplot2::ylim(ylim.icc[1], ylim.icc[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol = layout.col) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = main.size),
        axis.title = ggplot2::element_text(size = lab.size),
        axis.text = ggplot2::element_text(size = axis.size)
      ) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = strip.size, face = "bold"))

    print(p)
  }

  ## -------------------------------------------------------------------------
  if (type == "sr") {
    # (2) draw standardized residual plots
    if (missing(xlab.text)) xlab.text <- expression(theta)
    if (missing(ylab.text)) ylab.text <- "Standardized Residual"
    if (missing(main.text)) main.text <- paste0("Standardized Residuals: ", item_meta$id)
    point.color <- c("blue", "red")
    point.shape <- c(4, 1)

    # find a maximum value of the absolute standardized residuals
    # use the maximum value as the ylim value
    if (ylim.sr.adjust) {
      y.max <- round(max(abs(resid_df$std.rsd)) + 1, 0)
      ylim.sr <- c(-y.max, y.max)
    }

    p <-
      ggplot2::ggplot(data = resid_df, mapping = ggplot2::aes(x = .data$theta, y = .data$std.rsd), ...) +
      ggplot2::geom_point(mapping = ggplot2::aes(color = .data$lessSR, shape = .data$lessSR), size = point.size) +
      ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
      ggplot2::ylim(ylim.sr[1], ylim.sr[2]) +
      ggplot2::xlim(range.score[1], range.score[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol = layout.col) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = main.size),
        axis.title = ggplot2::element_text(size = lab.size),
        axis.text = ggplot2::element_text(size = axis.size)
      ) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = strip.size, face = "bold")) +
      ggplot2::scale_shape_manual(values = point.shape) +
      ggplot2::scale_colour_manual(values = point.color) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = c(-overSR, overSR), linetype = "dashed")

    print(p)
  }

  ## -------------------------------------------------------------------------
  if (type == "both") {
    xlab.text <- expression(theta)
    if (missing(ylab.text)) ylab.text <- c("Probability", "Standardized Residual")
    if (missing(main.text)) {
      main.text <- c(
        paste0("Raw Residuals: ", item_meta$id),
        paste0("Standardized Residuals: ", item_meta$id)
      )
    }

    # (1) draw ICC plots
    p1 <-
      ggplot2::ggplot(data = icc_df, mapping = ggplot2::aes(x = .data$theta, y = .data$icc), ...) +
      ggplot2::geom_line(color = "blue", linewidth = line.size) +
      ggplot2::geom_point(
        data = obs_df, mapping = ggplot2::aes(x = .data$theta, y = .data$obs.prop),
        shape = 4, size = point.size, color = "red"
      ) +
      ggplot2::geom_segment(
        data = obs_df,
        mapping = ggplot2::aes(x = .data$theta, y = .data$upper, xend = .data$theta, yend = .data$lower),
        lineend = "square", size = 0.7
      ) +
      ggplot2::labs(title = main.text[1], x = xlab.text, y = ylab.text[1]) +
      ggplot2::ylim(ylim.icc[1], ylim.icc[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol = 1) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = main.size),
        axis.title = ggplot2::element_text(size = lab.size),
        axis.text = ggplot2::element_text(size = axis.size)
      ) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = strip.size, face = "bold"))

    # (2) draw standardized residual plots
    point.color <- c("blue", "red")
    point.shape <- c(4, 1)

    # find a maximum value of the absolute standardized residuals
    # use the maximum value as the ylim value
    if (ylim.sr.adjust) {
      y.max <- round(max(abs(resid_df$std.rsd)) + 1, 0)
      ylim.sr <- c(-y.max, y.max)
    }

    p2 <-
      ggplot2::ggplot(data = resid_df, mapping = ggplot2::aes(x = .data$theta, y = .data$std.rsd), ...) +
      ggplot2::geom_point(mapping = ggplot2::aes(color = .data$lessSR, shape = .data$lessSR), size = point.size) +
      ggplot2::labs(title = main.text[2], x = xlab.text, y = ylab.text[2]) +
      ggplot2::ylim(ylim.sr[1], ylim.sr[2]) +
      ggplot2::xlim(range.score[1], range.score[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol = 1) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = main.size),
        axis.title = ggplot2::element_text(size = lab.size),
        axis.text = ggplot2::element_text(size = axis.size)
      ) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = strip.size, face = "bold")) +
      ggplot2::scale_shape_manual(values = point.shape) +
      ggplot2::scale_colour_manual(values = point.color) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = c(-overSR, overSR), linetype = "dashed")

    # combine the two plots in a window
    gridExtra::grid.arrange(p1, p2, ncol = 2) # For multitple plots
  }

  # return the contingency table
  if (show.table) {
    return(ctg.tb)
  }
}
