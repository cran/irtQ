#' Plot Item and Test Characteristic Curves
#'
#' This method visualizes item characteristic curves (ICCs), item score curves,
#' or the test characteristic curve (TCC) using the \pkg{ggplot2} package.
#' ICCs or item score curves can be plotted for one or more selected items,
#' while the TCC is plotted for the entire test form.
#'
#' @inheritParams plot.info
#' @param x x An object of class `traceline` obtained from [irtQ::traceline()].
#' @param item.loc A numeric vector specifying the position(s) of the item(s) to
#'   plot. If `NULL` (default), the test characteristic curve (TCC) for the
#'   entire test form is plotted.
#' @param score.curve Logical. If `TRUE`, plots the item score curve, defined
#'   as the weighted sum of category probabilities across score categories, in
#'   a panel.
#'
#'   If `FALSE`, plots item characteristic curves (ICCs) for all score categories,
#'   either in separate panels or in a single panel depending on the `overlap`
#'   setting.
#'
#'   For dichotomous items, the item score curve is equivalent to the ICC for
#'   score category 1. Ignored when `item.loc = NULL`. Default is `FALSE`.
#' @param overlap Logical. Determines how multiple curves are displayed when
#'   plotting ICCs or item score curves.
#'
#'   If `TRUE`, curves are overlaid in a single panel using different colors.
#'   If `FALSE`, each curve is drawn in a separate panel—either one panel per
#'   item or per score category, depending on the setting of `score.curve`.
#' @param layout.col An integer value indicating the number of columns in the
#'   plot when displaying multiple panels. Used only when `overlap = FALSE`.
#'   Default is 2.
#' @param strip.size Numeric. Font size of facet labels when ICCs are plotted.
#'
#' @details
#' All plots are generated using the \pkg{ggplot2} package.
#' If `item.loc = NULL`, the test characteristic curve (TCC) for the entire test
#' form is plotted. If `item.loc` is specified, it should be a vector of
#' positive integers indicating the position(s) of the items to be plotted.
#' For example, if the test form includes ten items and you wish to plot the
#' score curves of the 1st, 2nd, and 3rd items, set `item.loc = 1:3`.
#'
#' @return This method displays item characteristic curves (ICCs), item score
#' curves, or the test characteristic curve (TCC), depending on the specified
#' arguments.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::traceline()]
#'
#' @examples
#' ## Example using a "-prm.txt" file exported from flexMIRT
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read the item parameters and convert them to item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Define a sequence of theta values
#' theta <- seq(-3, 3, 0.1)
#'
#' # Compute item category probabilities and item/test characteristic functions
#' x <- traceline(x = test_flex, theta, D = 1)
#'
#' # Plot the test characteristic curve (TCC) for the full test form
#' plot(x, item.loc = NULL)
#'
#' # Plot ICCs for the first item (dichotomous),
#' # with a separate panel for each score category
#' plot(x, item.loc = 1, score.curve = FALSE, layout.col = 2)
#'
#' # Plot ICCs for the first item in a single panel
#' # (all score categories overlaid)
#' plot(x, item.loc = 1, score.curve = FALSE, overlap = TRUE)
#'
#' # Plot ICCs for multiple items (both dichotomous and polytomous),
#' # with each item's ICCs shown in a single panel
#' plot(x, item.loc = c(1:3, 53:55), score.curve = FALSE, overlap = TRUE)
#'
#' # Plot the item score curve for the first item (dichotomous)
#' plot(x, item.loc = 1, score.curve = TRUE)
#'
#' # Plot item score curves for the first six dichotomous items
#' # using multiple panels
#' plot(x, item.loc = 1:6, score.curve = TRUE, overlap = FALSE)
#'
#' # Plot item score curves for the first six dichotomous items
#' # overlaid in a single panel
#' plot(x, item.loc = 1:6, score.curve = TRUE, overlap = TRUE)
#'
#' # Plot ICCs for the last item (polytomous),
#' # with each score category in a separate panel
#' plot(x, item.loc = 55, score.curve = FALSE, layout.col = 2)
#'
#' # Plot the item score curve for the last item (polytomous)
#' plot(x, item.loc = 55, score.curve = TRUE)
#'
#' # Plot item score curves for the last three polytomous items
#' # using multiple panels
#' plot(x, item.loc = 53:55, score.curve = TRUE, overlap = FALSE)
#'
#' # Plot item score curves for the last three polytomous items
#' # overlaid in a single panel
#' plot(x, item.loc = 53:55, score.curve = TRUE, overlap = TRUE)
#'
#' @import ggplot2 dplyr
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
plot.traceline <- function(x,
                           item.loc = NULL,
                           score.curve = FALSE,
                           overlap = FALSE,
                           layout.col = 2,
                           xlab.text,
                           ylab.text,
                           main.text,
                           lab.size = 15,
                           main.size = 15,
                           axis.size = 15,
                           line.color,
                           line.size = 1,
                           strip.size = 12,
                           ...) {

  ## ----------------------------------------------------------
  if (length(item.loc) > 1 & score.curve == FALSE & overlap == FALSE) {
    stop(paste("When 'score.curve = FALSE' and 'overlap = FALSE',",
               "    item characteristic curves can be plotted for a single item only.",
               "    Try to set 'overlap = TRUE'.",
               sep = "\n"), call. = FALSE)
  }

  # extract theta values for x-axis
  theta <- x$theta

  # 1. plot TCCs
  if (is.null(item.loc)) {
    # data manipulation for plotting
    tcc.trace <- x$tcc
    df_tcc <- data.frame(tcc = tcc.trace, theta = theta)

    # plot
    # Set plot conditions
    if (missing(xlab.text)) xlab.text <- expression(theta)
    if (missing(ylab.text)) ylab.text <- "Expected Score"
    if (missing(main.text)) main.text <- "Test Characteristic Curve"
    if (missing(line.color)) line.color <- "#F8766D" else line.color <- line.color
    max.score <-
      purrr::map_dbl(x$prob.cat, .f = function(k) ncol(k) - 1) %>%
      sum()

    # draw a plot
    p <-
      df_tcc %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$tcc)) +
      ggplot2::geom_line(linewidth = line.size, color = line.color, ...) +
      ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
      ggplot2::ylim(0, max.score) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = main.size),
        axis.title = ggplot2::element_text(size = lab.size),
        axis.text = ggplot2::element_text(size = axis.size)
      )
  }

  # 2. plot ICCs
  if (!is.null(item.loc)) {
    if (!score.curve) {

      if (overlap) {

        # a data.frame including the ICC across all score categories for multiple items
        icc_all <-
          purrr::map(.x = item.loc,
                     .f = function(i) {
                       icc_df <-
                         data.frame(theta = theta, x$prob.cat[[i]]) %>%
                         dplyr::rename_all(.funs = function(k) gsub(pattern = "score.", replacement = "", x = k)) %>%
                         reshape2::melt(id.vars = "theta", variable.name = "score", value.name = "icc")
                       icc_df$score <- gsub(pattern = "resp.", replacement = "", x = icc_df$score)
                       icc_df$Item <- names(x$prob.cat[i])
                       icc_df
                     }) %>%
          dplyr::bind_rows()

        ## -------------------------------------------------------------------------
        # draw ICC plots
        if (missing(xlab.text)) xlab.text <- expression(theta)
        if (missing(ylab.text)) ylab.text <- "Probability"
        if (length(item.loc) == 1) {
          if (missing(main.text)) main.text <- paste0("Item Characteristic Curve: ", names(x$prob.cat[item.loc]))
        } else if (length(item.loc) > 1) {
          if (missing(main.text)) main.text <- "Item Characteristic Curve"
        }
        if (missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

        p <-
          icc_all %>%
          dplyr::mutate_at(.vars = "Item",
                           .funs = ~{factor(x = .x, levels = names(x$prob.cat[item.loc]))}) %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$icc)) +
          ggplot2::geom_line(mapping = ggplot2::aes(color = .data$score),
                             linewidth = line.size, ...) +
          ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text,
                        color = "Category") +
          ggplot2::ylim(0, 1) +
          ggplot2::theme_bw() +
          ggplot2::facet_wrap(~Item, ncol = layout.col) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = main.size),
            axis.title = ggplot2::element_text(size = lab.size),
            axis.text = ggplot2::element_text(size = axis.size)
          ) +
          ggplot2::theme(strip.text.x =
                           ggplot2::element_text(size = strip.size, face = "bold"))


      } else {

        # a data.frame including the ICC across all score categories
        icc_df <-
          data.frame(theta = theta, x$prob.cat[[item.loc]]) %>%
          dplyr::rename_all(.funs = function(k) gsub(pattern = "score.", replacement = "", x = k)) %>%
          reshape2::melt(id.vars = "theta", variable.name = "score", value.name = "icc")
        icc_df$score <- gsub(pattern = "^*", replacement = "Score: ", x = icc_df$score)

        ## -------------------------------------------------------------------------
        # draw ICC plots
        if (missing(xlab.text)) xlab.text <- expression(theta)
        if (missing(ylab.text)) ylab.text <- "Probability"
        if (missing(main.text)) main.text <- paste0("Item Characteristic Curve: ",
                                                    names(x$prob.cat[item.loc]))
        if (missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

        p <-
          ggplot2::ggplot(data = icc_df, mapping = ggplot2::aes(x = .data$theta, y = .data$icc)) +
          ggplot2::geom_line(color = line.color, linewidth = line.size, ...) +
          ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
          ggplot2::ylim(0, 1) +
          ggplot2::theme_bw() +
          ggplot2::facet_wrap(~score, ncol = layout.col) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = main.size),
            axis.title = ggplot2::element_text(size = lab.size),
            axis.text = ggplot2::element_text(size = axis.size)
          ) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size = strip.size, face = "bold"))

      }

    }

    if (score.curve) {
      # check the number of score categories
      # cats <- ncol(x$prob.cat[[item.loc]])
      cats <- purrr::map_dbl(.x = x$prob.cat[item.loc], ncol)

      # data manipulation for plotting
      score.trace <- x$icc[, item.loc, drop = FALSE]
      df_score <-
        data.frame(score.trace, theta = theta) %>%
        reshape2::melt(variable.name = "item", id.vars = "theta", value.name = "icc")

      # data manipulation for plotting
      # df_info$item <- as.numeric(df_info$item)

      # plot
      # Set plot conditions
      if (missing(xlab.text)) xlab.text <- expression(theta)
      if (missing(ylab.text)) ylab.text <- "Expected Score"
      if (length(cats) == 1) {
        if (missing(main.text)) main.text <- paste0("Item Score Curve: ", names(x$prob.cat[item.loc]))
      } else if (length(cats) > 1) {
        if (missing(main.text)) main.text <- "Item Score Curve"
      }
      if (missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

      # draw a plot
      if (!overlap) {
        p <-
          df_score %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$icc)) +
          ggplot2::geom_line(linewidth = line.size, color = line.color, ...) +
          ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
          ggplot2::ylim(0, (max(cats) - 1)) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = main.size),
            axis.title = ggplot2::element_text(size = lab.size),
            axis.text = ggplot2::element_text(size = axis.size)
          ) +
          ggplot2::facet_wrap(~item, ncol = layout.col) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size = strip.size, face = "bold"))
      } else {
        p <-
          df_score %>%
          dplyr::rename("Item" = "item") %>%
          dplyr::mutate_at(.vars = "Item", as.factor) %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$icc)) +
          ggplot2::geom_line(mapping = ggplot2::aes(color = .data$Item),
                             linewidth = line.size, ...) +
          ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
          ggplot2::ylim(0, (max(cats) - 1)) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = main.size),
            axis.title = ggplot2::element_text(size = lab.size),
            axis.text = ggplot2::element_text(size = axis.size)
          )
      }
    }
  }

  print(p)
}
