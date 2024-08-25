#' Plot ICC and TCC
#'
#' @description This method function plots item or test characteristic curve using the ggplot2 package. The item characteristic
#' (or category) curve (ICC) or item score curve is drawn for an individual item. The test characteristic curve (TCC) is drawn
#' based on a total test form.
#'
#' @param x An object of class \code{\link{traceline}}.
#' @param item.loc A numeric value indicating that the \emph{n}th item (or the location of item) is plotted.
#' If NULL, the TCC based on a total test form is drawn. Default is NULL.
#' @param score.curve Logical value. If TRUE, item score curve (i.e., a weighted sum of item category probabilities over the item scores) is plotted
#' in a panel. Otherwise, ICCs for all score categories are plotted in separate panels. For a dichotomous item, the item score curve is the same as
#' the ICC of score category 1. Ignored when \code{item.loc = NULL}. Default is FALSE.
#' @param overlap Logical value indicating whether multiple item score curves are plotted in one panel.
#' If FALSE, the multiple item score curves are displayed with multiple panels, one for each.
#' @param layout.col An integer value indicating the number of columns in the panel when displaying ICCs for an item or
#' when displaying multiple item scores with multiple panels.
#' @param xlab.text,ylab.text A title for the x and y axes.
#' @param main.text An overall title for the plot.
#' @param lab.size The size of xlab and ylab. Default is 15.
#' @param main.size The size of \code{main.text}. Default is 15.
#' @param axis.size The size of labels along the x and y axes. Default is 15.
#' @param line.color A character string specifying the color for a line. See \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/} for more details
#' about colors used in ggplot2.
#' @param line.size The size of lines. Default is 1.
#' @param strip.size The size of facet labels when ICCs for an item are plotted.
#' @param ... Further arguments passed from the function \code{geom_line()} in the \pkg{ggplot2} package.
#'
#' @details All of the plots are drawn using the ggplot2 package.
#' If \code{item.loc = NULL}, the TCC based on the total test form is plotted. In the argument \code{item.loc},
#' a vector of positive integer values should be specified to indicate the \emph{n}th items among the total test form. For example,
#' if there are ten items in the test form and the score curves of the 1st, 2nd, and 3rd items should be plotted, then \code{item.loc = 1:3}.
#'
#' @return This method function displays ICC or TCC plots of the studied item(s).
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{traceline}}
#'
#' @examples
#' ## example
#' ## using a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # read item parameters and transform them to item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-3, 3, 0.1)
#'
#' # compute the item category probabilities and item/test
#' # characteristic functions given the theta values
#' x <- traceline(x = test_flex, theta, D = 1)
#'
#' # plot TCC based on the total test form
#' plot(x, item.loc = NULL)
#'
#' # plot ICCs for the first item (dichotomous item)
#' plot(x, item.loc = 1, score.curve = FALSE, layout.col = 2)
#'
#' # plot item score curve for the first item (dichotomous item)
#' plot(x, item.loc = 1, score.curve = TRUE)
#'
#' # plot item score curves for the first six dichotomous items
#' # with multiple panels
#' plot(x, item.loc = 1:6, score.curve = TRUE, overlap = FALSE)
#'
#' # plot item score curve for the first six dichotomous items
#' # in one panel
#' plot(x, item.loc = 1:6, score.curve = TRUE, overlap = TRUE)
#'
#' # plot ICCs for the last item (polytomous item)
#' plot(x, item.loc = 55, score.curve = FALSE, layout.col = 2)
#'
#' # plot item score curve for the last item (polytomous item)
#' plot(x, item.loc = 55, score.curve = TRUE)
#'
#' # plot item score curves for the last three polytomous items
#' # with multiple panels
#' plot(x, item.loc = 53:55, score.curve = TRUE, overlap = FALSE)
#'
#' # plot item score curves for the last three polytomous items
#' # in one panel
#' plot(x, item.loc = 53:55, score.curve = TRUE, overlap = TRUE)
#'
#' @import ggplot2 dplyr
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
plot.traceline <- function(x, item.loc = NULL,
                           score.curve = FALSE, overlap = FALSE, layout.col = 2,
                           xlab.text, ylab.text, main.text, lab.size = 15, main.size = 15, axis.size = 15,
                           line.color, line.size = 1, strip.size = 12, ...) {
  ## ----------------------------------------------------------
  if (length(item.loc) > 1 & score.curve == FALSE) {
    stop("To plot score curves for multiple items, set 'score.curve = TRUE'.", call. = FALSE)
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
      if (missing(main.text)) main.text <- paste0("Item Characteristic Curve: ", names(x$prob.cat[item.loc]))
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
          ggplot2::geom_line(linewidth = line.size, color = line.color) +
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
          ggplot2::geom_line(mapping = ggplot2::aes(color = .data$Item), linewidth = line.size) +
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
