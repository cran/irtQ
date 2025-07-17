#' Plot Item and Test Information Functions
#'
#' This method plots item or test information functions for a specified set of
#' theta values. It can also display the conditional standard error of
#' estimation (CSEE) at the test level.
#'
#' @param x x An object of class `info` obtained from [irtQ::info()].
#' @param item.loc A numeric vector indicating which item information functions to plot,
#'   specified by item position (e.g., 1 for the first item). If `NULL` (default),
#'   the test information function for the entire test form is plotted.
#' @param overlap Logical. If `TRUE`, multiple item information functions are
#'   plotted in a single panel. If `FALSE` (default), each item information
#'   function is displayed in a separate panel.
#' @param csee Logical. If `TRUE`, plots the conditional standard error of
#'   estimation (CSEE) at the test level. Note that the CSEE plot is only
#'   available at the test level, not for individual items. If `FALSE`
#'   (default), item or test information functions are plotted.
#' @param xlab.text,ylab.text Character strings specifying the labels for
#'   the x and y axes, respectively.
#' @param main.text Character string specifying the overall title of the plot.
#' @param lab.size Numeric value specifying the font size of axis titles.
#'   Default is 15.
#' @param main.size Numeric value specifying the font size of the plot title.
#'   Default is 15.
#' @param axis.size Numeric value specifying the font size of axis tick labels.
#'   Default is 15.
#' @param line.color A character string specifying the color of the plot lines.
#'   See <http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/> for available
#'   color names.
#' @param line.size Numeric value specifying the thickness of plot lines.
#'   Default is 1.
#' @param layout.col Integer. Number of columns to use when faceting multiple
#'   item information functions. Used only when `overlap = FALSE`. Default is 4.
#' @param strip.size Numeric value specifying the font size of facet labels
#'   when multiple items are displayed.
#' @param ... Additional arguments passed to [ggplot2::geom_line()] from
#'   the \pkg{ggplot2} package.
#'
#' @details All of the plots are drawn using the \pkg{ggplot2} package. The
#'   object of class `info` can be obtained from the function [irtQ::info()].
#'
#' @return This method function displays the item or test information function
#'   plot. When `csee = TRUE`, the CSEE is returned at the test level.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::info()]
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
#' theta <- seq(-4, 4, 0.1)
#'
#' # Compute item and test information values for the given theta values
#' x <- info(x = test_flex, theta = theta, D = 1, tif = TRUE)
#'
#' # Plot the test information function
#' plot(x)
#'
#' # Plot the item information function for the second item
#' plot(x, item.loc = 2)
#'
#' # Plot multiple item information functions, each in a separate panel
#' plot(x, item.loc = 1:8, overlap = FALSE)
#'
#' # Plot multiple item information functions in a single panel
#' plot(x, item.loc = 1:8, overlap = TRUE)
#'
#' # Plot the conditional standard error of estimation (CSEE) at the test level
#' plot(x, csee = TRUE)
#'
#' @import ggplot2 dplyr
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
plot.info <- function(x,
                      item.loc = NULL,
                      overlap = FALSE,
                      csee = FALSE,
                      xlab.text,
                      ylab.text,
                      main.text,
                      lab.size = 15,
                      main.size = 15,
                      axis.size = 15,
                      line.color,
                      line.size = 1,
                      layout.col = 4,
                      strip.size = 12,
                      ...) {

  if (!csee) {
    # 1. plot test infomation
    if (is.null(item.loc)) {
      # data manipulation for plotting
      if (is.null(x$tif)) {
        stop("The test information (tif) is NULL in the provided object 'x'. ", call. = FALSE)
      }
      df_info <- data.frame(theta = x$theta, info = x$tif)

      # plot
      # Set plot conditions
      if (missing(xlab.text)) xlab.text <- expression(theta)
      if (missing(ylab.text)) ylab.text <- "Information"
      if (missing(main.text)) main.text <- "Test Information"
      if (missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

      # draw a plot
      p <-
        df_info %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$info)) +
        ggplot2::geom_line(linewidth = line.size, color = line.color, ...) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = main.size),
          axis.title = ggplot2::element_text(size = lab.size),
          axis.text = ggplot2::element_text(size = axis.size)
        )
    }

    # 2. plot item information
    if (!is.null(item.loc)) {
      # data manipulation for plotting
      df_info <-
        data.frame(t(x$iif[item.loc, , drop = FALSE]), theta = x$theta) %>%
        reshape2::melt(variable.name = "item", id.vars = "theta", value.name = "info")

      # plot
      # Set plot conditions
      if (missing(xlab.text)) xlab.text <- expression(theta)
      if (missing(ylab.text)) ylab.text <- "Information"
      if (length(item.loc) == 1) {
        if (missing(main.text)) main.text <- paste0("Item Information: ", unique(df_info$item))
      } else if (length(item.loc) > 1) {
        if (missing(main.text)) main.text <- "Item Information"
      }
      if (missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

      if (!overlap) {
        p <-
          df_info %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$info)) +
          ggplot2::geom_line(linewidth = line.size, color = line.color, ...) +
          ggplot2::theme_bw() +
          ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = main.size),
            axis.title = ggplot2::element_text(size = lab.size),
            axis.text = ggplot2::element_text(size = axis.size)
          ) +
          ggplot2::facet_wrap(~item, ncol = layout.col) +
          ggplot2::theme(strip.text.x = ggplot2::element_text(size = strip.size, face = "bold"))
      } else {
        p <-
          df_info %>%
          dplyr::rename("Item" = "item") %>%
          dplyr::mutate_at(.vars = "Item", as.factor) %>%
          ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$info)) +
          ggplot2::geom_line(mapping = ggplot2::aes(color = .data$Item), linewidth = line.size, ...) +
          ggplot2::theme_bw() +
          ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = main.size),
            axis.title = ggplot2::element_text(size = lab.size),
            axis.text = ggplot2::element_text(size = axis.size)
          )
      }
    }
  } else {
    # Plot only test level csee infomation
    # data manipulation for plotting
    df_csee <- data.frame(theta = x$theta, csee = 1 / sqrt(x$tif))

    # plot
    # Set plot conditions
    if (missing(xlab.text)) xlab.text <- expression(theta)
    if (missing(ylab.text)) ylab.text <- "Standard Error"
    if (missing(main.text)) main.text <- "Conditional Standard Error of Estimation"
    if (missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

    # draw a plot
    p <-
      df_csee %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$theta, y = .data$csee)) +
      ggplot2::geom_line(linewidth = line.size, color = line.color, ...) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = main.size),
        axis.title = ggplot2::element_text(size = lab.size),
        axis.text = ggplot2::element_text(size = axis.size)
      )
  }

  p
}
