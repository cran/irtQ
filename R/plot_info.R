#' Plot Item and Test Information Functions
#'
#' @description This method function plots item or test information function given a specified theta values. In addition,
#' it displays the conditional standard errors at a test level.
#'
#' @param x x An object of class \code{\link{info}}.
#' @param item.loc A vector of numeric values indicating that the item information functions of the \emph{n}th items
#' (or the location of items in a test form) are plotted. If NULL, the test information function for the total test form is drawn.
#' Default is NULL.
#' @param overlap Logical value indicating whether multiple item information functions are plotted in one panel.
#' If FALSE, multiple item information functions are displayed in multiple panels, one for each.
#' @param csee Logical value indicating whether the function displays the conditional standard error of estimation (CSEE) at a test level.
#' If FALSE, item/test information function is plotted. Note that the CSEE plot is displayed only at a test level.
#' @param xlab.text,ylab.text A title for the x and y axes.
#' @param main.text An overall title for the plot.
#' @param lab.size The size of xlab and ylab. Default is 15.
#' @param main.size The size of \code{main.text}. Default is 15.
#' @param axis.size The size of labels along the x and y axes. Default is 15.
#' @param line.color A character string specifying a color for the line. See \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}
#' for more details about colors used in ggplot2.
#' @param line.size The size of lines. Default is 1.
#' @param layout.col An integer value indicating the number of columns in the panel when displaying the item information functions of
#' the multiple items. Default is 4.
#' @param strip.size The size of facet labels when the item information functions of the multiple items are drawn.
#' @param ... Further arguments passed from the function \code{geom_line()} in the \pkg{ggplot2} package.
#'
#' @details All of the plots are drawn using the ggplot2 package.
#' The object of class \code{\link{info}} can be obtained from the function \code{\link{info}}.
#'
#' @return This method function displays the item or test information function plot. When \code{csee = TRUE},
#' the conditional standard error is returned at the test level.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{info}}
#'
#' @examples
#' ## the use of a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # read item parameters and transform them to item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-4, 4, 0.1)
#'
#' # compute item and test information values given the theta values
#' x <- info(x = test_flex, theta = theta, D = 1, tif = TRUE)
#'
#' # draw a plot of the test information function
#' plot(x)
#'
#' # draw a plot of the item information function for the second item
#' plot(x, item.loc = 2)
#'
#' # draw a plot of multiple item information functions across the multiple panels
#' plot(x, item.loc = 1:8, overlap = FALSE)
#'
#' # draw a plot of multiple item information functions in one panel
#' plot(x, item.loc = 1:8, overlap = TRUE)
#'
#' # draw a plot of conditional standard error at a test level
#' plot(x, csee = TRUE)
#'
#' @import ggplot2 dplyr
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
plot.info <- function(x, item.loc = NULL, overlap = FALSE, csee = FALSE, xlab.text, ylab.text, main.text,
                      lab.size = 15, main.size = 15, axis.size = 15, line.color, line.size = 1, layout.col = 4,
                      strip.size = 12, ...) {
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
