#' Log-Likelihood of Ability Parameters
#'
#' This function computes the log-likelihood values for a set of ability
#' parameters, given item parameters and response data
#'
#' @inheritParams catsib
#' @inheritParams est_score
#' @param theta A numeric vector of ability values at which to evaluate the
#'   log-likelihood function.
#' @param method A character string specifying the estimation method. Available
#'   options include:
#'   - `"ML"`: Maximum Likelihood Estimation
#'   - `"MLF"`: Maximum Likelihood Estimation with Fences
#'   - `"MAP"`: Maximum A Posteriori Estimation
#'   Default is `"ML"`.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal prior distribution (used only when `method
#'   = "MAP"`). Default is `c(0, 1)`. Ignored for `"ML"` and `"MLF"`.
#'
#' @details This function evaluates the log-likelihood of a given ability
#' (`theta`) for one or more examinees, based on item parameters (`x`) and item
#' response data (`data`).
#'
#' If `method = "MLF"` is selected, the function appends two virtual "fence"
#' items to the item pool with fixed parameters. These artificial items help
#' avoid unstable likelihood functions near the boundaries of the ability scale.
#'
#' For example, to compute the log-likelihood curves of two examinees' responses
#' to the same test items, supply a 2-row matrix to `data` and a vector of
#' ability values to `theta`.
#'
#' @return A data frame of log-likelihood values.
#' - Each **row** corresponds to an ability value (`theta`).
#' - Each **column** corresponds to an examinee’s response pattern.
#'
#' @examples
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameters and convert them to item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Generate ability values from N(0, 1)
#' set.seed(10)
#' score <- rnorm(5, mean = 0, sd = 1)
#'
#' # Simulate response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' # Specify ability values for log-likelihood evaluation
#' theta <- seq(-3, 3, 0.5)
#'
#' # Compute log-likelihood values (using MLE)
#' llike_score(x = x, data = data, theta = theta, D = 1, method = "ML")
#'
#' @export
#' @importFrom reshape2 melt
#' @import dplyr
llike_score <- function(x,
                        data,
                        theta,
                        D = 1,
                        method = "ML",
                        norm.prior = c(0, 1),
                        fence.a = 3.0,
                        fence.b = NULL,
                        missing = NA) {

  # check if the data set is a vector of an examinee
  if (is.vector(data)) {
    data <- rbind(data)
  }

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of examinees
  nstd <- nrow(data)

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # add two more items and data responses when "ML" with Fences method is used
  if (method == "MLF") {
    # when fence.b = NULL, use the range argument as the fence.b argument
    if (is.null(fence.b)) {
      fence.b <- range
    }

    # add two more response columns for the two fence items
    data <- cbind(data, f.lower = 1, f.upper = 0)

    # create a new item metadata for the two fence items
    x.fence <- shape_df(
      par.drm = list(a = rep(fence.a, 2), b = fence.b, g = rep(0, 2)),
      item.id = c("fence.lower", "fence.upper"), cats = 2,
      model = "3PLM"
    )

    # create the new item metadata by adding two fence items
    x <- dplyr::bind_rows(x, x.fence)
  }

  # check the maximum score category across all items
  max.cats <- max(x$cats)

  # count the number of columns of the item metadata
  max.col <- ncol(x)

  # create an empty elm_item object
  elm_item <- list(pars = NULL, model = NULL, cats = NULL)

  # reshape the data
  data <-
    data.frame(x, t(data), check.names = FALSE) %>%
    reshape2::melt(
      id.vars = 1:max.col,
      variable.name = "std",
      value.name = "resp",
      na.rm = TRUE,
      factorsAsStrings = FALSE
    )
  data$resp <- factor(data$resp, levels = (seq_len(max.cats) - 1))
  data$std <- as.numeric(data$std)

  # create a score table to contain the scoring results
  rst <- data.frame(std = 1:nstd)

  # compute the log-likelihood values for all the discrete theta values
  lls <-
    purrr::map(
      .x = 1:nstd,
      .f = function(x) {
        exam_dat <-
          data %>%
          dplyr::filter(.data$std == x)
        llike_score_one(
          exam_dat = exam_dat, theta = theta,
          elm_item = elm_item, max.col = max.col, D = D,
          method = method, norm.prior = norm.prior
        )
      }
    ) %>%
    bind.fill(type = "cbind") %>%
    data.frame() %>%
    dplyr::rename_all(.f = ~ {
      paste0("Resp.", 1:nstd)
    })

  # return results
  lls
}


# This function computes the loglikelihood values for one test data set
llike_score_one <- function(exam_dat, theta, elm_item, max.col, D = 1, method = "ML",
                            norm.prior = c(0, 1)) {
  # extract the required objects from the individual exam data
  elm_item$pars <- data.matrix(exam_dat[, 4:max.col])
  elm_item$model <- exam_dat$model
  elm_item$cats <- exam_dat$cats
  resp <- exam_dat$resp
  n.resp <- length(resp)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # calculate the score categories
  tmp.id <- 1:n.resp
  freq.cat <-
    matrix(
      stats::xtabs(~ tmp.id + resp,
        na.action = stats::na.pass, addNA = FALSE
      ),
      nrow = n.resp
    )

  # compute the negative log-likelihood values for all the discrete theta values
  ll_val <- ll_score(
    theta = theta, elm_item = elm_item, freq.cat = freq.cat, method = method,
    idx.drm = idx.drm, idx.prm = idx.prm, D = D, norm.prior = norm.prior,
    logL = TRUE
  )

  # return the log-likelihood values
  return(-ll_val)
}
