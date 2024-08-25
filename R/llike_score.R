#' Loglikelihood of Ability Parameters
#'
#' @description This function computes the loglikelihood of ability parameters
#' given the item parameters and response data.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number
#' of score categories, models). This can be created easily using the \code{\link{shape_df}}
#' function. See \code{\link{est_irt}}, \code{\link{irtfit}}, \code{\link{info}},
#' or \code{\link{simdat}} for more details about the item metadata.
#' @param data A matrix representing examinees' response data for the items in \code{x}.
#' Each row and column corresponds to an examinee and an item, respectively.
#' @param theta A numeric vector of ability parameters for which the loglikelihood
#' values will be computed.
#' @param D A scaling factor in IRT models that adjusts the logistic function to
#' approximate the normal ogive function (set to 1.7). The default is 1.
#' @param method A character string specifying the scoring method. Options include
#' "ML" for maximum likelihood estimation,  "MLF" for maximum likelihood estimation
#' with fences, and "MAP" for maximum a posteriori estimation. The default method is "MLE".
#' @param norm.prior A numeric vector of two elements indicating the mean and standard
#' deviation of the normal prior distribution. These parameters are used to obtain
#' the Gaussian quadrature points and corresponding weights from the normal distribution.
#' Default is c(0,1). This parameter is ignored if \code{method} is "ML" or "MLF".
#' @param fence.a A numeric value defining the item slope parameter (a-parameter) for
#' the two imaginary items in the MLF method. Default is 3.0.
#' @param fence.b A numeric vector of two elements specifying the lower and upper fences
#' of item difficulty parameters (b-parameters) for the two imaginary items in the MLF method.
#' If \code{fence.b = NULL}, the \code{range} values are used to set the fences.
#' The default is NULL.
#' @param missing A value used to denote missing values in the response data set. Default is NA.
#'
#' @details The function computes the loglikelihood value of the ability parameter given
#' the item parameters and response data for each item. As an example, to assess the loglikelihoods
#' of abilities for two examinees who have taken the same test items specified in \code{x}, supply
#' their item response data matrix with two rows in \code{data} and a vector of ability values
#' for which loglikelihood needs to be computed in \code{theta}.
#'
#' @return A data frame of loglikelihood values. Each row indicates the ability parameter
#' for which the loglikelihood was computed, and each column represents a response pattern.
#'
#' @examples
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameters and transform them into item metadata
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' # Generate examinees' abilities from N(0, 1)
#' set.seed(10)
#' score <- rnorm(5, mean=0, sd=1)
#'
#' # Simulate the response data
#' data <- simdat(x=x, theta=score, D=1)
#'
#' # Specify the ability values for which the loglikelihood values will be computed
#' theta <- seq(-3, 3, 0.5)
#'
#' # Compute the loglikelihood values (using the MLE method)
#' llike_score(x=x, data=data, theta=theta, D=1, method="ML")
#'
#' @export
#' @importFrom reshape2 melt
#' @import dplyr
llike_score <- function(x, data, theta, D = 1, method = "ML", norm.prior = c(0, 1),
                        fence.a = 3.0, fence.b = NULL, missing = NA) {
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
