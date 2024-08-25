#' Simulated Response Data
#'
#' @description This function generates a simulated response data for a single- or a mixed-format test forms. For dichotomous
#' item response data, the IRT 1PL, 2PL, and 3PL models are available. For polytomous item response data, the graded response model,
#' the partial credit model, and the generalized partial credit model are available.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...). This data frame
#' can be easily obtained using the function \code{\link{shape_df}}. See below for details.
#' @param theta A vector of theta values.
#' @param a.drm A vector of item discrimination (or slope) parameters for dichotomous response IRT models.
#' @param b.drm A vector of item difficulty (or threshold) parameters for dichotomous response IRT models.
#' @param g.drm A vector of item guessing parameters for dichotomous IRT models.
#' @param a.prm A vector of item discrimination (or slope) parameters for polytomous response IRT models.
#' @param d.prm A list containing vectors of item threshold (or step) parameters for polytomous response IRT models.
#' @param cats A vector containing the number of score categories for items.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param pr.model A vector of character strings specifying the polytomous model with which response data are simulated.
#' For each polytomous model, "GRM" for the graded response model or "GPCM" for the (generalized) partial credit model can be
#' specified.
#'
#' @details There are two ways of generating the simulated response data.
#' The first way is by using the argument \code{x} to read in a data frame of item metadata. In the data frame, the first column should have item IDs,
#' the second column should contain unique score category numbers of the items, and the third column should include IRT models being fit to the items.
#' The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous item data, and "GRM" and "GPCM" for polytomous item data.
#' Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
#' response model and (generalized) partial credit model, respectively. The next columns should include the item parameters of the fitted IRT models.
#' For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
#' item guessing parameters, respectively. When "1PLM" and "2PLM" are specified in the third column, NAs should be inserted in the sixth column
#' for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be included in the
#' fourth column and the item difficulty (or threshold) parameters of category boundaries should be contained from the fifth to the last columns.
#' When the number of unique score categories differs between items, the empty cells of item parameters should be filled with NAs.
#' In the \pkg{irtQ} package, the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as
#' the item location (or overall difficulty) parameter subtracted by the threshold parameter for unique score categories of the item.
#' Note that when an GPCM item has \emph{K} unique score categories, \emph{K-1} item difficulty parameters are necessary because
#' the item difficulty parameter for the first category boundary is always 0. For example, if an GPCM item has five score categories,
#' four item difficulty parameters should be specified. An example of a data frame with a single-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
#'   ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
#'   ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
#' }
#' And an example of a data frame for a mixed-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
#'   ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
#'   ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
#'   ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
#'   ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
#'   ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
#' }
#' See \code{IRT Models} section in the page of \code{\link{irtQ-package}} for more details about the IRT models used in the \pkg{irtQ} package.
#' An easier way to create a data frame for the argument \code{x} is by using the function \code{\link{shape_df}}.
#'
#' The second way is by directly specifying item parameters for each item for which response data should be simulated
#' (i.e., without using a data frame, as shown in the examples that follow). In addition to item parameters,
#' \code{theta}, \code{cats}, \code{pr.model}, and  \code{D} should be specified as well. \code{g.drm} does not need to be specified when only
#' the 1PL and 2PL models are used for dichotomous item response data. For dichotomous items, 2s should be specified in \code{cats}.
#' For polytomous items, the number of unique score categories should be specified in \code{cats}. When a response data set is generated with
#' a mixed-format test, it is important to clearly specify \code{cats} according to the order of items in the test form. Suppose that the response
#' data of ten examinees are simulated with five items, including three dichotomous items and two polytomous items with three categories.
#' Also, suppose that the second and the forth items are the polytomous items. Then, \code{cats = c(2, 3, 2, 3, 2)} should be used.
#' Additionally, among those two polytomous items, if the first and second item response data are simulated from the graded response model
#' and generalized partial credit model, respectively, then \code{pr.model = c('GRM', 'GPCM')}.
#'
#' @return This function returns a vector or a matrix. When a matrix is returned, rows indicate theta values and columns represent items.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{drm}}, \code{\link{prm}}
#'
#' @examples
#' ## example 1.
#' ## simulates response data with a mixed-format test.
#' ## for the first two polytomous items, the generalized partial credit model is used
#' ## for the last polytomous item, the graded response model is used
#' # 100 examinees are sampled
#' theta <- rnorm(100)
#'
#' # set item parameters for three dichotomous items with the 3PL model
#' a.drm <- c(1, 1.2, 1.3)
#' b.drm <- c(-1, 0, 1)
#' g.drm <- rep(0.2, 3)
#'
#' # set item parameters for three polytomous item parameters
#' # note that 4, 4, and 5 categories are used for polytomous items
#' a.prm <- c(1.3, 1.2, 1.7)
#' d.prm <- list(c(-1.2, -0.3, 0.4), c(-0.2, 0.5, 1.6), c(-1.7, 0.2, 1.1, 2.0))
#'
#' # create a numeric vector of score categories for both dichotomous and polytomous item data
#' # this score category vector is used to specify the location of the polytomous items
#' cats <- c(2, 2, 4, 4, 5, 2)
#'
#' # create a character vector of the IRT model for the polytomous items
#' pr.model <- c("GPCM", "GPCM", "GRM")
#'
#' # simulate the response data
#' simdat(
#'   theta = theta, a.drm = a.drm, b.drm = b.drm, g.drm = NULL,
#'   a.prm = a.prm, d.prm = d.prm, cats = cats, pr.model = pr.model, D = 1
#' )
#'
#'
#' ## example 2.
#' ## simulates response data with a single-format test with the 2PL model.
#' # create a numeric vector of score categories for the three 2PL model items
#' cats <- rep(2, 3)
#'
#' # simulate the response data
#' simdat(theta = theta, a.drm = a.drm, b.drm = b.drm, cats = cats, D = 1)
#'
#' ## example 3.
#' ## the use of a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # read item parameters and transform them to item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # simulate the response data
#' simdat(x = test_flex, theta = theta, D = 1) # use a data.frame of item meta information
#'
#' @importFrom stats na.exclude
#' @export
#'
simdat <- function(x = NULL, theta, a.drm, b.drm, g.drm = NULL, a.prm, d.prm, cats, pr.model, D = 1) {
  if (!is.null(x)) { # if the item metadata is inserted in the x argument

    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # break down the item metadata into several components
    comp <- breakdown(x)

    # set conditions
    nstd <- length(theta)
    nitem <- nrow(x)

    # create an empty matrix
    # res <- matrix(NA, nrow=nstd, ncol=nitem)
    res <- array(NA, c(nstd, nitem))

    # set initial numbers
    idx.all <- idxfinder(comp)
    idx.drm <- idx.all$idx.drm
    idx.prm <- idx.all$idx.prm

    # Simulate data
    # (1) for DRM items
    if (!is.null(idx.drm)) {
      res.drm <- simdat_drm(theta,
        a = comp$par[idx.drm, 1],
        b = comp$par[idx.drm, 2], g = comp$par[idx.drm, 3], D = D
      )
      res[, idx.drm] <- res.drm
    }

    # (2) for PRM items
    if (!is.null(idx.prm)) {
      for (i in idx.prm) {
        par.tmp <- stats::na.exclude(comp$par[i, ])
        res[, i] <- simdat_prm(theta,
          a = par.tmp[1], d = par.tmp[-1], D = D,
          pr.model = comp$model[i]
        )
      }
    }
  } else {
    # check whether argument is correctly specified
    if (missing(cats)) stop("Category of each item is missing", call. = FALSE)

    # Set conditions
    nstd <- length(theta)
    nitem <- length(cats)

    # Create an empty matrix
    # res <- matrix(NA, nrow=nstd, ncol=nitem)
    res <- array(NA, c(nstd, nitem))

    # Set initial numbers
    idx.drm <- which(cats == 2)
    idx.prm <- which(cats > 2)

    # Simulate data
    # (1) for DRM items
    if (any(cats == 2)) {
      res.drm <- simdat_drm(theta, a = a.drm, b = b.drm, g = g.drm, D = D)
      res[, idx.drm] <- res.drm
    }

    # (2) for PRM items
    if (any(cats > 2)) {
      for (i in 1:length(idx.prm)) {
        res[, idx.prm[i]] <- simdat_prm(theta, a = a.prm[i], d = d.prm[[i]], D = D, pr.model = pr.model[i])
      }
    }
  }

  res
}


# A function for generating binary data for one item
simdat_drm <- function(theta, a, b, g, D) {
  # Number of examinees
  nstd <- length(theta)

  # Number of items
  nitem <- length(a)

  # check the item guessing parameters
  if (is.null(g)) g <- rep(0, nitem)

  # calculate probability of correct answer
  z <- (D * a) * Rfast::Outer(x = theta, y = b, oper = "-")
  sim <- t(g + (1 - g) / (1 + exp(-z)))

  # Sample random variables from uniform dist
  tmp <- stats::runif(nstd * nitem, 0, 1)
  rv_unif <- array(tmp, c(nstd, nitem))

  # Simulated Response data for one item
  sim[sim >= rv_unif] <- 1
  sim[sim < rv_unif] <- 0

  # return the results
  sim
}


# A function for generating categorical data for one item
#' @importFrom Rfast colCumSums
#' @importFrom Rfast rowsums
simdat_prm <- function(theta, a, d, D, pr.model) {
  # Number of examinees
  nstd <- length(theta)

  # Calculate true probability for each category
  fit <- prm(theta, a, d, D, pr.model)

  # Sample random variables from uniform dist
  rv_unif <- stats::runif(nstd, 0, 1)

  # Simulated Response data for one item
  cumprob <- t(Rfast::colCumSums(t(fit)))
  matTF <- rv_unif >= cumprob
  res <- Rfast::rowsums(matTF)

  # return the results
  res
}
