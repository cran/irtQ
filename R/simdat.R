#' Simulated Response Data
#'
#' This function generates simulated response data for single-format or
#' mixed-format test forms. For dichotomous item response data, the IRT 1PL,
#' 2PL, and 3PL models are supported. For polytomous item response data, the
#' graded response model (GRM), the partial credit model (PCM), and the
#' generalized partial credit model (GPCM) are supported.
#'
#' @param x A data frame containing item metadata. This metadata is required to
#'   retrieve essential information for each item (e.g., number of score
#'   categories, IRT model type, etc.) necessary for calibration. You can create
#'   an empty item metadata frame using the function [irtQ::shape_df()]. See
#'   **below** for more details. Default is `NULL`.
#' @param theta A numeric vector of ability (theta) values.
#' @param a.drm A numeric vector of item discrimination (slope) parameters for
#'   dichotomous IRT models.
#' @param b.drm A numeric vector of item difficulty parameters for dichotomous
#'   IRT models.
#' @param g.drm A numeric vector of guessing parameters for dichotomous IRT
#'   models.
#' @param a.prm A numeric vector of item discrimination (slope) parameters for
#'   polytomous IRT models.
#' @param d.prm A list of numeric vectors, where each vector contains difficulty
#'   (threshold) parameters for a polytomous item.
#' @param cats A numeric vector indicating the number of score categories for
#'   each item.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.7 is commonly
#'   used for this purpose. Default is 1.
#' @param pr.model A character vector specifying the polytomous IRT model used
#'   to simulate responses for each polytomous item. Each element should be
#'   either "GRM" (graded response model) or "GPCM" (generalized partial credit
#'   model).
#'
#' @details There are two ways to generate simulated response data. The first is
#'   by providing a data frame of item metadata using the argument `x`. This
#'   data frame must follow a specific structure: the first column should
#'   contain item IDs, the second column should contain the number of unique
#'   score categories for each item, and the third column should specify the IRT
#'   model to be fitted to each item. Available IRT models are:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, and `"DRM"` for dichotomous item data
#'   - `"GRM"` and `"GPCM"` for polytomous item data
#'
#'   Note that `"DRM"` serves as a general label covering all dichotomous IRT
#'   models (i.e., `"1PLM"`, `"2PLM"`, and `"3PLM"`), while `"GRM"` and `"GPCM"`
#'   represent the graded response model and (generalized) partial credit model,
#'   respectively.
#'
#'   The subsequent columns should contain the item parameters for the specified
#'   models. For dichotomous items, the fourth, fifth, and sixth columns
#'   represent item discrimination (slope), item difficulty, and item guessing
#'   parameters, respectively. When `"1PLM"` or `"2PLM"` is specified in the
#'   third column, `NA`s must be entered in the sixth column for the guessing
#'   parameters.
#'
#'   For polytomous items, the item discrimination (slope) parameter should
#'   appear in the fourth column, and the item difficulty (or threshold)
#'   parameters for category boundaries should occupy the fifth through the last
#'   columns. When the number of unique score categories differs across items,
#'   unused parameter cells should be filled with `NA`s.
#'
#'   In the \pkg{irtQ} package, the threshold parameters for GPCM items are
#'   expressed as the item location (or overall difficulty) minus the threshold
#'   values for each score category. Note that when a GPCM item has *K* unique
#'   score categories, *K - 1* threshold parameters are required, since the
#'   threshold for the first category boundary is always fixed at 0. For
#'   example, if a GPCM item has five score categories, four threshold
#'   parameters must be provided.
#'
#'   An example of a data frame for a single-format test is shown below:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
#'   ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
#'   ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
#' }
#'
#'   An example of a data frame for a mixed-format test is shown below:
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
#'
#'   See the *IRT Models* section in the [irtQ-package] documentation for more
#'   details about the IRT models used in the \pkg{irtQ} package. A convenient
#'   way to create a data frame for the argument `x` is by using the function
#'   [irtQ::shape_df()].
#'
#'   The second approach is to simulate response data by directly specifying
#'   item parameters, instead of providing a metadata data frame via the `x`
#'   argument (see examples below). In this case, the following arguments must
#'   also be specified: `theta`, `cats`, `pr.model`, and `D`.
#'
#'   The `g.drm` argument is only required when simulating dichotomous item
#'   responses under the 3PL model. It can be omitted entirely if all
#'   dichotomous items follow the 1PL or 2PL model. However, if the test
#'   includes a mixture of 1PL, 2PL, and 3PL items, the `g.drm` vector must be
#'   specified for all items, using `NA` for non-3PL items. For example, if a
#'   test consists of four dichotomous items where the first two follow the 3PL
#'   model and the third and fourth follow the 1PL and 2PL models respectively,
#'   then `g.drm = c(0.2, 0.1, NA, NA)` should be used.
#'
#'   For dichotomous items, each element in `cats` should be set to 2. For
#'   polytomous items, the number of unique score categories should be specified
#'   in `cats`. When simulating data for a mixed-format test, it is important to
#'   specify `cats` in the correct item order. For example, suppose responses
#'   are simulated for 10 examinees across 5 items, including 3 dichotomous
#'   items and 2 polytomous items (each with 3 categories), where the second and
#'   fourth items are polytomous. In this case, `cats = c(2, 3, 2, 3, 2)` should
#'   be used.
#'
#'   Furthermore, if the two polytomous items are modeled using the graded
#'   response model and the generalized partial credit model, respectively, then
#'   `pr.model = c("GRM", "GPCM")`.
#'
#' @return A matrix or vector of simulated item responses.
#'   If a matrix is returned, rows correspond to examinees (theta values) and
#'   columns to items.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::drm()], [irtQ::prm()]
#'
#' @examples
#' ## Example 1:
#' ## Simulate response data for a mixed-format test.
#' ## The first two polytomous items use the generalized partial credit model (GPCM),
#' ## and the last polytomous item uses the graded response model (GRM).
#' # Generate theta values for 100 examinees
#' theta <- rnorm(100)
#'
#' # Set item parameters for three dichotomous items under the 3PL model
#' a.drm <- c(1, 1.2, 1.3)
#' b.drm <- c(-1, 0, 1)
#' g.drm <- rep(0.2, 3)
#'
#' # Set item parameters for three polytomous items
#' # These items have 4, 4, and 5 response categories, respectively
#' a.prm <- c(1.3, 1.2, 1.7)
#' d.prm <- list(c(-1.2, -0.3, 0.4), c(-0.2, 0.5, 1.6), c(-1.7, 0.2, 1.1, 2.0))
#'
#' # Specify the number of score categories for all items
#' # This vector also determines the location of polytomous items
#' cats <- c(2, 2, 4, 4, 5, 2)
#'
#' # Specify the IRT models for the polytomous items
#' pr.model <- c("GPCM", "GPCM", "GRM")
#'
#' # Simulate the response data
#' simdat(
#'   theta = theta, a.drm = a.drm, b.drm = b.drm, g.drm = NULL,
#'   a.prm = a.prm, d.prm = d.prm, cats = cats, pr.model = pr.model, D = 1
#' )
#'
#' ## Example 2:
#' ## Simulate response data for a single-format test using the 2PL model
#' # Specify score categories (2 for each dichotomous item)
#' cats <- rep(2, 3)
#'
#' # Simulate the response data
#' simdat(theta = theta, a.drm = a.drm, b.drm = b.drm, cats = cats, D = 1)
#'
#' ## Example 3:
#' ## Simulate response data using a "-prm.txt" file exported from flexMIRT
#' # Load the flexMIRT parameter file
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Convert the flexMIRT parameters to item metadata
#' test_flex <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Simulate the response data using the item metadata
#' simdat(x = test_flex, theta = theta, D = 1)
#'
#' @importFrom stats na.exclude
#' @export
#'
simdat <- function(x = NULL,
                   theta,
                   a.drm,
                   b.drm,
                   g.drm = NULL,
                   a.prm,
                   d.prm,
                   cats,
                   pr.model,
                   D = 1) {

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
