#' Dichotomous Response Model (DRM) Probabilities
#'
#' @description This function computes the probability of correct answers for multiple items
#' for a given set of theta values using the IRT 1PL, 2PL, and 3PL models.
#'
#' @param theta A vector of ability values.
#' @param a A vector of item discrimination (or slope) parameters.
#' @param b A vector of item difficulty (or threshold) parameters.
#' @param g A vector of item guessing parameters.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible
#' to the normal ogive function (if set to 1.7). Default is 1.
#'
#' @details \code{g} does not need to be specified when the response probabilities of
#' the 1PL and 2PL models are computed.
#'
#' @return This function returns a matrix where a row indicates the ability and a column
#' represents the item.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{prm}}
#'
#' @examples
#' ## when vectors are used for both theta values and item parameters (3PLM)
#' drm(c(-0.1, 0.0, 1.5), a = c(1, 2), b = c(0, 1), g = c(0.2, 0.1), D = 1)
#'
#' ## when vectors are only used for item parameters (2PLM)
#' drm(0.0, a = c(1, 2), b = c(0, 1), D = 1)
#'
#' ## when vectors are only used for theta values (3PLM)
#' drm(c(-0.1, 0.0, 1.5), a = 1, b = 1, g = 0.2, D = 1)
#'
#' @importFrom Rfast Outer
#'
#' @export
drm <- function(theta, a, b, g = NULL, D = 1) {
  # count the numbers of items
  nitem <- length(a)

  # check the item guessing parameters
  if (is.null(g)) g <- rep(0, nitem)

  # calculate probability of correct answer
  z <- (D * a) * Rfast::Outer(x = theta, y = b, oper = "-")
  P <- t(g + (1 - g) / (1 + exp(-z)))

  # prevent that the probabilities are equal to 1L or g (or 0 in case of 1PLM and 2PLM)
  P[P > 9999999999e-10] <- 9999999999e-10
  lg.lessg <- P < (g + 1e-10)
  P[lg.lessg] <- P[lg.lessg] + 1e-10

  # return the probability matrix
  P
}


#' Polytomous Response Model (PRM) Probabilities (GRM and GPCM)
#'
#' @description This function computes the probability of selecting a specific category for an item
#' for a given set of theta values using the graded response model and (generalized) partial credit model.
#'
#' @param theta A vector of ability values.
#' @param a A numeric value of item discrimination (or slope) parameter.
#' @param d A vector of item difficulty (or threshold) parameters.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal
#' ogive function  (if set to 1.7). Default is 1.
#' @param pr.model A character string indicating the polytomous model being used. Available models are "GRM" for
#' the the graded response model and "GPCM" for the (generalized) partial credit model.
#'
#' @details When the category probabilities are computed for an item with the partial credit model, provide \code{a = 1} for that item.
#' When \code{model = "GPCM"}, \code{d} should include the item difficulty (or threshold) parameters. In the \pkg{irtQ} package,
#' the item difficulty (or threshold) parameters of category boundaries for GPCM are expressed as the item location (or overall difficulty)
#' parameter subtracted by the threshold parameter for unique score categories of the item. Note that when an GPCM item has \emph{K}
#' unique score categories, \emph{K-1} item difficulty parameters are necessary because the item difficulty parameter for the first category
#' boundary is always 0. For example, if an GPCM item has five score categories, four item difficulty parameters should be specified.
#' For more details about the parameterization of the (generalized) partial credit model, See \code{IRT Models} section
#' in the page of \code{\link{irtQ-package}}.
#'
#' @return This function returns a matrix where a row indicates the ability and a column represents
#' score categories of the item.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{drm}}, \code{\link{irtfit}}
#'
#' @examples
#' ## Category probabilities for an item with four categories
#' ## using a generalized partial credit model
#' prm(theta = c(-0.2, 0, 0.5), a = 1.4, d = c(-0.2, 0, 0.5), D = 1, pr.model = "GPCM")
#'
#' ## Category probabilities for an item with five categories
#' ## using a graded response model
#' prm(theta = c(-0.2, 0, 0.5), a = 1.2, d = c(-0.4, -0.2, 0.4, 1.5), D = 1, pr.model = "GRM")
#'
#' @export
prm <- function(theta, a, d, D = 1, pr.model = c("GRM", "GPCM")) {
  pr.model <- toupper(pr.model)
  if (pr.model == "GRM") {
    P <- grm(theta = theta, a = a, d = d, D = D)
  } else {
    P <- gpcm(theta = theta, a = a, d = d, D = D)
  }

  # return
  P
}

# IRT GPC model
#' @importFrom Rfast Outer colCumSums rowsums
gpcm <- function(theta, a, d, D = 1) {
  # add 0 of the the first category step (threshold) parameter
  # to the step parameter vector
  d <- c(0, d)

  # calculate category probabilities
  z <- (D * a) * (Rfast::Outer(x = theta, y = d, oper = "-"))
  cumsum_z <- t(Rfast::colCumSums(z))
  if (any(cumsum_z > 700)) {
    cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
  }
  if (any(cumsum_z < -700)) {
    cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
  }
  numer <- exp(cumsum_z) # numerator
  denom <- Rfast::rowsums(numer) # denominator
  P <- numer / denom

  # return
  P
}

# IRT GRM model
grm <- function(theta, a, d, D = 1) {
  # count the number of d parameters
  # m <- length(d)

  # calculate all the probabilities greater than equal to each threshold
  allP <- drm(theta = theta, a = a, b = d, g = 0, D = D)

  # calculate category probabilities
  P <- cbind(1, allP) - cbind(allP, 0)
  P[P > 9999999999e-10] <- 9999999999e-10
  P[P < 1e-10] <- 1e-10

  # return
  P
}
