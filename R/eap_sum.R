# Compute EAP Summed Score
#
# @description This function computes the expected a posterior (EAP) summed score (Thissen et al., 1995; Thissen & Orlando, 2001)
# for each examinee. The EAP summed score is the mean of the posterior density for the summed score (or observed score) given
# the item parameter estimates.
#
# @param x A data.frame containing the item metadata (e.g., item parameters, number of categories, models ...).
# See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item metadata.
# This data.frame can be easily obtained using the function \code{\link{shape_df}}.
# @param data A matrix containing examinees' response data for the test items. A row and column indicate examinees and items, respectively.
# @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
# These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution.
# Default is c(0,1).
# @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
# @param weights A two-column matrix or data.frame containing the theta values (in the first column) and the weights (in the second column)
# for the prior distribution. If missing, default values are used (see \code{norm.prior} and \code{nquad}).
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
# Default is 1.
#
# @details
#
# @return A list of the EAP summed scoring results
#
# @author Hwanggyu Lim \email{hglim83@@gmail.com}
# @export
# @examples
# ## the use of a "-prm.txt" file obtained from a flexMIRT
# flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
# x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#
# # simulate the item responses of examinees
# set.seed(15)
# theta <- rnorm(500)
# data <- simdat(x, theta, D=1)
#
# # estimate the abilities
# eap_sum(x, data, norm.prior=c(0, 1), nquad=41, D=1)
#
#' @importFrom Rfast colsums rowsums
eap_sum <- function(x, data, norm.prior = c(0, 1), nquad = 41, weights = NULL, D = 1) {
  ## ------------------------------------------------------------------------------------------------
  # check missing data
  # replace NAs with 0
  na.lg <- is.na(data)
  if (any(na.lg)) {
    data[na.lg] <- 0
    memo <- "Any missing responses are replaced with 0s. \n"
    warning(memo, call. = FALSE)
  }

  # generate quad nodes and weights
  if (is.null(weights)) {
    weights <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    weights <- data.frame(weights)
  }

  # estimate likelihoods using lord-wingersky algorithm
  lkhd <- lwrc(x = x, theta = weights[, 1], prob = NULL, D = D)

  # estimate EAP for sum scores
  ss.prob <- c(lkhd %*% weights[, 2])
  post <- t((t(lkhd) * weights[, 2])) / ss.prob
  tr_post <- t(post)
  eap.est <- Rfast::colsums(tr_post * weights[, 1])
  eap.est2 <- Rfast::colsums(tr_post * weights[, 1]^2)
  se.est <- sqrt(eap.est2 - eap.est^2)

  # assign the EAP summed scores to each examinee
  obs.score <- 0:(length(eap.est) - 1)
  names(eap.est) <- obs.score
  names(se.est) <- obs.score
  sumScore <- Rfast::rowsums(data)
  est_score <- eap.est[as.character(sumScore)]
  est_se <- se.est[as.character(sumScore)]
  score_table <-
    data.frame(
      sum.score = obs.score,
      est.theta = eap.est, se.theta = se.est,
      stringsAsFactors = FALSE
    )
  rownames(score_table) <- NULL

  # create a data frame for the estimated scores
  est.par <- data.frame(
    sum.score = sumScore,
    est.theta = est_score,
    se.theta = est_se
  )
  rownames(est.par) <- NULL

  # return results
  rst <- list(
    est.par = est.par,
    score.table = score_table
  )
  rst
}
