# This function computes the posterior distribution of abilities for examinees
# given their item response data, item parameters, and population distribution
posterior <- function(likehd, weights, idx.std = NULL) {
  # check if MG-calibration is implemented
  if (is.null(idx.std)) {
    # count the number of quad points
    ntheta <- nrow(weights)
  } else {
    # count the number of quad points
    ntheta <- nrow(weights[[1]])

    # divide the likelihood matrix into several groups if idx.std is not NULL
    # this is only for MG-calibration
    likehd.gr <- purrr::map(.x = idx.std, ~ {
      likehd[.x, ]
    })
  }

  # create the joint likelihood matrix of likelihoods and population distribution
  if (is.null(idx.std)) {
    joint_like <- sweep(likehd, 2, weights[, 2], FUN = "*")
  } else {
    joint_like <- array(0, c(nrow(likehd), ntheta))
    for (g in 1:length(weights)) {
      joint_like[idx.std[[g]], ] <- sweep(likehd.gr[[g]], 2, weights[[g]][, 2], FUN = "*")
    }
  }

  # denominator of the posterior distribution
  denom <- Rfast::rowsums(joint_like)

  # compute the posterior distribution of examinees across the node values
  # a row and column indicate the examinees and nodes, respectively
  posterior <- joint_like / denom

  # return results
  posterior
}
