# This function computes a matrix of log-likelihoods for obtaining the observed item responses across examinees given theta values
# This function is used for MMLE-EM algorithm
#' @importFrom Matrix tcrossprod
likelihood <- function(elm_item, idx.drm = NULL, idx.prm = NULL,
                       data_drm = NULL, data_prm = NULL, theta, D = 1) {
  # compute the log-likelihoods given the quadrature points
  # when there are DRM items
  if (!is.null(idx.drm)) {
    # compute the probabilities of answering correctly on items
    ps <- drm(
      theta = theta, a = elm_item$pars[idx.drm, 1], b = elm_item$par[idx.drm, 2],
      g = elm_item$par[idx.drm, 3], D = D
    )

    # compute the probabilities of answering incorrectly on items
    qs <- 1 - ps

    # compute the log of the probabilities
    log_qps <- log(cbind(qs, ps))

    # compute the loglikelihood values for all examinees at each quadrature point
    # a row indicate the examinee and column indicate the quad point
    llike_drm <- Matrix::tcrossprod(x = data_drm, y = log_qps)
  } else {
    llike_drm <- 0L
  }

  # when there are PRM items
  if (!is.null(idx.prm)) {
    # compute the category probabilities of items
    n.prm <- length(idx.prm)
    prob.prm <- vector("list", n.prm)
    for (k in 1:n.prm) {
      par.tmp <- stats::na.exclude(elm_item$par[idx.prm[k], ])
      prob.prm[[k]] <- prm(
        theta = theta, a = par.tmp[1], d = par.tmp[-1],
        D = D, pr.model = elm_item$model[idx.prm[k]]
      )
    }
    prob.prm <- do.call(what = "cbind", prob.prm)

    # to prevent that log(prob.prm) have NaN or -Inf values
    log_prob.prm <- log(prob.prm)

    # compute the loglikelihood values for all examinees at each quadrature point
    llike_prm <- Matrix::tcrossprod(x = data_prm, y = log_prob.prm)
  } else {
    llike_prm <- 0L
  }

  # sum of log-likelihood matrix
  LL <- base::as.matrix(llike_drm + llike_prm)

  # transform to likelihood matrix
  L <- exp(LL)

  # return results
  rst <- list(L = L, LL = LL)
  rst
}
