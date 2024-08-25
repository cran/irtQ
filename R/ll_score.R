# This function computes a negative log likelihood or likelihood for scoring
#' @importFrom Rfast rowprods rowsums colsums
ll_score <- function(theta, elm_item, freq.cat, method = c("ML", "WL", "MAP", "MLF"),
                     idx.drm, idx.prm, D = 1, norm.prior = c(0, 1),
                     logL = TRUE) {
  method <- match.arg(method)

  # compute the probabilities for DRM items
  if (!is.null(idx.drm)) {
    # extract the response and item parameters
    r_i <- freq.cat[idx.drm, 2]
    a <- elm_item$pars[idx.drm, 1]
    b <- elm_item$pars[idx.drm, 2]
    g <- elm_item$pars[idx.drm, 3]

    # compute the probabilities of correct answers
    p.drm <-
      info_drm(
        theta = theta, a = a, b = b, g = g,
        D = D, one.theta = FALSE, r_i = r_i,
        grad = FALSE, info = FALSE
      )$P

    # compute 1 - p
    q.drm <- 1 - p.drm

    # probabilities of the responses
    prob.rp.drm <- q.drm * freq.cat[idx.drm, 1] + p.drm * freq.cat[idx.drm, 2]
  } else {
    # probabilities of the responses
    prob.rp.drm <- array(1, c(1, 1))
  }

  # compute the probabilities for PRM items
  if (!is.null(idx.prm)) {
    # count the number of examinees
    nstd <- length(theta)

    # check what poly models were used
    pr.mod <- unique(elm_item$model[idx.prm])
    prob.rp.prm <- NULL
    for (mod in pr.mod) {
      # extract the response, model, and item parameters
      lg.prm <- elm_item$model == mod
      par.tmp <- elm_item$par[lg.prm, , drop = FALSE]
      r_i <- freq.cat[lg.prm, , drop = FALSE]
      a <- par.tmp[, 1]
      d <- par.tmp[, -1, drop = FALSE]

      # compute the probabilities of endorsing each score category
      P.all <-
        info_prm(
          theta = theta, a = a, d = d, D = D, pr.model = mod,
          r_i = r_i, grad = FALSE, info = FALSE
        )$P
      P.sel <- do.call(
        rbind,
        replicate(n = nstd, expr = r_i, simplify = FALSE)
      ) * P.all
      Ps <-
        matrix(Rfast::rowsums(P.sel), ncol = nstd, byrow = FALSE)
      prob.rp.prm <- rbind(prob.rp.prm, Ps)
    }
  } else {
    # probabilities of the responses
    prob.rp.prm <- array(1, c(1, 1))
  }

  # compute the log-likelihood
  if (logL) {
    # sum of the loglikelihood
    rst <- Rfast::colsums(log(prob.rp.drm)) + Rfast::colsums(log(prob.rp.prm))

    # compute the negative loglikelihood values
    rst <- switch(method,
      ML = -rst,
      WL = -rst,
      MLF = -rst,
      MAP = -(rst + stats::dnorm(theta, mean = norm.prior[1], sd = norm.prior[2], log = TRUE))
    )
  } else {
    # compute the likelihood
    rst <- Rfast::colprods(prob.rp.drm) * Rfast::colprods(prob.rp.prm)
  }

  # return the results
  rst
}
