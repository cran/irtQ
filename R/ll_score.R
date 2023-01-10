# This function computes a negative log likelihood or likelihood for scoring
#' @importFrom Rfast rowprods
#' @importFrom stats na.exclude
ll_score <- function(theta, elm_item, freq.cat, method=c("ML", "MAP", "MLF"),
                     idx.drm, idx.prm, D=1, norm.prior=c(0, 1), one.theta=TRUE,
                     logL=TRUE) {

  method <- match.arg(method)

  # when one theta value was used (i.e., length(theta) = 1)
  if(one.theta) {

    # compute the probabilities for DRM items
    if(!is.null(idx.drm)) {

      # compute the probabilities of correct answers
      p.drm <- drm(theta=theta, a=elm_item$pars[idx.drm, 1],
                   b=elm_item$pars[idx.drm, 2], g=elm_item$pars[idx.drm, 3],
                   D=D)

      # compute 1 - p
      q.drm <- 1 - p.drm

      # probabilities of the responses
      prob.rp.drm <- c(q.drm * freq.cat[idx.drm, 1] + p.drm * freq.cat[idx.drm, 2])

    } else {

      # probabilities of the responses
      prob.rp.drm <- NULL

    }

    # compute the probabilities for PRM items
    if(!is.null(idx.prm)) {

      # probabilities of responses across all items
      n.prm <- length(idx.prm)
      prob.rp.prm <- rep(NA, n.prm)
      for(k in 1:n.prm) {
        par.tmp <- stats::na.exclude(elm_item$par[idx.prm[k], ])
        ps.tmp <- prm(theta=theta, a=par.tmp[1], d=par.tmp[-1],
                      D=D, pr.model=elm_item$model[idx.prm[k]])
        cats.tmp <- elm_item$cats[idx.prm[k]]
        prob.rp.prm[k] <- ps.tmp[, freq.cat[idx.prm[k], c(1:cats.tmp)] == 1]
      }

    } else {

      # probabilities of the responses
      prob.rp.prm <- NULL

    }

    # compute the log-likelihood
    if(logL) {

      # sum of loglikelihood
      rst <- sum(log(c(prob.rp.drm, prob.rp.prm)))

      # compute the negative loglikelihood values
      rst <- switch(method,
                    ML = -rst,
                    MLF = -rst,
                    MAP = -(rst + stats::dnorm(theta, mean=norm.prior[1], sd=norm.prior[2], log=TRUE))
      )

    } else {

      # compute the likelihood
      rst <- prod(prob.rp.drm) * prod(prob.rp.prm)

    }

  } else {

    # when multiple theta value were used (i.e., length(theta) != 1)
    # compute the probabilities for DRM items
    if(!is.null(idx.drm)) {

      # compute the probabilities of correct answers
      p.drm <- drm(theta=theta, a=elm_item$pars[idx.drm, 1],
                   b=elm_item$pars[idx.drm, 2], g=elm_item$pars[idx.drm, 3],
                   D=D)

      # compute 1 - p
      q.drm <- 1 - p.drm

      # probabilities of the responses
      prob.rp.drm <-
        cbind(q.drm[, freq.cat[idx.drm, 1] == 1, drop = FALSE],
              p.drm[, freq.cat[idx.drm, 2] == 1, drop = FALSE])


    } else {

      # probabilities of the responses
      prob.rp.drm <- array(1, c(1, 1))

    }

    # compute the probabilities for PRM items
    if(!is.null(idx.prm)) {

      # probabilities of responses across all items
      n.prm <- length(idx.prm)
      prob.rp.prm <- vector('list', n.prm)
      for(k in 1:n.prm) {
        par.tmp <- stats::na.exclude(elm_item$par[idx.prm[k], ])
        ps.tmp <- prm(theta=theta, a=par.tmp[1], d=par.tmp[-1],
                      D=D, pr.model=elm_item$model[idx.prm[k]])
        cats.tmp <- elm_item$cats[idx.prm[k]]
        prob.rp.prm[[k]] <- ps.tmp[, freq.cat[idx.prm[k], c(1:cats.tmp)] == 1]
      }
      prob.rp.prm <- do.call(what='cbind', prob.rp.prm)

    } else {

      # probabilities of the responses
      prob.rp.prm <- array(1, c(1, 1))

    }


    # compute the log-likelihood
    if(logL) {

      # sum of loglikelihood
      rst <- Rfast::rowsums(log(prob.rp.drm)) + Rfast::rowsums(log(prob.rp.prm))

      # compute the negative loglikelihood values
      rst <- switch(method,
                    ML = -rst,
                    MLF = -rst,
                    MAP = -(rst + stats::dnorm(theta, mean=norm.prior[1], sd=norm.prior[2], log=TRUE))
      )

    } else {

      # compute the likelihood
      rst <- Rfast::rowprods(prob.rp.drm) * Rfast::rowprods(prob.rp.prm)

    }
  }

  # return the results
  rst

}
