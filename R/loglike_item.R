# Negetive Loglikelihood of GPCM and GRM items
# @description This function computes the negative loglikelihood of an item with the
# polytomous IRT model
# @param item_par A vector of item parameters. The first element is the item discrimination (or slope)
# parameter. From the second elements, all all parameters are item threshold (or step) parameters.
# @param r_i A matrix of the frequencies for each score categories
# @param theta A vector of theta for an item.
# @param pr.mod A vector of character strings specifying the polytomous model with which response data are simulated.
# For each polytomous model, "GRM" for the graded response model or "GPCM" for the (generalized) partial credit model can be
# specified.
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal
# ogive function (if set to 1.7). Default is 1.
#
# @return A negative log-likelihood sum
loglike_prm <- function(item_par, r_i, theta, pr.mod = c("GRM", "GPCM"), D = 1, nstd,
                        fix.a = FALSE, a.val = 1,
                        aprior = list(dist = "lnorm", params = c(1, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        use.aprior = FALSE, use.bprior = FALSE) {
  ## -------------------------------------------------------------------------
  if (!fix.a) {
    # compute category probabilities for all thetas
    ps <- prm(theta, a = item_par[1], d = item_par[-1], D = D, pr.model = pr.mod)

    # compute log-likelihood
    log_ps <- log(ps)

    # log-likelihood
    llike <- sum(r_i * log_ps)

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + sum(ln.bprior)
    }
  } else {
    # compute category probabilities for all thetas
    ps <- prm(theta, a = a.val, d = item_par, D = D, pr.model = pr.mod)

    # compute loglikelihood
    log_ps <- log(ps)

    # log-likelihood
    llike <- sum(r_i * log_ps)

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par, is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + sum(ln.bprior)
    }
  }

  # return a negative log-likelihood
  -llike
}

# Negetive Loglikelihood of dichotomous item
# @description This function computes the negative loglikelihood of an item with the
# dichotomous IRT model
# @param item_par A vector of item parameters. The first element is the item discrimination (or slope)
# parameter, the second element is the item difficulty parameter, and the third element is the item guessing
# parameter. The third element is necessary only when the 3PLM is used.
# @param f_i A vector of the frequencies for correct answer + incorrect answer
# @param r_i A vector of the frequencies of correct answer
# @param s_i A vector of the frequencies of incorrect answer
# @param theta A vector of theta for an item.
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal
# ogive function (if set to 1.7). Default is 1.
# @param use.prior A logical value. If TRUE, a prior distribution specified in the argument \code{prior} is used when
# estimating item parameters of the IRT 3PLM. Default is TRUE.
#
# @return A negative log-likelihood sum
loglike_drm <- function(item_par, f_i, r_i, s_i, theta, mod = c("1PLM", "2PLM", "3PLM"), D = 1,
                        nstd, fix.a = FALSE, fix.g = FALSE, a.val = 1, g.val = .2, n.1PLM = NULL,
                        aprior = list(dist = "lnorm", params = c(1, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        gprior = list(dist = "beta", params = c(5, 17)),
                        use.aprior = FALSE, use.bprior = FALSE, use.gprior = TRUE) {
  # compute log-likelihood
  # (1) 1PLM: the slope parameters are constrained to be equal across the 1PLM items
  if (!fix.a & mod == "1PLM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # compute the negative log-likelihood values for all 1PLM items
    llike <- llike_drm(a = a, b = b, g = 0, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, D = D)

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + sum(ln.bprior)
    }
  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if (fix.a & mod == "1PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = a.val, b = item_par, g = 0,
      f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, D = D
    )

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par, is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }
  }

  # (3) 2PLM
  if (mod == "2PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = item_par[1], b = item_par[2], g = 0,
      f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, D = D
    )

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[2], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }
  }

  # (4) 3PLM
  if (!fix.g & mod == "3PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = item_par[1], b = item_par[2], g = item_par[3],
      f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, D = D
    )

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[2], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }

    # when the guessing parameter prior is used
    if (use.gprior) {
      ln.gprior <- logprior(
        val = item_par[3], is.aprior = FALSE, D = NULL, dist = gprior$dist,
        par.1 = gprior$params[1], par.2 = gprior$params[2]
      )
      llike <- llike + ln.gprior
    }
  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if (fix.g & mod == "3PLM") {
    # sum of log-likelihood
    llike <- llike_drm(
      a = item_par[1], b = item_par[2], g = g.val,
      f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, D = D
    )

    # when the slope parameter prior is used
    if (use.aprior) {
      ln.aprior <- logprior(
        val = item_par[1], is.aprior = TRUE, D = D, dist = aprior$dist,
        par.1 = aprior$params[1], par.2 = aprior$params[2]
      )
      llike <- llike + ln.aprior
    }

    # when the difficulty parameter prior is used
    if (use.bprior) {
      ln.bprior <- logprior(
        val = item_par[2], is.aprior = FALSE, D = NULL, dist = bprior$dist,
        par.1 = bprior$params[1], par.2 = bprior$params[2]
      )
      llike <- llike + ln.bprior
    }
  }

  # return a negative loglikelihood value
  -llike
}

# compute a sum of the log-likelihood value for each dichotomous item
llike_drm <- function(a, b, g, f_i, r_i, s_i, theta, D = 1) {
  # compute the probability of correct answers
  p <- drm(theta, a = a, b = b, g = g, D = D)

  # compute 1 - p
  q <- 1 - p

  # compute the log-likelihood
  log_p <- log(p)
  log_q <- log(q)

  # log-likelihood
  L <- r_i * log_p + s_i * log_q

  # sum of log-likelihood
  llike <- sum(L)

  # return
  llike
}
