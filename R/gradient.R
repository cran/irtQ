# This function analytically computes a gradient for scoring
grad_score <- function(theta, elm_item, freq.cat, idx.drm, idx.prm,
                       method = c("ML", "MAP", "MLF"), D = 1, norm.prior = c(0, 1)) {
  # For DRM items
  if (!is.null(idx.drm)) {
    # extract the response and item parameters
    r_i <- freq.cat[idx.drm, 2]
    a <- elm_item$pars[idx.drm, 1]
    b <- elm_item$pars[idx.drm, 2]
    g <- elm_item$pars[idx.drm, 3]

    # compute the gradient
    grad_drm <-
      grad_score_drm(theta = theta, a = a, b = b, g = g, r_i = r_i, D = D)
  } else {
    # assign 0 to the gradient
    grad_drm <- 0
  }

  # For PRM items
  if (!is.null(idx.prm)) {
    # count the number of the PRM items
    n.prm <- length(idx.prm)

    # compute the gradients for all items
    grad_prm <- c()
    for (i in 1:n.prm) {
      # extract the response, model, and item parameters
      par.tmp <- stats::na.exclude(elm_item$par[idx.prm[i], ])
      a <- par.tmp[1]
      d <- par.tmp[-1]
      pr.mod <- elm_item$model[idx.prm[i]]
      cats.tmp <- elm_item$cats[idx.prm[i]]
      r_i <- freq.cat[idx.prm[i], c(1:cats.tmp)]

      # compute the gradient
      grad_prm[i] <-
        grad_score_prm(theta = theta, a = a, d = d, r_i = r_i, pr.mod = pr.mod, D = D)
    }

    # sum of gradients for all PRM items
    grad_prm <- sum(grad_prm)
  } else {
    # assign 0 to the gradient
    grad_prm <- 0
  }

  # sum of the gradients for DRM and PRM items
  grad <- sum(grad_drm, grad_prm)

  # extract the gradient vector when MAP method is used
  if (method == "MAP") {
    # compute a gradient of prior distribution
    rst.prior <-
      logprior_deriv(
        val = theta, is.aprior = FALSE, D = NULL, dist = "norm",
        par.1 = norm.prior[1], par.2 = norm.prior[2]
      )

    # extract the gradient
    grad.prior <- attributes(rst.prior)$gradient

    # add the gradient
    grad <- sum(grad, grad.prior)
  }

  # return results
  grad
}


# This function computes a gradient of DRM items
grad_score_drm <- function(theta, a, b, g, r_i, D = 1) {
  # compute the probabilities of correct answers
  p <- drm(theta = theta, a = a, b = b, g = g, D = D)

  # compute the gradient
  grad <- -D * sum((a * (p - g) * (r_i - p)) / ((1 - g) * p))

  # return the gradient
  grad
}


# This function computes a gradient of PRM items
grad_score_prm <- function(theta, a, d, r_i, pr.mod, D = 1) {
  if (pr.mod == "GRM") {
    # count the number of d parameters
    m <- length(d)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta = theta, a = a, b = d, g = 0, D = D)
    allPst[allPst > 9999999999e-10] <- 9999999999e-10
    allPst[allPst < 1e-10] <- 1e-10
    allQst <- 1 - allPst[, , drop = FALSE]

    # calculate category probabilities
    P <- cbind(1, allPst) - cbind(allPst, 0)
    P[P > 9999999999e-10] <- 9999999999e-10
    P[P < 1e-10] <- 1e-10

    # compute the component values to get a gradient
    Da <- D * a
    deriv_Pstth <- (Da * allPst) * allQst
    deriv_Pth <- c(0, deriv_Pstth) - c(deriv_Pstth, 0)
    frac_rp <- r_i / P

    # compute the gradient
    grad <- -sum(frac_rp * deriv_Pth)
  }

  if (pr.mod == "GPCM") {
    # include zero for the step parameter of the first category
    d <- c(0, d)

    # check the number of step parameters
    m <- length(d) - 1

    # calculate category probabilities
    Da <- D * a
    z <- Da * (theta - d)
    cumsum_z <- t(Rfast::colCumSums(z))
    if (any(cumsum_z > 700)) {
      cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
    }
    if (any(cumsum_z < -700)) {
      cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
    }
    numer <- exp(cumsum_z) # numerator
    denom <- sum(numer) # denominator
    P <- numer / denom

    # compute the component values to get a gradient
    frac_rp <- r_i / P
    denom2 <- denom^2
    d1th_z <- Da * (1:(m + 1))
    d1th_denom <- sum(numer * d1th_z)
    deriv_Pth <- (numer / denom2) * (d1th_z * denom - d1th_denom)

    # compute the gradients
    grad <- -sum(frac_rp * deriv_Pth)
  }

  # return the results
  grad
}


# This function analytically computes a gradient vector of dichotomous item parameters
grad_item_drm <- function(item_par, f_i, r_i, s_i, theta, mod = c("1PLM", "2PLM", "3PLM"), D = 1,
                          nstd, fix.a = FALSE, fix.g = TRUE, a.val = 1, g.val = .2, n.1PLM = NULL,
                          aprior = list(dist = "lnorm", params = c(1, 0.5)),
                          bprior = list(dist = "norm", params = c(0.0, 1.0)),
                          gprior = list(dist = "beta", params = c(5, 17)),
                          use.aprior = FALSE,
                          use.bprior = FALSE,
                          use.gprior = TRUE) {
  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  # (1) 1PLM: the slope parameters are constrained to be equal across the 1PLM items
  if (!fix.a & mod == "1PLM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    r_p <- r_i - (f_i * p)
    theta_b <- -Rfast::Outer(x = b, y = theta, oper = "-")

    # compute the gradients of a and bs parameters
    gr_a <- -D * sum(theta_b * r_p)
    gr_b <- (D * a[1]) * Rfast::colsums(r_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a[1], is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[-1] <- attributes(rst.bprior)$gradient
    }
  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if (fix.a & mod == "1PLM") {
    # assign a and b parameters
    a <- a.val
    b <- item_par

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    r_p <- r_i - (f_i * p)

    # compute the gradients of a and bs parameters
    gr_b <- -as.numeric(-D * a * sum(r_p))

    # combine all gradients into a vector
    grad <- gr_b

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior <- attributes(rst.bprior)$gradient
    }
  }

  # (3) 2PLM
  if (mod == "2PLM") {
    # assign a and b parameters
    a <- item_par[1]
    b <- item_par[2]

    # compute the probabilities of correct and incorrect
    # p <- drm(theta=theta, a=a, b=b, g=0, D=D) + 1e-8
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    Da <- D * a
    r_p <- r_i - (f_i * p)
    theta_b <- theta - b

    # compute the gradients of a and bs parameters
    gr_a <- -D * sum(theta_b * r_p)
    gr_b <- Da * sum(r_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }
  }

  # (4) 3PLM
  if (!fix.g & mod == "3PLM") {
    # assign a, b, g parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- item_par[3]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = g, D = D)

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - (f_i * p)
    theta_b <- theta - b
    r_p_p <- r_p / p
    u1 <- p_g * r_p_p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and bs parameters
    gr_a <- -u2 * sum(theta_b * u1)
    gr_b <- u3 * sum(u1)
    gr_g <- -(1 / g_1) * sum(r_p_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b, gr_g)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }

    # extract the gradient vector when the guessing parameter prior is used
    if (use.gprior) {
      rst.gprior <-
        logprior_deriv(
          val = g, is.aprior = FALSE, D = NULL, dist = gprior$dist,
          par.1 = gprior$params[1], par.2 = gprior$params[2]
        )

      # extract the gradient vector
      grad.prior[3] <- attributes(rst.gprior)$gradient
    }
  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if (fix.g & mod == "3PLM") {
    # assign a, b, g parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- g.val

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = g, D = D)

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - (f_i * p)
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and bs parameters
    gr_a <- -u2 * sum(theta_b * u1)
    gr_b <- u3 * sum(u1)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = b, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }
  }

  # add the prior gradient vector
  grad <- grad + grad.prior

  # return results
  grad
}


# This function analytically computes a gradient vector of polytomous item parameters
grad_item_prm <- function(item_par, r_i, theta, pr.mod, D = 1, nstd, fix.a = FALSE, a.val = 1,
                          aprior = list(dist = "lnorm", params = c(1, 0.5)),
                          bprior = list(dist = "norm", params = c(0.0, 1.0)),
                          use.aprior = FALSE,
                          use.bprior = FALSE) {
  #   # count the number of item parameters to be estimated
  n.par <- length(item_par)

  ## -------------------------------------------------------------------------
  # compute the gradients
  # (1) GRM
  if (pr.mod == "GRM") {
    # assign a, b parameters
    a <- item_par[1]
    d <- item_par[-1]

    # check the number of step parameters
    m <- length(d)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta = theta, a = rep(a, m), b = d, g = 0, D = D)
    allQst <- 1 - allPst[, , drop = FALSE]

    # calculate category probabilities
    P <- (cbind(1, allPst) - cbind(allPst, 0))
    P[P > 9999999999e-10] <- 9999999999e-10
    P[P < 1e-10] <- 1e-10

    # compute the component values to get gradients
    Da <- D * a
    pq_st <- allPst * allQst
    q_p_st <- allQst - allPst
    w1 <- D * (-Rfast::Outer(x = d, y = theta, oper = "-"))
    w2 <- w1 * pq_st
    w3 <- cbind(0, w2) - cbind(w2, 0)
    frac_rp <- r_i / P
    w4 <- -Rfast::coldiffs(frac_rp)

    # compute the gradients of a and bs parameters
    gr_a <- -sum(frac_rp * w3)
    gr_b <- -Da * Rfast::colsums(pq_st * w4)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if (use.aprior) {
      rst.aprior <-
        logprior_deriv(
          val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
          par.1 = aprior$params[1], par.2 = aprior$params[2]
        )

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if (use.bprior) {
      rst.bprior <-
        logprior_deriv(
          val = d, is.aprior = FALSE, D = NULL, dist = bprior$dist,
          par.1 = bprior$params[1], par.2 = bprior$params[2]
        )

      # extract the gradient vector
      grad.prior[-1] <- attributes(rst.bprior)$gradient
    }
  }

  # (2) GPCM and PCM
  if (pr.mod == "GPCM") {
    if (!fix.a) {
      # For GPCM
      # assign a parameter
      a <- item_par[1]

      # include zero for the step parameter of the first category
      d <- c(0, item_par[-1])

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities
      Da <- D * a
      theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
      z <- Da * theta_d
      cumsum_z <- t(Rfast::colCumSums(z))
      if (any(cumsum_z > 700)) {
        cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
      }
      if (any(cumsum_z < -700)) {
        cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
      }
      numer <- exp(cumsum_z) # numerator
      denom <- Rfast::rowsums(numer)
      P <- (numer / denom)

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      frac_rp <- r_i / P
      w1 <- D * theta_d
      w1cum <- t(Rfast::colCumSums(w1))
      denom2 <- denom^2
      d1a_denom <- Rfast::rowsums(numer * w1cum)
      deriv_Pa <- (numer * (w1cum * denom - d1a_denom)) / denom2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_a <- -sum(frac_rp * deriv_Pa)
      gr_b <- c()
      for (k in 1:m) {
        gr_b[k] <- -sum(frac_rp * (-numer * (DaDenom_vec %*% dsmat[k + 1, , drop = FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # combine all gradients into a vector
      grad <- c(gr_a, gr_b)

      # create a prior gradient vector
      grad.prior <- rep(0, n.par)

      # extract the gradient vector when the slope parameter prior is used
      if (use.aprior) {
        rst.aprior <-
          logprior_deriv(
            val = a, is.aprior = TRUE, D = D, dist = aprior$dist,
            par.1 = aprior$params[1], par.2 = aprior$params[2]
          )

        # extract the gradient vector
        grad.prior[1] <- attributes(rst.aprior)$gradient
      }

      # extract the gradient vector when the difficulty parameter prior is used
      if (use.bprior) {
        rst.bprior <-
          logprior_deriv(
            val = d[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
            par.1 = bprior$params[1], par.2 = bprior$params[2]
          )

        # extract the gradient vector
        grad.prior[-1] <- attributes(rst.bprior)$gradient
      }
    } else {
      # for PCM
      # assign a parameter
      a <- a.val

      # include zero for the step parameter of the first category
      d <- c(0, item_par)

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities
      Da <- D * a
      theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
      z <- Da * theta_d
      cumsum_z <- t(Rfast::colCumSums(z))
      if (any(cumsum_z > 700)) {
        cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
      }
      if (any(cumsum_z < -700)) {
        cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
      }
      numer <- exp(cumsum_z) # numerator
      denom <- Rfast::rowsums(numer)
      P <- (numer / denom)

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      frac_rp <- r_i / P
      denom2 <- denom^2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_b <- c()
      for (k in 1:m) {
        gr_b[k] <- -sum(frac_rp * (-numer * (DaDenom_vec %*% dsmat[k + 1, , drop = FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # combine all gradients into a vector
      grad <- gr_b

      # create a prior gradient vector
      grad.prior <- rep(0, n.par)

      # extract the gradient vector when the difficulty parameter prior is used
      if (use.bprior) {
        rst.bprior <-
          logprior_deriv(
            val = d[-1], is.aprior = FALSE, D = NULL, dist = bprior$dist,
            par.1 = bprior$params[1], par.2 = bprior$params[2]
          )

        # extract the gradient vector
        grad.prior <- attributes(rst.bprior)$gradient
      }
    }
  }

  # add the prior gradient vector
  grad <- grad + grad.prior

  # return results
  grad
}



# This function is used to compute the standard errors of item parameter estimates using the cross-product method.
grad_llike <- function(item_par, f_i, r_i, theta, n.theta, mod = c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"),
                       D = 1, fix.a.1pl = TRUE, fix.a.gpcm = FALSE, fix.g = FALSE, a.val.1pl = 1, a.val.gpcm = 1,
                       g.val = .2, n.1PLM = NULL) {
  # for dichotomous models
  if (mod %in% c("1PLM", "2PLM", "3PLM")) {
    if (!is.null(n.1PLM)) {
      # 1PLM: when the item slope parameters are not constrained to be equal across all items
      # compute the gradient vectors
      grad <- grad_item_drm_se(
        item_par = item_par, f_i = f_i, r_i = r_i, theta = theta,
        mod = mod, D = D, fix.a = fix.a.1pl, n.1PLM = n.1PLM
      )
    } else {
      # for all other dichotomous models
      # compute the gradient vectors
      grad <- grad_item_drm_se(
        item_par = item_par, f_i = f_i, r_i = r_i, theta = theta,
        mod = mod, D = D, fix.a = fix.a.1pl, fix.g = fix.g,
        a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL
      )
    }
  } else {
    # for polytomous models
    # compute the gradient vectors
    grad <- grad_item_prm_se(
      item_par = item_par, r_i = r_i, theta = theta, n.theta = n.theta,
      pr.mod = mod, D = D, fix.a = fix.a.gpcm, a.val = a.val.gpcm
    )
  }

  # return results
  grad
}


# This function analytically computes a matrix of gradients of dichotomous item parameters across all examinees
# This function is used to compute the cross-product information matrix
grad_item_drm_se <- function(item_par, f_i, r_i, theta, mod = c("1PLM", "2PLM", "3PLM", "DRM"), D = 1,
                             fix.a = FALSE, fix.g = TRUE, a.val = 1, g.val = .2, n.1PLM = NULL) {
  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  # (1) 1PLM: the slope parameters are constrained to be equal across the 1PLM items
  if (!is.null(n.1PLM)) {
    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    r_p <- r_i - f_i * p
    theta_b <- -Rfast::Outer(x = b, y = theta, oper = "-")

    # compute the gradients of a and bs parameters
    gr_a <- -D * Rfast::rowsums(theta_b * r_p)
    gr_b <- (D * a[1]) * r_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if (is.null(n.1PLM)) {
    # assign a and b parameters
    a <- a.val
    b <- item_par

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    r_p <- r_i - f_i * p

    # compute the gradients of b parameter
    gr_b <- (D * a) * r_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_b)
  }

  # (3) 2PLM
  if (mod == "2PLM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    b <- item_par[2]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = 0, D = D)

    # compute the component values
    Da <- D * a
    r_p <- r_i - f_i * p
    theta_b <- theta - b

    # compute the gradients of a and b parameters
    gr_a <- (-D * theta_b) * r_p
    gr_b <- Da * r_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # (4) 3PLM
  if (!fix.g & mod == "3PLM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    b <- item_par[2]
    g <- item_par[3]

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = g, D = D)

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - f_i * p
    theta_b <- theta - b
    r_p_p <- r_p / p
    u1 <- p_g * r_p_p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a, b, and g parameters
    gr_a <- (-u2 * theta_b) * u1
    gr_b <- u3 * u1
    gr_g <- -(1 / g_1) * r_p_p

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b, gr_g)
  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if (fix.g & mod == "3PLM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    b <- item_par[2]
    g <- g.val

    # compute the probabilities of correct and incorrect
    p <- drm(theta = theta, a = a, b = b, g = g, D = D)

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i - f_i * p
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and b parameters
    gr_a <- (-u2 * theta_b) * u1
    gr_b <- u3 * u1

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # return results
  grad
}

# This function analytically computes a matrix of gradients of polytomous item parameters across all examinees
# This function is used to compute the cross-product information matrix
grad_item_prm_se <- function(item_par, r_i, theta, n.theta, pr.mod, D = 1, fix.a = FALSE, a.val = 1) {
  # transform item parameters as numeric values
  # item_par <- as.numeric(item_par)

  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  ## -------------------------------------------------------------------------
  # compute the gradients
  # (1) GRM
  if (pr.mod == "GRM") {
    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    d <- item_par[-1]

    # check the number of step parameters
    m <- length(d)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta = theta, a = rep(a, m), b = d, g = 0, D = D)
    allQst <- 1 - allPst[, , drop = FALSE]

    # calculate category probabilities
    P <- (cbind(1, allPst) - cbind(allPst, 0))
    P[P > 9999999999e-10] <- 9999999999e-10
    P[P < 1e-10] <- 1e-10

    # compute the component values to get gradients
    Da <- D * a
    pq_st <- allPst * allQst
    q_p_st <- allQst - allPst
    w1 <- D * (-Rfast::Outer(x = d, y = theta, oper = "-"))
    w2 <- w1 * pq_st
    w3 <- cbind(0, w2) - cbind(w2, 0)
    frac_rp <- r_i / P
    w4 <- -Rfast::coldiffs(frac_rp)

    # compute the gradients of a and bs parameters
    gr_a <- -Rfast::rowsums(frac_rp * w3)
    gr_b <- (-Da * pq_st) * w4

    # create a matrix of gradients across all examinees
    grad <- cbind(gr_a, gr_b)
  }

  # (2) GPCM and PCM
  if (pr.mod == "GPCM") {
    if (!fix.a) {
      # For GPCM
      # make vectors of a and b parameters for all 1PLM items
      a <- item_par[1]

      # include zero for the step parameter of the first category
      d <- c(0, item_par[-1])

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities
      Da <- D * a
      theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
      z <- Da * theta_d
      cumsum_z <- t(Rfast::colCumSums(z))
      if (any(cumsum_z > 700)) {
        cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
      }
      if (any(cumsum_z < -700)) {
        cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
      }
      numer <- exp(cumsum_z) # numerator
      denom <- Rfast::rowsums(numer)
      P <- (numer / denom)

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      # tr_dsmat <- t(dsmat)
      frac_rp <- r_i / P
      w1 <- D * theta_d
      w1cum <- t(Rfast::colCumSums(w1))
      denom2 <- denom^2
      d1a_denom <- Rfast::rowsums(numer * w1cum)
      deriv_Pa <- (numer * (w1cum * denom - d1a_denom)) / denom2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_a <- -Rfast::rowsums(frac_rp * deriv_Pa)
      gr_b <- array(0, c(n.theta, m))
      for (k in 1:m) {
        gr_b[, k] <-
          -Rfast::rowsums(frac_rp *
            (-numer * (DaDenom_vec %*%
              dsmat[k + 1, , drop = FALSE] +
              d1b_denom[, k + 1]) / denom2))
      }

      # create a matrix of gradients across all examinees
      grad <- cbind(gr_a, gr_b)
    } else {
      # for PCM
      # make vectors of a and b parameters for all 1PLM items
      a <- a.val

      # include zero for the step parameter of the first category
      d <- c(0, item_par)

      # check the number of step parameters
      m <- length(d) - 1

      # calculate category probabilities
      Da <- D * a
      theta_d <- (Rfast::Outer(x = theta, y = d, oper = "-"))
      z <- Da * theta_d
      cumsum_z <- t(Rfast::colCumSums(z))
      if (any(cumsum_z > 700)) {
        cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
      }
      if (any(cumsum_z < -700)) {
        cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
      }
      numer <- exp(cumsum_z) # numerator
      denom <- Rfast::rowsums(numer)
      P <- (numer / denom)

      # compute the component values to get a gradient vector
      dsmat <- array(0, c((m + 1), (m + 1)))
      dsmat[-1, -1][upper.tri(x = dsmat[-1, -1], diag = TRUE)] <- 1
      frac_rp <- r_i / P
      denom2 <- denom^2
      d1b_denom <- tcrossprod(x = (-Da * numer), y = dsmat)
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_b <- array(0, c(n.theta, m))
      for (k in 1:m) {
        gr_b[, k] <-
          -Rfast::rowsums(frac_rp *
            (-numer * (DaDenom_vec %*%
              dsmat[k + 1, , drop = FALSE] +
              d1b_denom[, k + 1]) / denom2))
      }

      # create a matrix of gradients across all examinees
      grad <- gr_b
    }
  }

  # return results
  grad
}
