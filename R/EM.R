# E-step function when FIPC method is used
#' @importFrom Matrix crossprod
Estep_fipc <- function(elm_item1, elm_item2, idx.drm2, idx.prm2,
                       data_drm2, data_prm2, data_all1, weights,
                       D = 1, idx.std = NULL) {
  # compute the likelihood and log-likelihood matrix of the fixed items
  # (in the first iteration of EM) or all items (in the rest of the iteration of EM)
  # this (log) likelihood matrix is used only for computing the posterior density of ability
  if (is.null(idx.std)) {
    L_LL <- likelihood(elm_item2,
      idx.drm = idx.drm2, idx.prm = idx.prm2,
      data_drm = data_drm2, data_prm = data_prm2, theta = weights[, 1], D = D
    )
  } else {
    L_LL <- likelihood(elm_item2,
      idx.drm = idx.drm2, idx.prm = idx.prm2,
      data_drm = data_drm2, data_prm = data_prm2, theta = weights[[1]][, 1], D = D
    )
  }

  # posterior distribution
  post_dist <- posterior(likehd = L_LL$L, weights = weights, idx.std = idx.std)

  # compute the expected frequency of scores categories across all items
  # this is the conditional expectation of item responses with respect to posterior likelihood distribution
  freq.exp <- base::as.matrix(Matrix::crossprod(post_dist, data_all1))

  # return results
  rst <- list(
    elm_item = elm_item1, post_dist = post_dist, freq.exp = freq.exp,
    loglikehd = L_LL$LL, likehd = L_LL$L, idx.std = idx.std
  )
  rst
}

# E-step function
#' @importFrom Matrix crossprod
Estep <- function(elm_item, idx.drm = NULL, idx.prm = NULL,
                  data_drm, data_prm, data_all, weights, D = 1, idx.std = NULL) {
  # compute the likelihood and log-likelihood matrix
  if (is.null(idx.std)) {
    L_LL <- likelihood(elm_item,
      idx.drm = idx.drm, idx.prm = idx.prm,
      data_drm = data_drm, data_prm = data_prm, theta = weights[, 1], D = D
    )
  } else {
    L_LL <- likelihood(elm_item,
      idx.drm = idx.drm, idx.prm = idx.prm,
      data_drm = data_drm, data_prm = data_prm, theta = weights[[1]][, 1], D = D
    )
  }

  # posterior distribution
  post_dist <- posterior(likehd = L_LL$L, weights = weights, idx.std = idx.std)

  # compute the expected frequency of scores categories across all items
  # this is the conditional expectation of item responses with respect to posterior likelihood distribution
  freq.exp <- base::as.matrix(Matrix::crossprod(post_dist, data_all))

  # return results
  rst <- list(
    elm_item = elm_item, post_dist = post_dist, freq.exp = freq.exp,
    loglikehd = L_LL$LL, likehd = L_LL$L, idx.std = idx.std
  )
  rst
}


# implement M-step
#' @importFrom Rfast rowsums colsums
Mstep <- function(estep, id, cats, model, quadpt, n.quad, D = 1, cols.item = NULL, loc_1p_const, loc_else, idx4est = NULL,
                  n.1PLM = NULL, EmpHist, weights, fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl, a.val.gpcm, g.val,
                  use.aprior, use.bprior, use.gprior, aprior, bprior, gprior, group.mean, group.var, nstd,
                  Quadrature, control, iter = NULL, fipc = FALSE, reloc.par,
                  ref.group = NULL, free.group = NULL, parbd = NULL) {
  # extract the results of E-step
  elm_item <- estep$elm_item
  post_dist <- estep$post_dist
  freq.exp <- estep$freq.exp
  loglikehd <- estep$loglikehd
  likehd <- estep$likehd
  idx.std <- estep$idx.std

  ## ----------------------------------------------------------------------
  # (1) item parameter estimation
  ## ----------------------------------------------------------------------
  # create empty vectors to contain results
  est_par <- NULL
  est_pure <- NULL
  convergence <- NULL
  noconv_items <- NULL
  se <- NULL

  if (!is.null(elm_item)) {
    # the dichotomous items: 1PLM with constrained slope values
    if (!is.null(loc_1p_const)) {
      # prepare input files to estimate the 1PLM item parameters
      f_i <- array(0, c(n.quad, n.1PLM))
      for (k in 1:n.1PLM) {
        cols.tmp <- cols.item$cols.all[[loc_1p_const[k]]]
        f_i[, k] <- Rfast::rowsums(freq.exp[, cols.tmp])
      }
      s_i <- freq.exp[, cols.item$cols.1pl][, c(TRUE, FALSE), drop = FALSE]
      r_i <- freq.exp[, cols.item$cols.1pl][, c(FALSE, TRUE), drop = FALSE]

      # set the starting values
      startval <-
        set_startval(
          pars = elm_item$pars, item = loc_1p_const,
          use.startval = TRUE, mod = "1PLM",
          score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
          fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
        )

      # bounds of the item parameters
      lower <- parbd$drm.slc$lower
      upper <- parbd$drm.slc$upper

      # item parameter estimation
      est <- estimation2(
        f_i = f_i, r_i = r_i, s_i = s_i, quadpt = quadpt, mod = "1PLM", D = D, n.quad = n.quad,
        fix.a.1pl = FALSE, n.1PLM = n.1PLM, aprior = aprior, bprior = bprior,
        use.aprior = use.aprior, use.bprior = use.bprior,
        control = control, startval = startval, lower = lower, upper = upper,
        iter = iter
      )

      # extract the results
      # item parameter estimates
      a <- est$pars[1]
      b <- est$pars[-1]
      pars <- purrr::map(1:n.1PLM, .f = function(x) c(a, b[x], 0))
      est_par <- c(est_par, pars)

      # convergence indicator
      convergence <- c(convergence, est$convergence)
      if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_1p_const)
    }

    # all other items
    if (length(loc_else) >= 1) {
      for (i in 1:length(loc_else)) {
        # prepare information to estimate item parameters
        mod <- model[loc_else][i]
        score.cat <- cats[loc_else][i]

        # in case of a DRM item
        if (score.cat == 2) {
          cols.tmp <- cols.item$cols.all[[loc_else[i]]]
          f_i <- Rfast::rowsums(freq.exp[, cols.tmp])
          s_i <- freq.exp[, cols.tmp[1]]
          r_i <- freq.exp[, cols.tmp[2]]

          # set the starting values
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = TRUE, mod = mod,
              score.cat = score.cat, fix.a.1pl = TRUE, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )

          # bounds of the item parameters
          loc.tmp <- which(idx4est$drm.else == loc_else[i])
          lower <- parbd$drm.else[[loc.tmp]]$lower
          upper <- parbd$drm.else[[loc.tmp]]$upper

          # item parameter estimation
          est <- estimation2(
            f_i = f_i, r_i = r_i, s_i = s_i, quadpt = quadpt, mod = mod, D = D,
            n.quad = n.quad, fix.a.1pl = ifelse(mod == "1PLM", TRUE, FALSE),
            fix.g = fix.g, a.val.1pl = a.val.1pl, g.val = g.val, n.1PLM = NULL,
            aprior = aprior, bprior = bprior, gprior = gprior,
            use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
            control = control, startval = startval, lower = lower, upper = upper,
            iter = iter
          )

          # extract the results
          # item parameter estimates
          a <- ifelse(mod == "1PLM", a.val.1pl, est$pars[1])
          b <- ifelse(mod == "1PLM", est$pars[1], est$pars[2])
          g <- ifelse(mod == "3PLM", ifelse(fix.g, g.val, est$pars[3]), 0)
          pars <- c(a, b, g)
          est_par <- c(est_par, list(pars))

          # convergence indicator
          convergence <- c(convergence, est$convergence)
          if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])
        }

        # in case of a PRM item
        if (score.cat > 2) {
          cols.tmp <- cols.item$cols.all[[loc_else[i]]]
          r_i <- freq.exp[, cols.tmp]

          # set the starting values
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = TRUE, mod = mod,
              score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )

          # bounds of the item parameters
          loc.tmp <- which(idx4est$prm == loc_else[i])
          lower <- parbd$prm[[loc.tmp]]$lower
          upper <- parbd$prm[[loc.tmp]]$upper

          # item parameter estimation
          est <- estimation2(
            r_i = r_i, quadpt = quadpt, mod = mod, score.cat = score.cat, D = D,
            n.quad = n.quad, fix.a.gpcm = ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
            a.val.gpcm = a.val.gpcm, n.1PLM = NULL, aprior = aprior, bprior = bprior,
            use.aprior = use.aprior, use.bprior = use.bprior,
            control = control, startval = startval, lower = lower, upper = upper,
            iter = iter
          )

          # extract the results
          # item parameter estimates
          a <- ifelse(mod == "GRM", est$pars[1], ifelse(fix.a.gpcm, a.val.gpcm, est$pars[1]))
          if (mod == "GRM") {
            bs <- est$pars[-1]
          } else {
            if (fix.a.gpcm) {
              bs <- est$pars
            } else {
              bs <- est$pars[-1]
            }
          }
          pars <- c(a, bs)
          est_par <- c(est_par, list(pars))

          # convergence indicator
          convergence <- c(convergence, est$convergence)
          if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])
        }
      }
    }
  }


  ## ---------------------------------------------------------------
  if (!is.null(elm_item)) {
    # arrange the estimated item parameters and standard errors
    par_df <- bind.fill(est_par, type = "rbind")
    par_df <- cbind(loc = c(loc_1p_const, loc_else), par_df)
    par_df <- par_df[order(par_df[, 1]), , drop = FALSE]
    par_df <- par_df[, -1, drop = FALSE]

    # create a full data.frame for the item parameter estimates
    colnames(par_df) <- paste0("par.", 1:ncol(par_df))
  } else {
    par_df <- NULL
    convergence <- NULL
    noconv_items <- NULL
  }

  ## ----------------------------------------------------------------------
  # (2) update the prior ability distribution
  ## ----------------------------------------------------------------------
  # divide the posterior dist matrix into each group when idx.std is not NULL
  # this is only for MG-calibration
  if (!is.null(idx.std)) {
    post_dist <- purrr::map(.x = idx.std, ~ {
      post_dist[.x, ]
    })
  }

  if (EmpHist) {
    # update the prior frequencies
    if (is.null(idx.std)) {
      # column sum across all quad points
      prior_freq <- prior_freq2 <- unname(Rfast::colsums(post_dist))

      # prevent that the frequency has less than 1e-20
      prior_freq[prior_freq < 1e-20] <- 1e-20

      # normalize the updated prior frequency to obtain prior density function
      prior_dense <- prior_freq2 / nstd

      # update the prior densities
      if (fipc) {
        # when FIPC is used, no rescaling is applied
        weights <- data.frame(theta = quadpt, weight = prior_dense)
      } else {
        # rescale the prior density distribution using the same quadrature point by applying Woods (2007) method
        prior_dense2 <-
          scale_prior(
            prior_freq = prior_freq, prior_dense = prior_dense, quadpt = quadpt,
            scale.par = c(group.mean, group.var), Quadrature = Quadrature
          )
        weights <- data.frame(theta = quadpt, weight = prior_dense2)
      }
    } else {
      # column sum across all quad points
      prior_freq <- prior_freq2 <-
        purrr::map(.x = post_dist, ~ {
          unname(Rfast::colsums(.x))
        })

      # prevent that the frequency has less than 1e-20
      prior_freq <-
        purrr::map(
          .x = prior_freq,
          .f = function(x) {
            x[x < 1e-20] <- 1e-20
            x
          }
        )

      # normalize the updated prior frequency to obtain prior density function
      prior_dense <- purrr::map2(.x = prior_freq2, .y = nstd, ~ {
        .x / .y
      })

      # divide the prior frequencies and densities into the reference and free groups
      prior_freq_ref <- prior_freq[ref.group]
      prior_dense_ref <- prior_dense[ref.group]
      prior_freq_free <- prior_freq[free.group]
      prior_dense_free <- prior_dense[free.group]

      # rescale the prior density distribution using the same quadrature point by applying Woods (2007) method
      # rescale the distribution of the reference group first
      prior_dense_ref2 <-
        purrr::map2(
          .x = prior_freq_ref, .y = prior_dense_ref,
          ~ {
            scale_prior(
              prior_freq = .x, prior_dense = .y,
              quadpt = quadpt, scale.par = c(group.mean, group.var), Quadrature = Quadrature
            )
          }
        )

      # extract the rescaled densities and quadrature points
      weights.ref <- purrr::map(.x = prior_dense_ref2, ~ {
        data.frame(theta = quadpt, weight = .x)
      })

      # replace the old densities with the rescaled densities for the reference group
      weights[ref.group] <- weights.ref

      # extract the densities of free groups
      weights.free <- purrr::map(.x = prior_dense_free, ~ {
        data.frame(theta = quadpt, weight = .x)
      })

      # a list of density functions for all reference and free groups
      weights[free.group] <- weights.free
    }
  } else {
    if (is.null(idx.std)) {
      if (fipc) {
        # update the prior frequencies
        prior_freq <- unname(Rfast::colsums(post_dist))

        # normalize the updated prior frequency to obtain prior density function
        prior_dense <- prior_freq / nstd

        # compute the mean and sd of the updated prior distribution
        moments <- cal_moment(node = quadpt, weight = prior_dense)
        mu <- moments[1]
        sigma <- sqrt(moments[2])

        # obtain the updated prior densities from the normal distribution
        weights <- gen.weight(dist = "norm", mu = mu, sigma = sigma, theta = quadpt)
      } else {
        weights <- weights
      }
    } else {
      # column sum across all quad points
      prior_freq <- prior_freq2 <-
        purrr::map(.x = post_dist, ~ {
          unname(Rfast::colsums(.x))
        })

      # normalize the updated prior frequency to obtain prior density function
      prior_dense <- purrr::map2(.x = prior_freq2, .y = nstd, ~ {
        .x / .y
      })

      # divide the prior densities into the reference and free groups
      prior_dense_free <- prior_dense[free.group]

      # compute the mean and sd of the updated prior distribution for the free groups
      moments_free <- purrr::map(.x = prior_dense_free, ~ {
        cal_moment(node = quadpt, weight = .x)
      })
      weights.free <- purrr::map(.x = moments_free, ~ {
        gen.weight(dist = "norm", mu = .x[1], sigma = sqrt(.x[2]), theta = quadpt)
      })

      # replace the old densities with the new densities for the free groups
      weights[free.group] <- weights.free
    }
  }

  ## ---------------------------------------------------------------
  # compute the sum of loglikelihood values
  if (is.null(idx.std)) {
    llike <- sum(log(likehd %*% matrix(weights[, 2])))
  } else {
    likehd.gr <- purrr::map(.x = idx.std, ~ {
      likehd[.x, ]
    })
    llike <- purrr::map2(
      .x = likehd.gr, .y = weights,
      .f = ~ {
        sum(log(.x %*% matrix(.y[, 2])))
      }
    )
  }

  # update the item parameters in the elm_item object
  elm_item$pars <- par_df

  # organize the the results
  rst <- list(
    elm_item = elm_item, convergence = convergence, noconv_items = noconv_items,
    weights = weights, loglike = llike
  )

  # return the results
  rst
}
