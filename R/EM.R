# Internal E-step core shared by Estep() and Estep_fipc().  Both
# public wrappers compute the same three quantities -- a per-(
# examinee, theta) likelihood matrix, the posterior density over
# the ability quadrature, and the per-item expected frequencies
# of scored categories -- and differ only in argument naming and
# in whether the elm_item used to evaluate the likelihood is the
# same as the one returned in the result list.  Folding the
# duplicated body into this helper keeps the public signatures
# unchanged while removing ~30 lines of mirrored code.
#
# Args:
#   elm_item_likehd : list (output of breakdown()) whose $pars,
#                     $cats, $model drive the likelihood() call.
#                     For Estep() this is the only elm_item; for
#                     Estep_fipc() this is the "fixed-or-merged"
#                     elm_item (elm_item2 in the FIPC notation),
#                     used because the posterior must be evaluated
#                     against the items whose parameters are
#                     currently considered known.
#   elm_item_return : list passed back as the $elm_item slot of
#                     the result.  For Estep() this equals
#                     elm_item_likehd; for Estep_fipc() this is
#                     the "new items" elm_item (elm_item1) so
#                     downstream Mstep updates the right rows.
#   idx.drm, idx.prm : integer item-row indices into elm_item_likehd
#                     for dichotomous and polytomous items.
#   data_drm, data_prm, data_all : sparse-Matrix indicator inputs
#                     produced by divide_data().
#   weights : single-group case is a data.frame(theta, weight); the
#             multi-group case is a list of such data.frames, in
#             which case theta is read from weights[[1]][, 1] (all
#             groups share the same quadrature grid by construction).
#   D       : scaling constant (1 = logistic, 1.7 ~= normal-ogive).
#   idx.std : NULL for single group, or a list of per-group row
#             indices into the combined response matrix.
#' @importFrom Matrix crossprod
.Estep_core <- function(elm_item_likehd, elm_item_return,
                        idx.drm, idx.prm,
                        data_drm, data_prm, data_all,
                        weights, D, idx.std) {
  # quadrature points: weights[, 1] in single-group form, or
  # weights[[1]][, 1] in multi-group form (all groups share the
  # same grid in the current implementation)
  theta <- if (is.null(idx.std)) weights[, 1] else weights[[1]][, 1]

  # likelihood matrix L(theta_q | response_i) for each examinee
  # at every quadrature point; used only as the input to the
  # posterior() Bayes-rule normalization below
  likehd <- likelihood(elm_item_likehd,
    idx.drm = idx.drm, idx.prm = idx.prm,
    data_drm = data_drm, data_prm = data_prm,
    theta = theta, D = D
  )$L

  # posterior density of ability per examinee (rows sum to 1);
  # idx.std = NULL takes the single-group branch, otherwise the
  # multi-group branch with per-group prior weights
  post_dist <- posterior(likehd = likehd, weights = weights, idx.std = idx.std)

  # per-(quadrature, score-category) expected frequency = the
  # conditional expectation of the one-hot response indicator
  # given the posterior; sparse-aware crossprod exploits zero
  # entries in data_all (incl. NA-encoded missing responses)
  freq.exp <- base::as.matrix(Matrix::crossprod(post_dist, data_all))

  # the returned $elm_item is the caller-chosen one: same as the
  # likelihood-driving elm_item for Estep, the "new items" elm_item
  # for Estep_fipc
  list(
    elm_item = elm_item_return, post_dist = post_dist, freq.exp = freq.exp,
    likehd = likehd, idx.std = idx.std
  )
}

# E-step function when FIPC method is used.  Drives the likelihood
# evaluation with the "fixed (iter 1) or merged fixed+new (iter > 1)"
# elm_item (elm_item2) but returns the new-items elm_item (elm_item1)
# so the subsequent Mstep updates only the parameters being calibrated.
Estep_fipc <- function(elm_item1, elm_item2, idx.drm2, idx.prm2,
                       data_drm2, data_prm2, data_all1, weights,
                       D = 1, idx.std = NULL) {
  .Estep_core(
    elm_item_likehd = elm_item2, elm_item_return = elm_item1,
    idx.drm = idx.drm2, idx.prm = idx.prm2,
    data_drm = data_drm2, data_prm = data_prm2, data_all = data_all1,
    weights = weights, D = D, idx.std = idx.std
  )
}

# E-step function for the standard (non-FIPC) calibration path.
# Both the likelihood-driving elm_item and the returned elm_item are
# the same single argument.
Estep <- function(elm_item, idx.drm = NULL, idx.prm = NULL,
                  data_drm, data_prm, data_all, weights, D = 1, idx.std = NULL) {
  .Estep_core(
    elm_item_likehd = elm_item, elm_item_return = elm_item,
    idx.drm = idx.drm, idx.prm = idx.prm,
    data_drm = data_drm, data_prm = data_prm, data_all = data_all,
    weights = weights, D = D, idx.std = idx.std
  )
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
      # prepare input files to estimate the 1PLM item parameters.
      # cols.item$cols.1pl holds 2 * n.1PLM consecutive column indices
      # into freq.exp (one (incorrect, correct) pair per 1PLM item by
      # construction in cols4item()), so the odd positions extract the
      # "wrong" column (s_i) and the even positions extract the
      # "correct" column (r_i) for every 1PLM item simultaneously.
      s_i <- freq.exp[, cols.item$cols.1pl][, c(TRUE, FALSE), drop = FALSE]
      r_i <- freq.exp[, cols.item$cols.1pl][, c(FALSE, TRUE), drop = FALSE]
      # f_i (total responses per 1PLM item per quadrature point) is
      # the elementwise sum of s_i and r_i; this replaces the previous
      # n.1PLM-iteration for-loop that called Rfast::rowsums on each
      # 2-column slice of freq.exp -- the loop was redundant because
      # cats[k] == 2 for every 1PLM item, so the 2-col rowsum is
      # identical to a single elementwise add of the s_i/r_i matrices
      f_i <- s_i + r_i

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
          # extract the (incorrect, correct) column pair for this item;
          # cols.item$cols.all[[k]] always has length 2 for cats == 2
          cols.tmp <- cols.item$cols.all[[loc_else[i]]]
          s_i <- freq.exp[, cols.tmp[1]]
          r_i <- freq.exp[, cols.tmp[2]]
          # total responses per quadrature point = s_i + r_i; replaces
          # the previous Rfast::rowsums(freq.exp[, cols.tmp]) call,
          # which was a 2-column rowsum and is identical to the direct
          # add (avoids a function-call indirection per item per Mstep)
          f_i <- s_i + r_i

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
    # arrange the estimated item parameters into natural item order.
    # The estimation loop appends results in [loc_1p_const items,
    # then loc_else items] order, so we permute by the inverse of
    # that mapping -- which is exactly order(c(loc_1p_const,
    # loc_else)).  Replaces a 4-step copy chain (cbind a loc column
    # -> sort by it -> drop the column) that produced the same
    # final permutation but allocated three intermediate matrix
    # copies; the new path computes the integer permutation once
    # on a small length-nitem vector and applies it as a single
    # row index.
    par_df <- bind.fill(est_par, type = "rbind")
    par_df <- par_df[order(c(loc_1p_const, loc_else)), , drop = FALSE]

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
