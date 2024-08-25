# This function is an internal function used in the 'est_item' function.
estimation1 <- function(f_i, r_i, s_i, theta, mod = c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), score.cat, D = 1, nstd,
                        fix.a.1pl = TRUE, fix.a.gpcm = FALSE, fix.g = FALSE, a.val.1pl = 1, a.val.gpcm = 1, g.val = .2, n.1PLM = NULL,
                        aprior = list(dist = "lnorm", params = c(1, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        gprior = list(dist = "beta", params = c(5, 17)),
                        use.aprior = FALSE, use.bprior = FALSE, use.gprior = TRUE,
                        control, startval = NULL, lower, upper) {
  # estimation of item parameter & standard error
  if (!fix.a.1pl & mod == "1PLM") {
    # estimate the item parameters
    est <-
      suppressWarnings(
        tryCatch(
          {
            stats::nlminb(startval,
              objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
              nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = n.1PLM,
              aprior = aprior, bprior = bprior, gprior = gprior,
              use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
              gradient = grad_item_drm,
              hessian = hess_item_drm,
              control = control, lower = lower, upper = upper
            )
          },
          error = function(e) {
            NULL
          }
        )
      )
    if (is.null(est) || est$convergence > 0L) {
      # if error or non-convergence, only use the gradient
      est <- stats::nlminb(startval,
        objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = n.1PLM,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        gradient = grad_item_drm,
        # hessian=hess_item_drm,
        control = control, lower = lower, upper = upper
      )
    }

    # estimate the standard error of estimates
    hess <- hess_item_drm(est$par,
      f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D, nstd = nstd,
      fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = n.1PLM,
      aprior = aprior, bprior = bprior, gprior = gprior,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      adjust = FALSE
    )
  } else {
    # when the item slope parameters are not constrained to be across all items
    if (mod %in% c("1PLM", "2PLM", "3PLM")) {
      # initial estimation to find better starting values
      tmp_est <- vector("list", 3)
      tmp_est[[1]] <- stats::nlminb(startval,
        objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        gradient = grad_item_drm,
        control = list(eval.max = 100, iter.max = 50, trace = 0, step.min = 0.1, step.max = 1), lower = lower, upper = upper
      )
      tmp_est[[2]] <- stats::nlminb(startval,
        objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        gradient = grad_item_drm,
        control = list(eval.max = 100, iter.max = 50, trace = 0, step.min = 1, step.max = 2), lower = lower, upper = upper
      )
      tmp_est[[3]] <- stats::nlminb(startval,
        objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        gradient = grad_item_drm,
        control = list(eval.max = 100, iter.max = 50, trace = 0, step.min = 2, step.max = 3), lower = lower, upper = upper
      )
      tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))
      startval <- tmp_est[[tmp_num]]$par

      # estimate the item parameters
      est <-
        suppressWarnings(
          tryCatch(
            {
              stats::nlminb(startval,
                objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
                nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
                aprior = aprior, bprior = bprior, gprior = gprior,
                use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
                gradient = grad_item_drm,
                hessian = hess_item_drm,
                control = control, lower = lower, upper = upper
              )
            },
            error = function(e) {
              NULL
            }
          )
        )

      # second estimation with second alternative values when the first calibration fails
      if (is.null(est) || est$convergence > 0L) {
        tmp_num <- order(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))[2]
        startval <- tmp_est[[tmp_num]]$par
        est <-
          suppressWarnings(
            tryCatch(
              {
                stats::nlminb(startval,
                  objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
                  nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
                  aprior = aprior, bprior = bprior, gprior = gprior,
                  use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
                  gradient = grad_item_drm,
                  hessian = hess_item_drm,
                  control = control, lower = lower, upper = upper
                )
              },
              error = function(e) {
                NULL
              }
            )
          )
      }

      # third estimation by only using the gradient when the second calibration fails
      if (is.null(est) || est$convergence > 0L) {
        tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))
        startval <- tmp_est[[tmp_num]]$par

        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(startval,
          objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
          nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
          aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
          gradient = grad_item_drm,
          # hessian=hess_item_drm,
          control = control, lower = lower, upper = upper
        )
      }

      # estimate the standard error of estimates
      hess <- hess_item_drm(est$par,
        f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        adjust = FALSE
      )
    } else {
      # initial estimation to find better starting values
      tmp_est <- stats::nlminb(startval,
        objective = loglike_prm, r_i = r_i, theta = theta, pr.mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
        aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
        gradient = grad_item_prm,
        control = list(eval.max = 50, iter.max = 30, step.min = 0.1, step.max = 1, trace = 0), lower = lower, upper = upper
      )
      startval <- tmp_est$par

      # estimate the item parameters
      est <-
        suppressWarnings(
          tryCatch(
            {
              stats::nlminb(startval,
                objective = loglike_prm, r_i = r_i, theta = theta, pr.mod = mod, D = D,
                nstd = nstd, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
                aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
                gradient = grad_item_prm,
                hessian = hess_item_prm,
                control = control, lower = lower, upper = upper
              )
            },
            error = function(e) {
              NULL
            }
          )
        )

      # second estimation by only using the gradient when the first calibration fails
      if (is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(startval,
          objective = loglike_prm, r_i = r_i, theta = theta, pr.mod = mod, D = D,
          nstd = nstd, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
          aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
          gradient = grad_item_prm,
          # hessian=hess_item_prm,
          control = control, lower = lower, upper = upper
        )
      }

      # estimate the standard error of estimates
      hess <- hess_item_prm(est$par,
        r_i = r_i, theta = theta, pr.mod = mod, D = D,
        nstd = nstd, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
        aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
        adjust = FALSE
      )
    }
  }

  # compute the variance-covariance matrix of the item parameter estimates, and
  # check if the hessian matrix can be inversed
  cov_mat <- suppressWarnings(tryCatch(
    {
      solve(hess, tol = 1e-200)
    },
    error = function(e) {
      NULL
    }
  ))

  # compute the standard errors of item parameter estimates
  # also, if the hessian matrix does not have an inverse matrix, then create a matrix with NaNs
  if (is.null(cov_mat)) {
    se <- rep(99999, length(diag(hess)))
    cov_mat <- array(NaN, c(length(startval), length(startval)))
  } else {
    se <- suppressWarnings(sqrt(diag(cov_mat)))
  }

  # prevent showing NaN values of standard errors
  if (any(is.nan(se))) {
    se[is.nan(se)] <- 99999
  }

  # set an upper bound of standard error
  se <- ifelse(se > 99999, 99999, se)

  # return results
  rst <- list(pars = est$par, se = se, covariance = cov_mat, convergence = est$convergence, objective = est$objective)
  rst
}


# This function is an internal function used in the 'est_irt' function.
estimation2 <- function(f_i, r_i, s_i, quadpt, mod = c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), score.cat, D = 1, n.quad,
                        fix.a.1pl = TRUE, fix.a.gpcm = FALSE, fix.g = FALSE, a.val.1pl = 1, a.val.gpcm = 1, g.val = .2, n.1PLM = NULL,
                        aprior = list(dist = "lnorm", params = c(1, 0.5)),
                        bprior = list(dist = "norm", params = c(0.0, 1.0)),
                        gprior = list(dist = "beta", params = c(5, 17)),
                        use.aprior = FALSE, use.bprior = FALSE, use.gprior = TRUE,
                        control, startval = NULL, lower, upper, iter = NULL) {
  # estimation of item parameter & standard error
  if (!fix.a.1pl & mod == "1PLM") {
    # when the item slope parameters are constrained to be equal across all 1PLM items
    # estimate the item parameters
    est <-
      suppressWarnings(
        tryCatch(
          {
            stats::nlminb(startval,
              objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt,
              mod = mod, D = D, nstd = n.quad, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl,
              g.val = g.val, n.1PLM = n.1PLM, aprior = aprior, bprior = bprior, gprior = gprior,
              use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
              gradient = grad_item_drm, hessian = hess_item_drm,
              control = control, lower = lower, upper = upper
            )
          },
          error = function(e) {
            NULL
          }
        )
      )
    if (is.null(est) || est$convergence > 0L) {
      # if error or non-convergence, only use the gradient
      est <- stats::nlminb(startval,
        objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt,
        mod = mod, D = D, nstd = n.quad, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl,
        g.val = g.val, n.1PLM = n.1PLM, aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        gradient = grad_item_drm,
        # hessian=hess_item_drm,
        control = control, lower = lower, upper = upper
      )
    }
  } else {
    # when the item slope parameters are not constrained to be across all items
    if (mod %in% c("1PLM", "2PLM", "3PLM")) {
      # initial estimation to find better starting values
      if (iter < 3) {
        tmp_est <- vector("list", 3)
        tmp_est[[1]] <- stats::nlminb(startval,
          objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt, mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
          aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
          gradient = grad_item_drm,
          control = list(eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5), trace = 0, step.min = 0.1, step.max = 1),
          lower = lower, upper = upper
        )
        tmp_est[[2]] <- stats::nlminb(startval,
          objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt, mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
          aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
          gradient = grad_item_drm,
          control = list(eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5), trace = 0, step.min = 1, step.max = 2),
          lower = lower, upper = upper
        )
        tmp_est[[3]] <- stats::nlminb(startval,
          objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt, mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
          aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
          gradient = grad_item_drm,
          control = list(eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5), trace = 0, step.min = 2, step.max = 3),
          lower = lower, upper = upper
        )
        tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))
        startval <- tmp_est[[tmp_num]]$par
      }

      # estimate the item parameters
      est <- tryCatch(
        {
          stats::nlminb(startval,
            objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt, mod = mod, D = D,
            nstd = n.quad, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
            aprior = aprior, bprior = bprior, gprior = gprior,
            use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
            gradient = grad_item_drm,
            hessian = hess_item_drm,
            control = control, lower = lower, upper = upper
          )
        },
        error = function(e) {
          NULL
        }
      )
      if (is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(startval,
          objective = loglike_drm, f_i = f_i, r_i = r_i, s_i = s_i, theta = quadpt, mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.1pl, fix.g = fix.g, a.val = a.val.1pl, g.val = g.val, n.1PLM = NULL,
          aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
          gradient = grad_item_drm,
          # hessian=hess_item_drm,
          control = control, lower = lower, upper = upper
        )
      }
    } else {
      # initial estimation to find better starting values
      if (iter < 3) {
        tmp_est <- vector("list", 3)
        tmp_est[[1]] <- stats::nlminb(startval,
          objective = loglike_prm, r_i = r_i, theta = quadpt, pr.mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
          aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
          gradient = grad_item_prm,
          control = list(eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5), step.min = 0.1, step.max = 1, trace = 0),
          lower = lower, upper = upper
        )
        tmp_est[[2]] <- stats::nlminb(startval,
          objective = loglike_prm, r_i = r_i, theta = quadpt, pr.mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
          aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
          gradient = grad_item_prm,
          control = list(eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5), step.min = 1, step.max = 2, trace = 0),
          lower = lower, upper = upper
        )
        tmp_est[[3]] <- stats::nlminb(startval,
          objective = loglike_prm, r_i = r_i, theta = quadpt, pr.mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
          aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
          gradient = grad_item_prm,
          control = list(eval.max = 40 - (iter * 10), iter.max = 25 - (iter * 5), step.min = 2, step.max = 3, trace = 0),
          lower = lower, upper = upper
        )
        tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective, tmp_est[[3]]$objective))
        startval <- tmp_est[[tmp_num]]$par
      }

      # estimate the item parameters
      est <- tryCatch(
        {
          stats::nlminb(startval,
            objective = loglike_prm, r_i = r_i, theta = quadpt, pr.mod = mod, D = D,
            nstd = n.quad, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
            aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
            gradient = grad_item_prm,
            hessian = hess_item_prm,
            control = control, lower = lower, upper = upper
          )
        },
        error = function(e) {
          NULL
        }
      )
      if (is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(startval,
          objective = loglike_prm, r_i = r_i, theta = quadpt, pr.mod = mod, D = D,
          nstd = n.quad, fix.a = fix.a.gpcm, a.val = a.val.gpcm,
          aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior,
          gradient = grad_item_prm,
          # hessian=hess_item_prm,
          control = control, lower = lower, upper = upper
        )
      }
    }
  }

  # return results
  rst <- list(pars = est$par, convergence = est$convergence, objective = est$objective)

  rst
}
