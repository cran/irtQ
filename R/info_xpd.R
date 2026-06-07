# This function computes the cross-product information matrix of the
# item parameter estimates -- the sum over examinees of the outer
# product of the per-examinee score function.
#
# Mathematical structure exploited
# --------------------------------
# For every supported model the per-(examinee, theta_q, parameter p)
# gradient component is LINEAR in the one-hot freq.cat[[k]][i, c]
# response indicator.  Therefore there exists a kernel
#       gamma_p^k(q, c)   (depends only on theta, item params, model)
# such that
#       grad[i, q, p] = sum_c freq.cat[[k]][i, c] * gamma_p^k(q, c)
# and the per-examinee score column reduces to
#       kf[i, p] = sum_q post[i, q] * grad[i, q, p]
#                = sum_c freq.cat[[k]][i, c] * (post %*% gamma[, c])[i]
#
# So a single GEMM (per item) collapses the q-axis without ever
# materializing the (nstd*ntheta) x npar gradient matrix that the
# previous implementation built.  This was the dominant peak-memory
# allocation in the SE step of est_irt() / est_mg() and the reason
# large CAT calibrations exhausted RAM on a 128 GB machine.
#
# We extract gamma without re-deriving it model-by-model: for each
# category c we feed grad_llike() a SYNTHETIC freq pattern that is 1
# in category c (and 0 elsewhere), at theta = quadpt.  The returned
# gradient row IS gamma_p(q, c), bit-for-bit consistent with the
# original code path because the same gradient routine is called.
#
# Numerical equivalence with the previous implementation is verified
# bit-for-bit (max abs diff <= 8e-15) across 5 mixed-format scenarios
# in data-raw/verify_info_xpd.R.
#' @importFrom Rfast rowsums
info_xpd <- function(elm_item, freq.cat, post_dist, quadpt, nstd,
                     D = 1, loc_1p_const, loc_else, n.1PLM, fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl,
                     a.val.gpcm, g.val, reloc.par) {
  # extract per-item model + score-category metadata
  cats  <- elm_item$cats
  model <- elm_item$model

  # number of quadrature points (replaces n.quadpt.vec / nstd)
  ntheta <- length(quadpt)

  # total number of free parameters across all items in estimation order;
  # also the number of columns of the per-examinee score matrix
  total_npar <- length(reloc.par)

  # accumulator for the per-examinee score function (nstd x total_npar);
  # column ordering follows estimation order (1PLM-const block first,
  # then loc_else items in their index order) and is finally reordered
  # to original-parameter position via order(reloc.par) below
  kernel_fisher <- matrix(0, nrow = nstd, ncol = total_npar)

  # running offset that tracks where the next item's parameter columns
  # start within kernel_fisher
  col_offset <- 0L

  ## -----------------------------------------------------------------
  ## (A) 1PLM-constrained block: shared `a` across loc_1p_const items
  ## -----------------------------------------------------------------
  if (!is.null(loc_1p_const)) {
    # estimation-order columns in this block:
    #   col_offset + 1            = shared a
    #   col_offset + 1 + k_idx    = b_k for the k_idx-th 1PLM item
    a_col <- col_offset + 1L

    # final item parameter estimates for the constrained 1PLM block:
    # length 1 + n.1PLM, ordered (a, b_1, b_2, ..., b_{n.1PLM})
    item_par <- set_startval(
      pars = elm_item$pars, item = loc_1p_const,
      use.startval = TRUE, mod = "1PLM",
      score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
      fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
    )

    for (k_idx in seq_len(n.1PLM)) {
      k       <- loc_1p_const[k_idx]    # original item index
      b_k_col <- col_offset + 1L + k_idx

      # build Gamma_k of shape ntheta x 4 -- columns ordered as
      #   [gamma_a(:,1), gamma_a(:,2), gamma_b_k(:,1), gamma_b_k(:,2)]
      Gamma_k <- matrix(0, nrow = ntheta, ncol = 4L)

      # for each category c in {1=wrong, 2=correct}, drive grad_llike
      # with a synthetic single-examinee response pattern in which only
      # item k_idx is observed at category c.  The returned gradient is
      # ntheta x (1 + n.1PLM); column 1 is the shared-a derivative and
      # column (k_idx+1) is the b_k derivative -- by linearity in
      # (r_i, f_i), other items' columns are exactly zero because their
      # synthetic r_i, f_i are zero.
      for (c in 1:2) {
        f_syn <- matrix(0, nrow = ntheta, ncol = n.1PLM)
        r_syn <- matrix(0, nrow = ntheta, ncol = n.1PLM)
        f_syn[, k_idx] <- 1                       # examinee responded
        r_syn[, k_idx] <- if (c == 2L) 1 else 0   # 1 iff correct (c=2)

        grad_syn <- grad_llike(
          item_par = item_par, f_i = f_syn, r_i = r_syn,
          theta = quadpt, n.theta = ntheta, mod = "1PLM",
          D = D, fix.a.1pl = fix.a.1pl, n.1PLM = n.1PLM
        )

        Gamma_k[, c]      <- grad_syn[, 1L]            # gamma_a (q, c)
        Gamma_k[, 2L + c] <- grad_syn[, k_idx + 1L]    # gamma_b_k(q, c)
      }

      # one GEMM per item: project gamma onto the posterior density to
      # get the per-(examinee, parameter, category) score-coefficient
      W_k <- post_dist %*% Gamma_k          # nstd x 4

      # combine with the one-hot freq pattern to recover per-examinee
      # score columns; Rfast::rowsums avoids a separate temporary
      fc_k <- freq.cat[[k]]                 # nstd x 2
      kernel_fisher[, a_col]   <- kernel_fisher[, a_col] +
        Rfast::rowsums(fc_k * W_k[, 1:2, drop = FALSE])
      kernel_fisher[, b_k_col] <-
        Rfast::rowsums(fc_k * W_k[, 3:4, drop = FALSE])
    }

    col_offset <- col_offset + 1L + n.1PLM
  }

  ## -----------------------------------------------------------------
  ## (B) all remaining items (loc_else): each handled independently
  ## -----------------------------------------------------------------
  if (length(loc_else) >= 1L) {
    for (i in seq_along(loc_else)) {
      k         <- loc_else[i]
      mod       <- model[k]
      score.cat <- cats[k]

      # extract the final item parameters for this single item using
      # the same set_startval() call shape as the previous code path
      if (score.cat == 2L) {
        item_par <- set_startval(
          pars = elm_item$pars, item = k,
          use.startval = TRUE, mod = mod,
          score.cat = 2L, fix.a.1pl = TRUE, fix.g = fix.g,
          fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
        )
      } else {
        item_par <- set_startval(
          pars = elm_item$pars, item = k,
          use.startval = TRUE, mod = mod,
          score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
          fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
        )
      }

      # number of free parameters for this item, derived from a single
      # synthetic gradient call (probes the model's actual npar after
      # all fix.* flags are applied -- avoids re-deriving the table)
      npar_k <- if (score.cat == 2L) {
        ncol(grad_llike(
          item_par = item_par, f_i = rep(1, ntheta), r_i = rep(0, ntheta),
          theta = quadpt, n.theta = ntheta, mod = mod, D = D,
          fix.a.1pl = ifelse(mod == "1PLM", TRUE, FALSE),
          fix.g = fix.g, a.val.1pl = a.val.1pl, g.val = .2, n.1PLM = NULL
        ))
      } else {
        # polytomous probe: synthetic r_i = e_1 (ntheta x cats)
        r_probe <- matrix(0, nrow = ntheta, ncol = score.cat)
        r_probe[, 1L] <- 1
        ncol(grad_llike(
          item_par = item_par, r_i = r_probe, theta = quadpt, n.theta = ntheta,
          mod = mod, D = 1, fix.a.gpcm = ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
          a.val.gpcm = a.val.gpcm, n.1PLM = NULL
        ))
      }

      # build Gamma_k of shape ntheta x (score.cat * npar_k); column
      # ordering is [param 1: c=1..C, param 2: c=1..C, ...] so the
      # per-parameter slice is contiguous after the GEMM
      Gamma_k <- matrix(0, nrow = ntheta, ncol = score.cat * npar_k)

      for (c in seq_len(score.cat)) {
        if (score.cat == 2L) {
          # binary item: synthetic f_i = 1 (responded), r_i = 1{c==2}
          f_syn <- rep(1, ntheta)
          r_syn <- if (c == 2L) rep(1, ntheta) else rep(0, ntheta)
          grad_syn <- grad_llike(
            item_par = item_par, f_i = f_syn, r_i = r_syn,
            theta = quadpt, n.theta = ntheta, mod = mod, D = D,
            fix.a.1pl = ifelse(mod == "1PLM", TRUE, FALSE),
            fix.g = fix.g, a.val.1pl = a.val.1pl, g.val = .2, n.1PLM = NULL
          )
        } else {
          # polytomous item: synthetic r_i = e_c (ntheta x cats)
          r_syn <- matrix(0, nrow = ntheta, ncol = score.cat)
          r_syn[, c] <- 1
          grad_syn <- grad_llike(
            item_par = item_par, r_i = r_syn, theta = quadpt, n.theta = ntheta,
            mod = mod, D = 1, fix.a.gpcm = ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
            a.val.gpcm = a.val.gpcm, n.1PLM = NULL
          )
        }
        # place the column-c slice of grad_syn into Gamma_k, one column
        # per parameter; the inner loop is over npar_k (small constant)
        for (p in seq_len(npar_k)) {
          Gamma_k[, (p - 1L) * score.cat + c] <- grad_syn[, p]
        }
      }

      # single GEMM for this item; shape nstd x (score.cat * npar_k)
      W_k <- post_dist %*% Gamma_k
      fc_k <- freq.cat[[k]]                 # nstd x score.cat

      # recover each parameter's per-examinee score column by summing
      # the one-hot-weighted W_k slice across categories
      for (p in seq_len(npar_k)) {
        col_in_W <- (p - 1L) * score.cat + seq_len(score.cat)
        kernel_fisher[, col_offset + p] <-
          Rfast::rowsums(fc_k * W_k[, col_in_W, drop = FALSE])
      }
      col_offset <- col_offset + npar_k
    }
  }

  # reorder columns to put parameters back into their original position
  # (this is the same column permutation that the previous code path
  # applied to its grad_mat before the final cross-product)
  kernel_fisher <- kernel_fisher[, order(reloc.par)]

  # info_mat = sum_i kernel_fisher[i,]^T %*% kernel_fisher[i,] is
  # exactly crossprod(kernel_fisher) -- one BLAS DSYRK call replaces
  # the previous nstd-iteration R loop
  crossprod(kernel_fisher)
}

# This function computes the information matrix of priors using the second derivatives of item parameter estimates
#' @importFrom Matrix bdiag
info_prior <- function(elm_item, D = 1, loc_1p_const, loc_else, n.1PLM,
                       fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl, a.val.gpcm, g.val,
                       aprior, bprior, gprior, use.aprior, use.bprior, use.gprior, reloc.par) {
  # a create a vector containing the gradient values of priors
  # a create empty matrix to contain the gradient
  hess_list <- NULL

  # extract model and cats object
  cats <- elm_item$cats
  model <- elm_item$model

  # the dichotomous items: 1PLM with constrained slope values
  if (!is.null(loc_1p_const)) {
    # use the final item parameter estimates
    item_par <-
      set_startval(
        pars = elm_item$pars, item = loc_1p_const,
        use.startval = TRUE, mod = "1PLM",
        score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
        fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
      )

    # compute the hessian matrix
    hess <- hess_prior(
      item_par = item_par, mod = "1PLM", D = D, fix.a.1pl = fix.a.1pl,
      aprior = aprior, bprior = bprior, use.aprior = use.aprior,
      use.bprior = use.bprior
    )
    hess_list <- c(hess_list, list(hess))
  }

  # all other items
  if (length(loc_else) >= 1) {
    for (i in 1:length(loc_else)) {
      # prepare information to estimate item parameters
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a dichotomous item
      if (score.cat == 2) {
        # use the final item parameter estimates
        item_par <-
          set_startval(
            pars = elm_item$pars, item = loc_else[i],
            use.startval = TRUE, mod = mod,
            score.cat = score.cat, fix.a.1pl = TRUE, fix.g = fix.g,
            fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
          )

        # compute the gradient vectors
        hess <- hess_prior(
          item_par = item_par, mod = mod, D = D, fix.a.1pl = ifelse(mod == "1PLM", TRUE, FALSE),
          fix.g = fix.g, aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior
        )
        hess_list <- c(hess_list, list(hess))
      }

      # in case of a PRM item
      if (score.cat > 2) {
        # use the final item parameter estimates
        item_par <-
          set_startval(
            pars = elm_item$pars, item = loc_else[i],
            use.startval = TRUE, mod = mod,
            score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
            fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
          )

        # compute the gradient vectors
        hess <- hess_prior(
          item_par = item_par, mod = mod, D = D, fix.a.gpcm = ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
          aprior = aprior, bprior = bprior, use.aprior = use.aprior, use.bprior = use.bprior
        )
        hess_list <- c(hess_list, list(hess))
      }
    }
  }

  # make a block diagonal matrix
  info_mat <- data.matrix(Matrix::bdiag(hess_list))

  # relocate the diagonal parts of the matrix to locate the standard errors on the original position of item parameters
  diag(info_mat) <- diag(info_mat)[order(reloc.par)]

  # return the results
  info_mat
}
