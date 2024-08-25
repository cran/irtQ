# This function computes the information matrix of the item parameter estimates using the cross-product method
#' @importFrom Rfast rowsums
info_xpd <- function(elm_item, freq.cat, post_dist, quadpt.vec, n.quadpt.vec, nstd,
                     D = 1, loc_1p_const, loc_else, n.1PLM, fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl,
                     a.val.gpcm, g.val, reloc.par) {
  # extract model and cats object
  cats <- elm_item$cats
  model <- elm_item$model

  # create an empty object to contain the kernel of fisher identity equation
  # across all item response patterns
  kernel_fisher <- NULL

  # create an empty list to contain the gradient matrix of joint log likelihood functions
  # across all item response patterns
  grad_list <- list()

  # the dichotomous items: 1PLM with constrained slope values
  if (!is.null(loc_1p_const)) {
    # prepare input files to estimate the 1PLM item parameters
    f_i <- r_i <- array(0, c(n.quadpt.vec, n.1PLM))
    for (k in 1:n.1PLM) {
      f_i[, k] <- Rfast::rowsums(freq.cat[loc_1p_const][[k]])
      r_i[, k] <- freq.cat[loc_1p_const][[k]][, 2]
    }

    # use the final item parameter estimates
    item_par <-
      set_startval(
        pars = elm_item$pars, item = loc_1p_const,
        use.startval = TRUE, mod = "1PLM",
        score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
        fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
      )

    # compute the gradient vectors
    grad <-
      unname(grad_llike(
        item_par = item_par, f_i = f_i, r_i = r_i, theta = quadpt.vec,
        n.theta = n.quadpt.vec, mod = "1PLM", D = D, fix.a.1pl = fix.a.1pl,
        n.1PLM = n.1PLM
      ))

    # add the gradient matrix
    grad_list <- c(grad_list, list(grad))
  }

  # all other items
  if (length(loc_else) >= 1) {
    for (i in 1:length(loc_else)) {
      # prepare information to estimate item parameters
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a dichotomous item
      if (score.cat == 2) {
        f_i <- Rfast::rowsums(freq.cat[loc_else][[i]])
        r_i <- freq.cat[loc_else][[i]][, 2]

        # use the final item parameter estimates
        item_par <-
          set_startval(
            pars = elm_item$pars, item = loc_else[i],
            use.startval = TRUE, mod = mod,
            score.cat = score.cat, fix.a.1pl = TRUE, fix.g = fix.g,
            fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
          )

        # compute the gradient vectors
        grad <-
          unname(grad_llike(
            item_par = item_par, f_i = f_i, r_i = r_i, theta = quadpt.vec,
            n.theta = n.quadpt.vec, mod = mod, D = D,
            fix.a.1pl = ifelse(mod == "1PLM", TRUE, FALSE),
            fix.g = fix.g, a.val.1pl = a.val.1pl, g.val = .2, n.1PLM = NULL
          ))

        # add the gradient matrix
        grad_list <- c(grad_list, list(grad))
      }

      # in case of a polytomous item
      if (score.cat > 2) {
        r_i <- freq.cat[loc_else][[i]]
        # r_i <- cbind(rep(1, (n.quadpt.vec/nstd))) %x% r_i
        # r_i <- r_i[rep(1:nrow(r_i), times = (n.quadpt.vec/nstd)),   ]
        r_i <- do.call(
          "rbind",
          replicate(
            n = (n.quadpt.vec / nstd),
            expr = r_i, simplify = FALSE
          )
        )

        # use the final item parameter estimates
        item_par <-
          set_startval(
            pars = elm_item$pars, item = loc_else[i],
            use.startval = TRUE, mod = mod,
            score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
            fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
          )

        # compute the gradient vectors
        grad <-
          unname(grad_llike(
            item_par = item_par, r_i = r_i, theta = quadpt.vec, n.theta = n.quadpt.vec,
            mod = mod, D = 1, fix.a.gpcm = ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
            a.val.gpcm = a.val.gpcm, n.1PLM = NULL
          ))

        # add the gradient matrix
        grad_list <- c(grad_list, list(grad))
      }
    }
  }

  # delete 'freq.cat' object
  rm(freq.cat, envir = environment(), inherits = FALSE)

  # cbind the grad matrices
  grad_mat <- do.call("cbind", grad_list)

  # delete 'grad_list' object
  rm(grad_list, envir = environment(), inherits = FALSE)

  # compute the kernel_fisher matrix
  kernel_fisher <- split.data.frame(cbind(grad_mat * c(post_dist)), quadpt.vec)
  kernel_fisher <- Reduce(f = "+", x = kernel_fisher)

  # relocate the columns of matrix to locate the standard errors on the original position of item parameters
  kernel_fisher <- kernel_fisher[, order(reloc.par)]

  # compute the observed information matrix using the cross-product methods
  info_mat <- 0L
  for (i in 1:nstd) {
    info_mat <- info_mat + crossprod(x = rbind(kernel_fisher[i, ]))
  }

  # return the results
  info_mat
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
