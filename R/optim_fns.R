# Factories that build the (objective, gradient, hessian) trio for
# nlminb-based item parameter estimation. The factory captures the data
# arguments once and returns three closures that share a single P(theta)
# cache, so each nlminb evaluation point computes drm()/prm() at most
# once instead of three times (once per closure). The cache is bit-exact
# — the cached P matrix is the unmodified return value of drm()/prm()
# for the active model branch.
#
# Cache hits in practice. nlminb (PORT) typically evaluates objective,
# then gradient, then hessian at the same accepted iterate, with line
# search trials only invoking objective. Whenever the most recent
# objective call landed on the iterate that gradient/hessian then
# evaluate, the cache hits and saves two drm()/prm() calls per
# accepted point. Line search rejections fill the cache with the
# rejected trial; the next gradient/hessian call simply misses,
# refreshes the cache, and proceeds — no correctness risk.

# Build the (objective, gradient, hessian) closure trio for dichotomous
# (1PLM/2PLM/3PLM) calibration with nlminb. All static arguments
# (response counts, theta grid, model flags, priors) are captured once
# in the factory's environment; each closure takes only `item_par` and
# delegates to loglike_drm / grad_item_drm / hess_item_drm with the
# shared cached P(theta) injected via `p_cache`.
make_drm_optim_fns <- function(f_i, r_i, s_i, theta, mod, D, nstd,
                               fix.a, fix.g, a.val, g.val, n.1PLM,
                               aprior, bprior, gprior,
                               use.aprior, use.bprior, use.gprior) {
  # cache state: last item_par for which P was computed and the
  # corresponding probability matrix returned by drm(). identical()
  # comparison is bit-exact, so a cache hit guarantees the returned
  # P is the same object loglike/grad/hess would have computed.
  cache_par <- NULL
  cache_p <- NULL

  # resolve the active branch's (a, b, g) from item_par exactly as
  # loglike_drm / grad_item_drm / hess_item_drm_inner do internally,
  # then call drm() once. Keeping this branch logic in one place avoids
  # subtle drift if the model dispatch in those functions changes.
  get_p <- function(item_par) {
    # cache hit: same item_par as the previous call
    if (!is.null(cache_par) && identical(cache_par, item_par)) {
      return(cache_p)
    }

    # (1) 1PLM with slope constrained equal across items
    if (!fix.a & mod == "1PLM") {
      a <- rep(item_par[1], n.1PLM)
      b <- item_par[-1]
      g <- 0
    }

    # (2) 1PLM with slope fixed to a constant (a.val)
    if (fix.a & mod == "1PLM") {
      a <- a.val
      b <- item_par
      g <- 0
    }

    # (3) 2PLM
    if (mod == "2PLM") {
      a <- item_par[1]
      b <- item_par[2]
      g <- 0
    }

    # (4) 3PLM with guessing estimated
    if (!fix.g & mod == "3PLM") {
      a <- item_par[1]
      b <- item_par[2]
      g <- item_par[3]
    }

    # (5) 3PLM with guessing fixed to a constant (g.val)
    if (fix.g & mod == "3PLM") {
      a <- item_par[1]
      b <- item_par[2]
      g <- g.val
    }

    # compute P(theta) once and store in the closure's environment so
    # the next gradient/hessian call at the same point gets a hit
    p <- drm(theta = theta, a = a, b = b, g = g, D = D)
    cache_par <<- item_par
    cache_p <<- p
    p
  }

  # the trio: each closure asks get_p() for the cached P, then forwards
  # all original arguments plus p_cache to the underlying function. The
  # `hessian` closure routes through hess_item_drm() so the singularity
  # adjustment loop (adjust = TRUE) keeps its current behavior inside
  # nlminb — the SE-only hess_item_drm() call sites elsewhere remain
  # unchanged because p_cache defaults to NULL.
  list(
    objective = function(item_par) {
      p_cache <- get_p(item_par)
      loglike_drm(
        item_par = item_par, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta,
        mod = mod, D = D, nstd = nstd,
        fix.a = fix.a, fix.g = fix.g, a.val = a.val, g.val = g.val, n.1PLM = n.1PLM,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        p_cache = p_cache
      )
    },
    gradient = function(item_par) {
      p_cache <- get_p(item_par)
      grad_item_drm(
        item_par = item_par, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta,
        mod = mod, D = D, nstd = nstd,
        fix.a = fix.a, fix.g = fix.g, a.val = a.val, g.val = g.val, n.1PLM = n.1PLM,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        p_cache = p_cache
      )
    },
    hessian = function(item_par) {
      p_cache <- get_p(item_par)
      hess_item_drm(
        item_par = item_par, f_i = f_i, r_i = r_i, s_i = s_i, theta = theta,
        mod = mod, D = D, nstd = nstd,
        fix.a = fix.a, fix.g = fix.g, a.val = a.val, g.val = g.val, n.1PLM = n.1PLM,
        aprior = aprior, bprior = bprior, gprior = gprior,
        use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        adjust = TRUE,
        p_cache = p_cache
      )
    }
  )
}


# Build the (objective, gradient, hessian) closure trio for polytomous
# (GRM / GPCM / PCM) calibration with nlminb. The factory captures all
# static data once; each closure shares a single prob_cache list so the
# expensive probability computation runs at most once per nlminb point.
#
# Cache structure per model:
#   GRM  → list(allPst, P)             — boundary and category probs
#   GPCM → list(theta_d, numer, denom, P) — full intermediate tensors
#                                           needed by both grad and hess
#
# The (theta_d, numer, denom) tuple is cached for GPCM because
# grad_item_prm and hess_item_prm_inner both use all three quantities
# beyond P; recomputing only P would still duplicate the Outer + cumsum
# + exp work.
make_prm_optim_fns <- function(r_i, theta, pr.mod, D, nstd,
                               fix.a, a.val,
                               aprior, bprior,
                               use.aprior, use.bprior) {
  # cache state: last item_par for which prob was computed
  cache_par <- NULL
  cache_prob <- NULL

  # compute the full probability objects for the active model and branch
  get_prob <- function(item_par) {
    # cache hit: return without any computation
    if (!is.null(cache_par) && identical(cache_par, item_par)) {
      return(cache_prob)
    }

    # determine a and d from item_par exactly as loglike/grad/hess do
    if (!fix.a) {
      a <- item_par[1]
      d <- item_par[-1]
    } else {
      a <- a.val
      d <- item_par
    }

    # (A) GRM: boundary probabilities from drm(), category probs derived
    if (pr.mod == "GRM") {
      m <- length(d)
      allPst <- drm(theta = theta, a = rep(a, m), b = d, g = 0, D = D)
      P <- cbind(1, allPst) - cbind(allPst, 0)
      P[P > 9999999999e-10] <- 9999999999e-10
      P[P < 1e-10] <- 1e-10
      prob <- list(allPst = allPst, P = P)
    }

    # (B) GPCM and PCM (fix.a): numerator / denominator tensors + P.
    # For PCM `d` is item_par itself and a = a.val.
    # The leading 0 is prepended to d exactly as the function bodies do.
    if (pr.mod == "GPCM") {
      d_full <- c(0, d)
      Da <- D * a
      theta_d <- Rfast::Outer(x = theta, y = d_full, oper = "-")
      z <- Da * theta_d
      cumsum_z <- t(Rfast::colCumSums(z))
      if (any(cumsum_z > 700)) {
        cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
      }
      if (any(cumsum_z < -700)) {
        cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
      }
      numer <- exp(cumsum_z)
      denom <- Rfast::rowsums(numer)
      P <- numer / denom
      prob <- list(theta_d = theta_d, numer = numer, denom = denom, P = P)
    }

    # store for the next call at this item_par
    cache_par <<- item_par
    cache_prob <<- prob
    prob
  }

  # the trio: each closure retrieves the cached prob list and forwards it
  # to the corresponding function via prob_cache. The hessian closure
  # routes through hess_item_prm() so the singularity-adjust loop keeps
  # its current behavior. SE-only hess_item_prm() call sites outside the
  # factory remain unchanged (prob_cache defaults to NULL there).
  list(
    objective = function(item_par) {
      prob_cache <- get_prob(item_par)
      loglike_prm(
        item_par = item_par, r_i = r_i, theta = theta, pr.mod = pr.mod, D = D, nstd = nstd,
        fix.a = fix.a, a.val = a.val,
        aprior = aprior, bprior = bprior,
        use.aprior = use.aprior, use.bprior = use.bprior,
        prob_cache = prob_cache
      )
    },
    gradient = function(item_par) {
      prob_cache <- get_prob(item_par)
      grad_item_prm(
        item_par = item_par, r_i = r_i, theta = theta, pr.mod = pr.mod, D = D, nstd = nstd,
        fix.a = fix.a, a.val = a.val,
        aprior = aprior, bprior = bprior,
        use.aprior = use.aprior, use.bprior = use.bprior,
        prob_cache = prob_cache
      )
    },
    hessian = function(item_par) {
      prob_cache <- get_prob(item_par)
      hess_item_prm(
        item_par = item_par, r_i = r_i, theta = theta, pr.mod = pr.mod, D = D, nstd = nstd,
        fix.a = fix.a, a.val = a.val,
        aprior = aprior, bprior = bprior,
        use.aprior = use.aprior, use.bprior = use.bprior,
        adjust = TRUE,
        prob_cache = prob_cache
      )
    }
  )
}
