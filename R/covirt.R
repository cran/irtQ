#' Asymptotic Variance-Covariance Matrices of Item Parameter Estimates
#'
#' This function computes the analytical asymptotic variance-covariance matrices
#' of item parameter estimates for dichotomous and polytomous IRT models,
#' without requiring examinee response data. Given a set of item parameter
#' estimates and the corresponding sample sizes, the function derives the
#' matrices using analytical formulas (e.g., Li & Lissitz, 2004; Thissen &
#' Wainer, 1982). The square roots of the diagonal elements (variances) provide
#' the asymptotic standard errors of the maximum likelihood estimates.
#'
#' @inheritParams catsib
#' @param nstd An integer or a vector of integers indicating the sample size(s).
#'   If a vector is provided, its length must match the number of items in the
#'   `x` argument. Default is 1,000. See **Details**.
#' @param pcm.loc A vector of integers indicating the positions of items
#'   calibrated using the partial credit model (PCM). For PCM items, the
#'   variance-covariance matrices are computed only for the item category
#'   difficulty parameters. Default is `NULL`. See the **Details** for more
#'   information.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal prior distribution. These values are used
#'   to generate the Gaussian quadrature points and weights when `weights = NULL`.
#'   Default is `c(0, 1)`.
#' @param nquad An integer indicating the number of Gaussian quadrature points
#'   to be generated from the normal prior distribution. The specified value is
#'   used when `weights` is not `NULL`. Default is 41.
#' @param weights An optional two-column data frame or matrix where the first
#'   column is the quadrature points (nodes) and the second column is the
#'   corresponding weights. This is typically used in quadrature-based IRT
#'   analysis.
#'
#' @details The standard errors obtained from this analytical approach are
#' generally considered lower bounds of the true standard errors (Thissen &
#' Wainer, 1982). Thus, they may serve as useful approximations for evaluating
#' the precision of item parameter estimates when empirical standard errors are
#' not reported in the literature or research reports.
#'
#' If the item parameters provided in the `x` argument were calibrated using
#' different sample sizes, a corresponding vector of sample sizes must be
#' specified via the `nstd` argument. For example, suppose you wish to compute
#' the variance-covariance matrices of five 3PLM items that were calibrated
#' using 500, 600, 1,000, 2,000, and 700 examinees, respectively. In this case,
#' set `nstd = c(500, 600, 1000, 2000, 700)`.
#'
#' Since the item metadata allows only `"GPCM"` to denote both the partial
#' credit model (PCM) and the generalized partial credit model (GPCM), PCM items
#' must be explicitly identified using the `pcm.loc` argument. This is necessary
#' because the category difficulty parameters of PCM items require special
#' handling when computing variance-covariance matrices. For instance, if you
#' wish to compute the matrices for five polytomous items and the last two were
#' calibrated using PCM, then specify `pcm.loc = c(4, 5)`.
#'
#' @return A named list with the following two components:
#'
#' - `cov`: A named list of variance-covariance matrices for item parameter estimates.
#' Each element corresponds to a single item and contains a square matrix whose
#' dimensions match the number of estimated parameters for that item. For
#' dichotomous items, this typically includes slopes and intercepts. For
#' polytomous items, it includes category difficulty parameters (for PCM) or
#' both slope and difficulty (or threshold) parameters (for GRM and GPCM).
#'
#' - `se`: A named list of vectors containing the asymptotic standard errors (ASEs)
#' of the item parameter estimates, computed as the square roots of the diagonal
#' elements of each corresponding variance-covariance matrix in `cov`.
#'
#' The names of the list elements in both `cov` and `se` correspond to the item
#' identifiers (e.g., item names or labels) as given in the first column of the
#' input `x`.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Li, Y. & Lissitz, R. (2004). Applications of the analytically
#'   derived asymptotic standard errors of item response theory item parameter
#'   estimates. *Journal of educational measurement, 41*(2), 85-117.
#'
#'   Thissen, D. & Wainer, H. (1982). Weighted likelihood estimation of ability
#'   in item response theory. *Psychometrika, 54*(3), 427-450.
#'
#' @seealso [irtQ::est_irt()], [irtQ::simdat()], [irtQ::shape_df()],
#'   [irtQ::gen.weight()]
#'
#' @export
#'
#' @examples
#' # Example using a "-prm.txt" file exported from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select the first two dichotomous items and the last polytomous item
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # Compute the variance-covariance matrices assuming a sample size of 2,000
#' covirt(x, D = 1, nstd = 2000, norm.prior = c(0, 1), nquad = 41)
#'
covirt <- function(x, D = 1, nstd = 1000, pcm.loc = NULL, norm.prior = c(0, 1), nquad = 41, weights = NULL) {
  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # break down the item metadata into several components
  elm_list <-
    purrr::map(.x = 1:nrow(x), ~ {
      breakdown(x[.x, ])
    })

  # specify the locations of PCM items
  pcm.a.logic <- rep(FALSE, nrow(x))
  if (!is.null(pcm.loc)) {
    pcm.a.logic[pcm.loc] <- TRUE
  }

  # generate quadrature points and weights
  if (is.null(weights)) {
    wts <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
  } else {
    wts <- data.frame(weights)
  }

  # sample size
  if (length(nstd) == 1L) nstd <- rep(nstd, nrow(x))
  if (length(nstd) != nrow(x)) {
    stop("Sample size must be an integer value or a vector with a length of total items.", call. = FALSE)
  }

  # compute the var-covariance matrix
  cov_par <- purrr::pmap(
    .l = list(x = elm_list, y = nstd, z = pcm.a.logic),
    .f = function(x, y, z) {
      cov_mat(elm_item = x, D = D, nstd = y, wts = wts, pcm.a = z)
    }
  )
  se_par <- purrr::map(cov_par, .f = function(x) sqrt(diag(x)))
  names(cov_par) <- x[, 1]
  names(se_par) <- x[, 1]

  # return results
  rst <- list(cov = cov_par, se = se_par)
  rst
}


# This function computes the analytical var-covariance matrix of item parameter estimates for an item
cov_mat <- function(elm_item, D = 1, nstd = 1000, wts, pcm.a = FALSE) {
  info_mat <- integrand(elm_item, theta = wts[, 1], dens = wts[, 2], D = D, pcm.a = pcm.a)

  rst <- tryCatch(
    {
      solve(info_mat * nstd)
    },
    error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    }
  )
  if (is.null(rst)) {
    rst <- matrix(NA, nrow = dim(info_mat)[1], ncol = dim(info_mat)[2])
  }

  rst
}


# This function computes the kernel of integration when computing
# the analytical var-covariance matrix of item parameter estimates
#' @importFrom stats na.exclude
#' @importFrom Rfast rowprods
integrand <- function(elm_item, theta, dens, D = 1, pcm.a = FALSE) {
  # assign items into DRM and PRM groups
  idx.all <- idxfinder(elm_item)
  idx.drm <- idx.all$idx.drm
  idx.prm <- idx.all$idx.prm

  # For a DRM item
  if (!is.null(idx.drm)) {
    # extract required information
    item_par <- c(
      elm_item$par[idx.drm, 1],
      elm_item$par[idx.drm, 2],
      elm_item$par[idx.drm, 2]
    )
    cats <- elm_item$cats[idx.drm]
    model <- elm_item$model[idx.drm]

    # compute a product of the probabilities p and q
    pq <- Rfast::rowprods(trace(elm_item = elm_item, theta = theta, D = D, tcc = FALSE)$prob.cats[[1]])

    # prevent that pq has zero probability
    pq[pq == 0L] <- 1e-20

    # create the gradient for score category probability
    funList <- equation_scocat(model = model, cats = NULL, hessian = FALSE, type = "item")

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for (i in 1:3) {
      args.pars[[i]] <- item_par[i]
    }
    args.list <- args.pars
    args.list$theta <- theta
    args.list$D <- D

    # compute the score category information
    # select a function for each score category
    fun.tmp <- funList[[1]]

    # implement the function for each score category
    tmp <- do.call("fun.tmp", args.list, envir = environment())

    # compute the kernel of integration
    p.d1 <- attributes(tmp)$gradient
    info.list <- purrr::map(
      .x = 1:length(theta),
      .f = function(i) (1 / pq[i]) * outer(p.d1[i, ], p.d1[i, ]) * dens[i]
    )
    info.mat <- Reduce(f = "+", x = info.list)
  }

  # For a PRM item
  if (!is.null(idx.prm)) {
    # extract required information
    cats <- elm_item$cats[idx.prm]
    model <- elm_item$model[idx.prm]
    item_par <- stats::na.exclude(elm_item$par[idx.prm, ])

    # compute a product of the probabilities p and q
    ps <- trace(elm_item = elm_item, theta = theta, D = D, tcc = FALSE)$prob.cats[[1]]
    pqs <- ps * (1 - ps)

    # prevent that pq has zero probability
    pqs[pqs == 0L] <- 1e-20

    # create the gradient for score category probability
    funList <- equation_scocat(model = model, cats = cats, fix.a.gpcm = pcm.a, hessian = FALSE, type = "item")

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for (i in 1:cats) {
      args.pars[[i]] <- item_par[i]
    }
    args.list <- args.pars
    args.list$theta <- theta
    args.list$D <- D

    # compute the score category information
    info.mat <- vector("list", cats)
    for (i in 1:cats) {
      # select a function for each score category
      fun.tmp <- funList[[i]]

      # implement the function for each score category
      tmp <- do.call("fun.tmp", args.list, envir = environment())

      # compute the kernel of integration
      ps.d1 <- attributes(tmp)$gradient
      info.list <- purrr::map(
        .x = 1:length(theta),
        .f = function(k) (1 / pqs[k, i]) * outer(ps.d1[k, ], ps.d1[k, ]) * dens[k]
      )
      info.mat[[i]] <- Reduce(f = "+", x = info.list)
    }

    info.mat <- Reduce(f = "+", x = info.mat)
  }

  info.mat
}
