#' Asymptotic variance-covariance matrices of item parameter estimates
#'
#' @description This function calculates the analytical asymptotic variance-covariance matrices (e.g., Li & Lissitz, 2004; Thissen & Wainer, 1982)
#' of item parameter estimates for dichotomous and polytomous IRT Models without examinee's responses to test items,
#' given a set of item parameter estimates and sample size. The square roots of variance terms in the matrices can be used as the asymptotic
#' standard errors of maximum likelihood item parameter estimates.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...).
#' See \code{\link{irtfit}}, \code{\link{info}}, or \code{\link{simdat}} for more details about the item metadata.
#' This data frame can be easily obtained using the function \code{\link{shape_df}}.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param nstd An integer value or a vector of integer values indicating a sample size. When a vector is specified, length of the vector must be
#' the same as the number of test items in the argument \code{x}. Default is 1,000. See below for details.
#' @param pcm.loc A vector of integer values indicating the locations of partial credit model (PCM) items. For the PCM items,
#' the variance-covariance matrices are computed only for the item category difficulty parameters. Default is NULL. See below for details.
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution.
#' Default is c(0,1).
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
#' @param weights A two-column matrix or data frame containing the theta values (in the first column) and the weights (in the second column)
#' for the prior distribution. The weights and theta values can be easily obtained using the function \code{\link{gen.weight}}.
#' If NULL, default values are used for the prior distribution (see the arguments of \code{norm.prior} and \code{nquad}). Default is NULL.
#'
#' @details
#' The standard errors obtained from the analytical approach are likely to represent lower bounds for the actual standard errors (Thissen & Wainer, 1982).
#' Therefore, they may be useful for assessing the degree of precision of a set of item parameter estimates when the corresponding standard errors of
#' the estimates are not presented in literature or research reports.
#'
#' Sometimes item parameters need to be estimated using different sample size. If the item parameters in the argument \code{x} were
#' calibrated with different number of examinees, a vector of different sample sizes should be specified in the argument \code{nstd}. Suppose
#' that you want to compute the variance-covariance matrices of five IRT 3PLM items and the five items were calibrated with 500, 600, 1,000, 2,000,
#' and 700 examinees, respectively. Then, \code{nstd = c(500, 600, 1000, 2000, 700)} must be specified.
#'
#' Because you can specify only "GPCM" for both the partial credit model (PCM) or the generalized partial credit model (GPCM) in the item metadata,
#' you must indicate which items are the PCM items through the argument \code{pcm.loc}. This is because the item category difficulty parameters are estimated
#' from the PCM, meaning that the variance-covariance of item parameter estimates must be computed for the item category difficulty parameters. Suppose
#' that you want to compute the variance-covariance matrices of five polytomous items and the last two items were calibrated with the PCM. Then,
#' \code{pcm.loc = c(4, 5)} must be specified.
#'
#' @return A list of two internal objects. The first internal object contains a list of the variance-covariance matrices of item parameter estimates.
#' The second internal object contains a list of the standard errors of item parameter estimates.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Li, Y. & Lissitz, R. (2004). Applications of the analytically derived asymptotic standard errors of item response theory
#' item parameter estimates. \emph{Journal of educational measurement, 41}(2), 85-117.
#'
#' Thissen, D. & Wainer, H. (1982). Weighted likelihood estimation of ability in item response theory.
#' \emph{Psychometrika, 54}(3), 427-450.
#'
#' @seealso \code{\link{irtfit}}, \code{\link{info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{gen.weight}}
#'
#' @export
#'
#' @examples
#' ## the use of a "-prm.txt" file obtained sfrom a flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # select the first two dichotomous items and last polytomous item
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # compute the var-covariance matrices with sample size of 2,000
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
