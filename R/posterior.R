# This function computes the posterior distribution of latent abilities for
# each examinee, given the likelihood of their observed responses at each
# quadrature point and the prior population density.
#
# Args:
#   likehd  : numeric matrix (nstd x ntheta) of likelihoods L(theta_q | resp_i)
#   weights : data.frame with columns (theta, weight) for the single-group
#             case; a list of such data.frames (one per group) for the
#             multi-group case
#   idx.std : NULL for single-group calibration; a list where idx.std[[g]]
#             contains the row indices of `likehd` belonging to group g
#             (used by est_mg() multi-group calibration)
#
# Returns:
#   numeric matrix (nstd x ntheta) of posterior densities P(theta_q | resp_i);
#   each row sums to 1.
#' @importFrom Rfast eachrow
posterior <- function(likehd, weights, idx.std = NULL) {

  ## --------------------------------------------------------------------
  ## (1) Single-group case
  ## --------------------------------------------------------------------
  if (is.null(idx.std)) {

    # extract the prior weights at each quadrature point (length = ntheta)
    w <- weights[, 2]

    # denominator of Bayes' rule per examinee:
    #   denom[i] = sum_q likehd[i, q] * w[q]
    # computed as a single BLAS GEMV (likehd %*% w), giving an nstd-vector;
    # avoids the sweep + rowsums two-pass pattern used previously
    denom <- as.vector(likehd %*% w)

    # numerator: column-scale likehd by w (Rfast::eachrow gives the
    # element-wise product likehd[i, q] * w[q] in C++); then row-divide by
    # denom (length-nstd vector recycles down columns), producing the
    # normalized posterior:
    #   posterior[i, q] = likehd[i, q] * w[q] / denom[i]
    Rfast::eachrow(likehd, w, oper = "*") / denom

  ## --------------------------------------------------------------------
  ## (2) Multi-group case (used by est_mg())
  ## --------------------------------------------------------------------
  } else {

    # number of quadrature points (same across groups by construction)
    ntheta <- nrow(weights[[1]])

    # allocate the full nstd x ntheta output matrix; per-group results are
    # written into the rows belonging to each group below
    out <- array(0, c(nrow(likehd), ntheta))

    # process each group independently: groups can have different prior
    # population densities (different `weights[[g]]`), so the GEMV
    # denominator must be computed per group
    for (g in seq_along(weights)) {

      # row indices of examinees belonging to group g
      idx <- idx.std[[g]]

      # likelihood submatrix for group g; drop = FALSE preserves the matrix
      # structure when group g has only one examinee
      lh_g <- likehd[idx, , drop = FALSE]

      # prior weights for group g at each quadrature point
      w_g <- weights[[g]][, 2]

      # group-g denominator via GEMV: denom_g[i] = sum_q lh_g[i, q] * w_g[q]
      denom_g <- as.vector(lh_g %*% w_g)

      # normalized posterior for group g (same column-scale + row-divide
      # construction as the single-group branch); written back into the
      # rows of `out` that belong to this group
      out[idx, ] <- Rfast::eachrow(lh_g, w_g, oper = "*") / denom_g
    }

    # return the assembled posterior matrix (rows aligned with the original
    # examinee ordering in likehd)
    out
  }
}
