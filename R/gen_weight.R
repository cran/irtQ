#' Generate Weights
#'
#' This function generates a set of normalized weights based on theta (ability)
#' values to be used in functions such as [irtQ::est_score()],
#' [irtQ::sx2_fit()], and [irtQ::covirt()].
#'
#' @param n An integer specifying the number of theta (node) values for which
#'   weights are to be generated. Default is 41.
#' @param dist A character string indicating the distribution type used to
#'   generate weights. Available options are `"norm"` for a normal distribution,
#'   `"unif"` for a uniform distribution, and `"emp"` for an empirical
#'   distribution.
#'   - If `dist = "norm"`, either `n` or `theta` must be provided.
#'   - If `dist = "unif"`, only `n` is applicable.
#'   - If `dist = "emp"`, only `theta` must be specified.
#' @param mu,sigma Mean and standard deviation of the normal distribution (used
#'   when `dist = "norm"`).
#' @param l,u Lower and upper bounds of the uniform distribution (used when
#'   `dist = "unif"`).
#' @param theta A numeric vector of empirical theta (node) values for which
#'   weights are generated.
#'
#' @details If `theta` is not specified, *n* equally spaced quadrature points
#' and corresponding weights are generated from either the normal or uniform
#' distribution:
#' - When `dist = "norm"`, Gaussian quadrature points and weights are computed
#' using `gauss.quad.prob()` from the \pkg{statmod} package.
#' - When `dist = "unif"`, equally spaced points are drawn from the specified
#' interval \[`l`, `u`\], and weights are proportional to the uniform density.
#'
#' If `theta` is specified:
#' - When `dist = "norm"`, the weights are proportional to the normal density
#' evaluated at each theta value and normalized to sum to 1.
#' - When `dist = "emp"`, equal weights are assigned to each provided theta value.
#'
#' @return A data frame with two columns:
#' - `theta`: The theta (node) values.
#' - `weight`: The corresponding normalized weights.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_score()], [irtQ::sx2_fit()], [irtQ::covirt()]
#'
#' @examples
#' ## Example 1:
#' ## Generate 41 Gaussian quadrature points and weights from the normal distribution
#' gen.weight(n = 41, dist = "norm", mu = 0, sigma = 1)
#'
#' ## Example 2:
#' ## Generate 41 theta values and weights from the uniform distribution,
#' ## given a minimum value of -4 and a maximum value of 4
#' gen.weight(n = 41, dist = "unif", l = -4, u = 4)
#'
#' ## Example 3:
#' ## Generate normalized weights from the standard normal distribution,
#' ## given a user-defined set of theta values
#' theta <- seq(-4, 4, by = 0.1)
#' gen.weight(dist = "norm", mu = 0, sigma = 1, theta = theta)
#'
#' ## Example 4:
#' ## Generate equal normalized weights for theta values
#' ## randomly sampled from the standard normal distribution
#' theta <- rnorm(100)
#' gen.weight(dist = "emp", theta = theta)
#'
#' @export
#' @importFrom statmod gauss.quad.prob
gen.weight <- function(n = 41,
                       dist = "norm",
                       mu = 0,
                       sigma = 1,
                       l = -4,
                       u = 4,
                       theta) {
  dist <- tolower(dist)

  if (missing(theta)) {
    if (dist == "emp") {
      stop("To use actual option in a distribution argument, theta values are needed.", call. = FALSE)
    }
    if (dist == "norm") {
      wts.nd <- statmod::gauss.quad.prob(n, dist = dist, mu = mu, sigma = sigma)
      nodes <- wts.nd$nodes
      whts <- wts.nd$weights
    }
    if (dist == "unif") {
      nodes <- seq(l, u, length.out = n)
      dens <- stats::dunif(nodes, min = l, max = u)
      whts <- dens / sum(dens)
    }
  }

  if (!missing(theta)) {
    if (dist == "emp") {
      nodes <- theta
      whts <- rep(1 / length(theta), length(theta))
    }
    if (dist == "norm") {
      nodes <- theta
      dens <- stats::dnorm(theta, mean = mu, sd = sigma)
      whts <- dens / sum(dens)
    }
  }

  return(data.frame(theta = nodes, weight = whts))
}
