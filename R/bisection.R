#' The Bisection Method to Find a Root
#'
#' This function is a modified version of the `bisection` function
#' in the \pkg{cmna} R package (Howard, 2017), designed to find a root of the function
#' `.fun` with respect to its first argument. Unlike the original `bisection()` in
#' \pkg{cmna}, this version allows additional arguments to be passed to `.fun`.
#'
#' @param .fun A function for which the root is to be found.
#' @param ... Additional arguments to be passed to `.fun`.
#' @param lb A numeric value specifying the lower bound of the search interval.
#' @param ub A numeric value specifying the upper bound of the search interval.
#' @param tol A numeric value specifying the tolerance for convergence. Default is 1e-4.
#' @param max.it An integer specifying the maximum number of iterations. Default is 100.
#'
#' @details The bisection method is a well-known root-finding numerical algorithm
#' that applies to any continuous function, provided that the function values at the lower (`lb`)
#' and upper (`ub`) bounds have opposite signs. The method repeatedly bisects the interval until
#' the absolute difference between successive estimates is smaller than the error tolerance
#' (`tol`) or the maximum number of iterations (`max.it`) is reached.
#'
#' @return A list with the following components:
#' - `root`: The estimated root of the function.
#' - `iter`: The number of iterations performed.
#' - `accuracy`: The final absolute difference between the last two interval
#'  points with opposite signs.
#'
#' @seealso [irtQ::est_score()]
#'
#' @references Howard, J. P. (2017). *Computational methods for numerical
#' analysis with R*. New York: Chapman and Hall/CRC.
#'
#' @examples
#' ## Example: Find the theta value corresponding to a given probability
#' ## of a correct response using the item response function of a 2PLM
#' ## (a = 1, b = 0.2)
#'
#' # Define a function of theta
#' find.th <- function(theta, p) {
#'   p - drm(theta = theta, a = 1, b = 0.2, D = 1)
#' }
#'
#' # Find the theta value corresponding to p = 0.2
#' bisection(.fun = find.th, p = 0.2, lb = -10, ub = 10)$root
#'
#' # Find the theta value corresponding to p = 0.8
#' bisection(.fun = find.th, p = 0.8, lb = -10, ub = 10)$root
#'
#' @export
bisection <- function(.fun, ..., lb, ub, tol = 1e-4, max.it = 100) {
  iter <- 0
  f.ub <- .fun(ub, ...)
  while (abs(lb - ub) > tol) {
    mb <- (lb + ub) / 2
    f.mb <- .fun(mb, ...)
    if (f.mb == 0) {
      ub <- mb + (tol / 2)
      lb <- mb - (tol / 2)
    } else if (f.ub * f.mb < 0) {
      lb <- mb
    } else {
      ub <- mb
      f.ub <- f.mb
    }
    iter <- iter + 1
    if (iter > max.it) {
      warning("The maximum number of iteration is reached. \n", call. = FALSE)
      break
    }
  }
  root <- (lb + ub) / 2
  list(root = root, iter = iter, delta = abs(lb - ub))
}
