#' The bisection method to find a root
#'
#' @description This function is a modified version of the \code{bisetion} function
#' in the \pkg{cmna} R package  (Howard, 2017) to find a root of the function \code{.funs}
#' with respect to its first argument. Unlike the \code{bisetion} of the \pkg{cmna},
#' this \code{bisetion} function accepts additional arguments of the function \code{.fun}.
#'
#' @param .fun A function for which the root is searched.
#' @param ... Additional arguments to be passed to \code{.fun}.
#' @param lb A lower bound of the interval to be searched.
#' @param ub An upper bound of the interval to be searched.
#' @param tol The tolerance of error. Default is 1e-4.
#' @param max.it The maximum number of iterations. Default is 100.
#'
#' @details
#' The bisection method is a well-known root finding numerical algorithm that works for any continuous
#' function when the lower (\code{lb}) and upper (\code{ub}) bounds with opposite signs are provided.
#' This method repeatedly bisects the defined interval by two values with opposite signs until the absolute
#' difference of two values becomes less than the error tolerance (\code{tol}) or the maximum
#' number of iterations (\code{max.it}) is reached.
#'
#' @return A list with three internal objects. The first object is the root found, the second object is
#' the number of iterations used, and the third object is the approximate accuracy of the root
#' (i.e., absolute difference between the final two values with opposite signs).
#'
#' @seealso \code{\link{est_score}}
#'
#' @references
#' Howard, J. P. (2017). \emph{Computational methods for numerical analysis with R}. New York:
#' Chapman and Hall/CRC.
#'
#' @examples
#' ## example: find a theta corresponding to the probability of
#' ## correct answer using the item response function of 2PLM
#' ## (a = 1, b = 0.2)
#'
#' # set a function of theta
#' find.th <- function(theta, p) {
#'   p - drm(theta = theta, a = 1, b = 0.2, D = 1)
#' }
#'
#' # find the theta corresponding to p = 0.2
#' bisection(.fun = find.th, p = 0.2, lb = -10, ub = 10)$root
#'
#' # find the theta corresponding to p = 0.8
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
