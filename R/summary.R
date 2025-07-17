#' Summary of Item Calibration Results
#'
#' This S3 method summarizes the IRT calibration results from an object of class `est_irt`,
#' `est_mg`, or `est_item`, which are returned by the functions [irtQ::est_irt()],
#' [irtQ::est_mg()], and [irtQ::est_item()], respectively.
#'
#' @param object An object of class `est_irt`, `est_mg`, or `est_item`.
#' @param ... Additional arguments passed to or from other methods (currently not used).
#'
#' @return A list of internal components extracted from the given object. In addition,
#' the summary method prints an overview of the IRT calibration results to the console.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::est_mg()], [irtQ::est_item()]
#' @examples
#' \donttest{
#' # Fit the 1PL model to LSAT6 data and constrain the slope parameters to be equal
#' fit.1pl <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2, fix.a.1pl = FALSE)
#'
#' # Display the calibration summary
#' summary(fit.1pl)
#' }
#'
#' @export
summary <- function(object, ...) UseMethod("summary")

#' @describeIn summary An object created by the function [irtQ::est_irt()].
#' @export
summary.est_irt <- function(object, ...) {
  call.expr <- deparse(object$call)
  nitem <- object$nitem
  ncase <- object$ncase
  MaxE <- object$MaxE
  Etol <- object$Etol
  weights <- object$weights
  npar.est <- object$npar.est
  fix.loc <- object$fix.loc
  niter <- object$niter
  maxpar.diff <- object$maxpar.diff
  EMtime <- object$EMtime
  SEtime <- object$SEtime
  TotalTime <- object$TotalTime
  test.1 <- object$test.1
  test.2 <- object$test.2
  var.note <- object$var.note
  loglikelihood <- object$loglikelihood
  aic <- object$aic
  bic <- object$bic
  estimates <- object$estimates
  group.par <- object$group.par

  out <- list(
    call.expr = call.expr, nitem = nitem, ncase = ncase, MaxE = MaxE, Etol = Etol,
    weights = weights, npar.est = npar.est, fix.loc = fix.loc, niter = niter,
    maxpar.diff = maxpar.diff, EMtime = EMtime, SEtime = SEtime, TotalTime = TotalTime,
    test.1 = test.1, test.2 = test.2, var.note = var.note, loglikelihood = loglikelihood,
    aic = aic, bic = bic, estimates = estimates, group.par = group.par
  )
  class(out) <- "summary.est_irt"
  out
}

#' @describeIn summary An object created by the function [irtQ::est_mg()].
#' @export
summary.est_mg <- function(object, ...) {
  call.expr <- deparse(object$call)
  nitem <- object$nitem
  ncase <- object$ncase
  MaxE <- object$MaxE
  Etol <- object$Etol
  weights <- object$weights
  npar.est <- object$npar.est
  niter <- object$niter
  maxpar.diff <- object$maxpar.diff
  EMtime <- object$EMtime
  SEtime <- object$SEtime
  TotalTime <- object$TotalTime
  test.1 <- object$test.1
  test.2 <- object$test.2
  var.note <- object$var.note
  loglikelihood <- object$loglikelihood
  aic <- object$aic
  bic <- object$bic
  estimates <- object$estimates
  group.par <- object$group.par
  group.name <- names(object$loglikelihood$group)
  ngroup <- length(group.name)
  fix.loc <- object$fix.loc
  if (is.null(fix.loc$overal)) {
    fix.loc$group <- purrr::map(.x = 1:ngroup, ~ {
      NULL
    })
  }
  names(fix.loc$group) <- group.name

  out <- list(
    call.expr = call.expr, nitem = nitem, ncase = ncase, MaxE = MaxE, Etol = Etol,
    weights = weights, npar.est = npar.est, fix.loc = fix.loc, niter = niter,
    maxpar.diff = maxpar.diff, EMtime = EMtime, SEtime = SEtime, TotalTime = TotalTime,
    test.1 = test.1, test.2 = test.2, var.note = var.note, loglikelihood = loglikelihood,
    aic = aic, bic = bic, estimates = estimates, group.par = group.par, group.name = group.name,
    ngroup = ngroup
  )
  class(out) <- "summary.est_mg"
  out
}


#' @describeIn summary An object created by the function [irtQ::est_item()].
#' @export
summary.est_item <- function(object, ...) {
  call.expr <- deparse(object$call)
  nitem <- object$nitem
  deleted.item <- object$deleted.item
  nitem.del <- length(deleted.item)
  npar.est <- object$npar.est
  n.response <- object$n.response
  TotalTime <- object$TotalTime
  convergence <- object$convergence
  loglikelihood <- object$loglikelihood
  estimates <- object$estimates
  group.par <- object$group.par

  out <- list(
    call.expr = call.expr, nitem = nitem, deleted.item = deleted.item, nitem.del = nitem.del, npar.est = npar.est,
    n.response = n.response, TotalTime = TotalTime, convergence = convergence, loglikelihood = loglikelihood,
    estimates = estimates, group.par = group.par
  )
  class(out) <- "summary.est_item"
  out
}
