#' Extract Components from 'est_irt', 'est_mg', or 'est_item' Objects
#'
#' @description Extracts internal components from an object of class `est_irt`
#' (from [irtQ::est_irt()]), `est_mg` (from [irtQ::est_mg()]), or `est_item`
#' (from [irtQ::est_item()]).
#'
#' @param x An object of class `est_irt`, `est_mg`, or `est_item` as returned by
#'   [irtQ::est_irt()], [irtQ::est_mg()], or [irtQ::est_item()], respectively.
#' @param what A character string specifying the name of the internal component
#'   to extract.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @details
#' The following components can be extracted from an object of class `est_irt`
#' created by [irtQ::est_irt()]:
#'
#' \describe{
#'   \item{estimates}{A data frame containing both the item parameter estimates
#'   and their corresponding standard errors.}
#'   \item{par.est}{A data frame containing only the item parameter estimates.}
#'   \item{se.est}{A data frame containing the standard errors of the item parameter
#'   estimates, calculated using the cross-product approximation method (Meilijson, 1989).}
#'   \item{pos.par}{A data frame indicating the position index of each estimated
#'   item parameter. This is useful when interpreting the variance-covariance matrix.}
#'   \item{covariance}{A variance-covariance matrix of the item parameter estimates.}
#'   \item{loglikelihood}{The total marginal log-likelihood value summed across all items.}
#'   \item{aic}{Akaike Information Criterion (AIC) based on the marginal log-likelihood.}
#'   \item{bic}{Bayesian Information Criterion (BIC) based on the marginal log-likelihood.}
#'   \item{group.par}{A data frame containing the mean, variance, and standard
#'   deviation of the latent variable's prior distribution.}
#'   \item{weights}{A two-column data frame containing quadrature points
#'   (first column) and corresponding weights (second column) of the (updated)
#'   latent trait prior.}
#'   \item{posterior.dist}{A matrix of normalized posterior densities for all
#'   response patterns at each quadrature point. Rows represent examinees, and
#'   columns represent quadrature points.}
#'   \item{data}{A data frame of the examinee response dataset used in estimation.}
#'   \item{scale.D}{The scaling constant (usually 1 or 1.7) used in the IRT model.}
#'   \item{ncase}{The number of unique response patterns.}
#'   \item{nitem}{The number of items included in the dataset.}
#'   \item{Etol}{The convergence criterion used for the E-step in the EM algorithm.}
#'   \item{MaxE}{The maximum number of E-steps allowed during EM estimation.}
#'   \item{aprior}{A list describing the prior distribution for item slope parameters.}
#'   \item{bprior}{A list describing the prior distribution for item difficulty
#'   (or threshold) parameters.}
#'   \item{gprior}{A list describing the prior distribution for item guessing parameters.}
#'   \item{npar.est}{The total number of parameters estimated.}
#'   \item{niter}{The number of EM cycles completed.}
#'   \item{maxpar.diff}{The maximum change in parameter estimates at convergence.}
#'   \item{EMtime}{Computation time (in seconds) for the EM algorithm.}
#'   \item{SEtime}{Computation time (in seconds) for estimating standard errors.}
#'   \item{TotalTime}{Total computation time (in seconds) for model estimation.}
#'   \item{test.1}{Result of the first-order test indicating whether the gradients
#'   were sufficiently close to zero.}
#'   \item{test.2}{Result of the second-order test indicating whether the
#'   information matrix was positive definite (a condition for maximum likelihood).}
#'   \item{var.note}{A note indicating whether the variance-covariance matrix was
#'   successfully derived from the information matrix.}
#'   \item{fipc}{Logical value indicating whether Fixed Item Parameter Calibration
#'   (FIPC) was applied.}
#'   \item{fipc.method}{The specific method used for FIPC.}
#'   \item{fix.loc}{An integer vector indicating the positions of fixed items used
#'   during FIPC.}
#' }
#'
#'
#'  Components that can be extracted from an object of class `est_mg` created by
#'  [irtQ::est_mg()] include:
#'
#' \describe{
#'   \item{estimates}{A list with two components: `overall` and `group`.
#'   - `overall`: A data frame containing item parameter estimates and their
#'   standard errors, based on the combined data set across all groups.
#'   - `group`: A list of group-specific data frames containing item parameter
#'   estimates and standard errors for each group.}
#'
#'   \item{par.est}{Same structure as `estimates`, but containing only the item
#'   parameter estimates (without standard errors).}
#'   \item{se.est}{Same structure as `estimates`, but containing only the standard
#'   errors of the item parameter estimates. The standard errors are computed
#'   using the cross-product approximation method (Meilijson, 1989).}
#'   \item{pos.par}{A data frame indicating the position index of each estimated
#'   parameter. This index is based on the combined item set across all groups
#'   and is useful when interpreting the variance-covariance matrix.}
#'   \item{covariance}{A variance-covariance matrix for the item parameter
#'   estimates based on the combined data from all groups.}
#'   \item{loglikelihood}{A list with `overall` and `group` components:
#'   - `overall`: The marginal log-likelihood summed over all unique items across all groups.
#'   - `group`: Group-specific marginal log-likelihood values.}
#'
#'   \item{aic}{Akaike Information Criterion (AIC) computed from the overall log-likelihood.}
#'   \item{bic}{Bayesian Information Criterion (BIC) computed from the overall log-likelihood.}
#'   \item{group.par}{A list of group-specific summary statistics (mean, variance,
#'   and standard deviation) of the latent trait prior distribution.}
#'   \item{weights}{A list of two-column data frames (one per group) containing
#'   the quadrature points (first column) and the corresponding weights (second column)
#'   for the updated prior distributions.}
#'   \item{posterior.dist}{A matrix of normalized posterior densities for all
#'   response patterns at each quadrature point. Rows correspond to individuals,
#'   and columns to quadrature points.}
#'   \item{data}{A list with `overall` and `group` components, each containing
#'   examinee response data.}
#'   \item{scale.D}{The scaling constant used in the IRT model (typically 1 or 1.7).}
#'   \item{ncase}{A list with `overall` and `group` components indicating the
#'   number of response patterns in each.}
#'   \item{nitem}{A list with `overall` and `group` components indicating the
#'   number of items in the respective response sets.}
#'   \item{Etol}{Convergence criterion used for the E-step in the EM algorithm.}
#'   \item{MaxE}{Maximum number of E-steps allowed in the EM algorithm.}
#'   \item{aprior}{A list describing the prior distribution for item slope parameters.}
#'   \item{gprior}{A list describing the prior distribution for item guessing parameters.}
#'   \item{npar.est}{Total number of parameters estimated across all unique items.}
#'   \item{niter}{Number of EM cycles completed.}
#'   \item{maxpar.diff}{Maximum change in item parameter estimates at convergence.}
#'   \item{EMtime}{Computation time (in seconds) for EM estimation.}
#'   \item{SEtime}{Computation time (in seconds) for estimating standard errors.}
#'   \item{TotalTime}{Total computation time (in seconds) for model estimation.}
#'   \item{test.1}{First-order condition test result indicating whether gradients
#'   converged sufficiently.}
#'   \item{test.2}{Second-order condition test result indicating whether the
#'   information matrix is positive definite.}
#'   \item{var.note}{A note indicating whether the variance-covariance matrix
#'   was successfully derived from the information matrix.}
#'   \item{fipc}{Logical value indicating whether Fixed Item Parameter Calibration
#'   (FIPC) was used.}
#'   \item{fipc.method}{The method used for FIPC.}
#'   \item{fix.loc}{A list with `overall` and `group` components specifying the
#'   locations of fixed items when FIPC was applied.}
#' }
#'
#'
#' Components that can be extracted from an object of class `est_item` created by
#' [irtQ::est_item()] include:
#'
#' \describe{
#'   \item{estimates}{A data frame containing both the item parameter estimates
#'   and their corresponding standard errors.}
#'   \item{par.est}{A data frame containing only the item parameter estimates.}
#'   \item{se.est}{A data frame containing the standard errors of the item parameter
#'   estimates, computed using observed information functions.}
#'   \item{pos.par}{A data frame indicating the position index of each estimated
#'   item parameter. This is useful when interpreting the variance-covariance matrix.}
#'   \item{covariance}{A variance-covariance matrix of the item parameter estimates.}
#'   \item{loglikelihood}{The sum of log-likelihood values across all items in
#'   the complete data set.}
#'   \item{data}{A data frame of examinee response data.}
#'   \item{score}{A numeric vector of examinees' ability values used as fixed
#'   effects during estimation.}
#'   \item{scale.D}{The scaling constant (typically 1 or 1.7) used in the IRT model.}
#'   \item{convergence}{A character string indicating the convergence status of
#'   the item parameter estimation.}
#'   \item{nitem}{The total number of items included in the response data.}
#'   \item{deleted.item}{Items that contained no response data and were excluded
#'   from estimation.}
#'   \item{npar.est}{The total number of estimated item parameters.}
#'   \item{n.response}{An integer vector indicating the number of responses used
#'   to estimate parameters for each item.}
#'   \item{TotalTime}{Total computation time (in seconds) for the estimation process.}
#' }
#'
#' See [irtQ::est_irt()], [irtQ::est_mg()], and [irtQ::est_item()] for more details.
#'
#' @return
#' The internal component extracted from an object of class `est_irt`, `est_mg`, or `est_item`,
#' depending on the input to the `x` argument.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::est_mg()], [irtQ::est_item()]
#'
#' @examples
#' \donttest{
#' # Fit a 2PL model to the LSAT6 data
#' mod.2pl <- est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2)
#'
#' # Extract item parameter estimates
#' (est.par <- getirt(mod.2pl, what = "par.est"))
#'
#' # Extract standard error estimates
#' (est.se <- getirt(mod.2pl, what = "se.est"))
#'
#' # Extract the variance-covariance matrix of item parameter estimates
#' (cov.mat <- getirt(mod.2pl, what = "covariance"))
#' }
#'
#' @export
getirt <- function(x, ...) UseMethod("getirt")

#' @describeIn getirt An object created by the function [irtQ::est_irt()].
#' @export
getirt.est_irt <- function(x, what, ...) {
  rst <- switch(what,
    estimates = x$estimates,
    par.est = x$par.est,
    se.est = x$se.est,
    pos.par = x$pos.par,
    covariance = x$covariance,
    loglikelihood = x$loglikelihood,
    aic = x$aic,
    bic = x$bic,
    group.par = x$group.par,
    weights = x$weights,
    data = x$data,
    ncase = x$ncase,
    nitem = x$nitem,
    Etol = x$Etol,
    MaxE = x$MaxE,
    aprior = x$aprior,
    bprior = x$bprior,
    gprior = x$gprior,
    npar.est = x$npar.est,
    niter = x$niter,
    maxpar.diff = x$maxpar.diff,
    EMtime = x$EMtime,
    SEtime = x$SEtime,
    TotalTime = x$TotalTime,
    test.1 = x$test.1,
    test.2 = x$test.2,
    var.note = x$var.note,
    fipc = x$fipc,
    fipc.method = x$fipc.method,
    fix.loc = x$fix.loc,
    stop(sprintf("Could not extract element \'%s\'", what), call. = FALSE)
  )

  rst
}

#' @describeIn getirt An object created by the function [irtQ::est_mg()].
#' @export
getirt.est_mg <- function(x, what, ...) {
  rst <- switch(what,
    estimates = x$estimates,
    par.est = x$par.est,
    se.est = x$se.est,
    pos.par = x$pos.par,
    covariance = x$covariance,
    loglikelihood = x$loglikelihood,
    aic = x$aic,
    bic = x$bic,
    group.par = x$group.par,
    weights = x$weights,
    data = x$data,
    ncase = x$ncase,
    nitem = x$nitem,
    Etol = x$Etol,
    MaxE = x$MaxE,
    aprior = x$aprior,
    bprior = x$bprior,
    gprior = x$gprior,
    npar.est = x$npar.est,
    niter = x$niter,
    maxpar.diff = x$maxpar.diff,
    EMtime = x$EMtime,
    SEtime = x$SEtime,
    TotalTime = x$TotalTime,
    test.1 = x$test.1,
    test.2 = x$test.2,
    var.note = x$var.note,
    fipc = x$fipc,
    fipc.method = x$fipc.method,
    fix.loc = x$fix.loc,
    stop(sprintf("Could not extract element \'%s\'", what), call. = FALSE)
  )

  rst
}


#' @describeIn getirt An object created by the function [irtQ::est_item()].
#' @export
getirt.est_item <- function(x, what, ...) {
  rst <- switch(what,
    estimates = x$estimates,
    par.est = x$par.est,
    se.est = x$se.est,
    pos.par = x$pos.par,
    covariance = x$covariance,
    loglikelihood = x$loglikelihood,
    group.par = x$group.par,
    data = x$data,
    score = x$score,
    scale.D = x$scale.D,
    convergence = x$convergence,
    nitem = x$nitem,
    deleted.item = x$deleted.item,
    npar.est = x$npar.est,
    n.response = x$n.response,
    TotalTime = x$TotalTime,
    stop(sprintf("Could not extract element \'%s\'", what), call. = FALSE)
  )

  rst
}
