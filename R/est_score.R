#' Estimate examinees' ability (proficiency) parameters
#'
#' @description This function estimates examinees' latent ability parameters. Available scoring methods are maximum likelihood estimation (ML),
#' maximum likelihood estimation with fences (MLF; Han, 2016), weighted likelihood estimation (Warm, 1989), maximum a posteriori estimation
#' (MAP; Hambleton et al., 1991), expected a posteriori estimation (EAP; Bock & Mislevy, 1982), EAP summed scoring
#' (Thissen et al., 1995; Thissen & Orlando, 2001), and inverse test characteristic curve (TCC) scoring
#' (e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996).
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...) or an object of
#' class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}. See \code{\link{irtfit}}, \code{\link{info}},
#' or \code{\link{simdat}} for more details about the item metadata. This data frame can be easily obtained using the function \code{\link{shape_df}}.
#' @param data A matrix or vector containing examinees' response data for the items in the argument \code{x}. When a matrix is used, a row and column indicate
#' the examinees and items, respectively. When a vector is used, it should contains the item response data for an examinee.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param method A character string indicating a scoring method. Available methods are "ML" for the maximum likelihood estimation,
#' "MLF" for the maximum likelihood estimation with fences, "WL" for the weighted likelihood estimation, "MAP" for the maximum a posteriori estimation,
#' "EAP" for the expected a posteriori estimation, "EAP.SUM" for the expected a posteriori summed scoring, and "INV.TCC" for the inverse TCC scoring.
#' Default method is "ML".
#' @param range A numeric vector of two components to restrict the range of ability scale for the ML, MLF, WL, and MAP scoring methods. Default is c(-5, 5).
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
#' c(0,1). Ignored if \code{method} is "ML", "MLF", "WL", or "INV.TCC".
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
#' Ignored if \code{method} is "ML", "MLF", "WL", "MAP", or "INV.TCC".
#' @param weights A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
#' using the function \code{\link{gen.weight}}. If NULL and \code{method} is "EAP" or "EAP.SUM", default values are used (see the arguments
#' of \code{norm.prior} and \code{nquad}). Ignored if \code{method} is "ML", "MLF", "WL", "MAP", or "INV.TCC".
#' @param fence.a A numeric value specifying the item slope parameter (i.e., \emph{a}-parameter) for the two imaginary items in MLF. See below for details.
#' Default is 3.0.
#' @param fence.b A numeric vector of two components specifying the lower and upper fences of item difficulty parameters (i.e., \emph{b}-parameters)
#' for the two imaginary items, respectively, in MLF. When \code{fence.b = NULL}, the \code{range} values were used to set the lower and upper fences of
#' item difficulty parameters. Default is NULL.
#' @param tol A numeric value of the convergent tolerance for the ML, MLF, WL, MAP, and inverse TCC scoring methods. For the ML, MLF, WL, and MAP,
#' Newton Raphson method is implemented for optimization. For the inverse TCC scoring, the bisection method is used. Default is 1e-4.
#' @param max.iter An positive integer value specifying the maximum number of iterations of Newton Raphson method. Default is 100.
#' @param stval.opt An positive integer value specifying the starting value option for the ML, MLF, WL, and MAP scoring methods.
#' Available options are 1 for the brute-force method, 2 for the observed sum score-based method, and 3 for setting to 0. Default is 1.
#' See below for details.
#' @param se A logical value. If TRUE, the standard errors of ability estimates are computed. However, if \code{method} is "EAP.SUM" or "INV.TCC", the standard
#' errors are always returned. Default is TRUE.
#' @param intpol A logical value. If TRUE and \code{method = "INV.TCC"}, linear interpolation method is used to approximate the ability estimates
#' corresponding to the observed sum scores in which ability estimates cannot be obtained using the TCC (e.g., observed sum scores less than
#' the sum of item guessing parameters). Default is TRUE. See below for details.
#' @param range.tcc A numeric vector of two components to be used as the lower and upper bounds of ability estimates when \code{method = "INV.TCC"}.
#' Default is c(-7, 7).
#' @param missing A value indicating missing values in the response data set. Default is NA. See below for details.
#' @param ncore The number of logical CPU cores to use. Default is 1. See below for details.
#' @param ... additional arguments to pass to \code{parallel::makeCluster}.
#'
#' @details For MAP scoring method, only the normal prior distribution is available for the population distribution.
#'
#' When there are missing data in the response data set, the missing value must be specified in \code{missing}. The missing data are taken into account
#' when either of ML, MLF, WL, MAP, and EAP is used. When "EAP.SUM" or "INV.TCC" is used, however, any missing responses are replaced with incorrect
#' responses (i.e., 0s).
#'
#' In the maximum likelihood estimation with fences (MLF; Han, 2016), two 2PLM imaginary items are necessary. The first imaginary item serves as the lower
#' fence and its difficulty parameter (i.e., \emph{b}-parameters) should be lower than any difficulty parameter values in the test form. Likewise, the second
#' imaginary item serves as the upper fence and its difficulty parameter should be greater than any difficulty parameter values in the test form. Also, the two
#' imaginary items should have a very high item slope parameter (i.e., \emph{a}-parameter) value. See Han (2016) for more details. When \code{fence.b = NULL} in MLF,
#' the function automatically sets the lower and upper fences of item difficulty parameters using the values
#' in the \code{range} argument.
#'
#' When "INV.TCC" method is used employing the IRT 3-parameter logistic model (3PLM) in a test, ability estimates for the observed sum scores less than the
#' sum of item guessing parameters are not attainable. In this case, linear interpolation can be applied by setting \code{intpol = TRUE}.
#' Let \eqn{\theta_{min}} and \eqn{\theta_{max}} be the minimum and maximum ability estimates and \eqn{\theta_{X}} be the ability estimate for
#' the smallest observed sum score, X, but greater than or equal to the sum of item guessing parameters. When linear interpolation method is used,
#' the first value of the \code{range.tcc} is set to \eqn{\theta_{min}}. Then, a linear line is constructed between two points of
#' (x=\eqn{\theta_{min}}, y=0) and (x=\eqn{\theta_{X}}, y=X). Also, the first value of the \code{range.tcc} is set to \eqn{\theta_{max}}, which is
#' the ability estimates corresponding to the maximum observed sum score. When it comes to the scoring method of "INV.TCC", the standard errors of ability
#' estimates are computed using an approach suggested by Lim, Davey, and Wells (2020). The code for the inverse TCC scoring was written by modifying
#' the function \code{irt.eq.tse} of the \pkg{SNSequate} R package (González, 2014).
#'
#' In terms of the starting value to be used for ML, MLF, WL, and MAP scoring methods, the brute-force method is used when \code{stval.opt = 1}. With this option,
#' the log-likelihood values were evaluated at the discrete theta values with increments of 0.1 given \code{range}. The theta node that has the largest
#' log-likelihood is used as the starting value. when \code{stval.opt = 2}, the starting value is obtained based on the observed sum score. For example,
#' if the maximum observed sum score (max.score) is 30 for a test and an examinee has an observed sum score of 20 (obs.score), then the starting value
#' is "log(obs.score / (max.score - obs.score))". For all incorrect response, the starting value is "log(1 / max.score)" and for all correct responses,
#' it is "log(max.score / 1)".
#'
#' To speed up the ability estimation for ML, MLF, WL, MAP, and EAP methods, this function applies a parallel process using multiple logical CPU cores.
#' You can set the number of logical CPU cores by specifying a positive integer value in the argument \code{ncore}. Default value is 1.
#'
#' Note that the standard errors of ability estimates are computed using the Fisher expected information for ML, MLF, WL, and MAP methods.
#'
#' To implement WL method, the \code{Pi}, \code{Ji}, and \code{Ii} functions of \pkg{catR}
#' (Magis & Barrada, 2017) were referred.
#'
#' @return When \code{method} is either of "ML", "MLF", "WL", "MAP", or "EAP", a two column data frame including the ability estimates (1st column)
#' and the standard errors of ability estimates (2nd column) is returned. When \code{method} is "EAP.SUM" or "INV.TCC", a list of two internal
#' objects are returned. The first object is a three column data frame including the observed sum scores (1st column), the ability estimates (2nd column),
#' and the standard errors of ability estimates (3rd column). The second object is a score table including the possible raw sum scores
#' and corresponding ability and standard error estimates.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{irtfit}}, \code{\link{info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{gen.weight}}
#'
#' @references
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Psychometrika, 35}, 179-198.
#'
#' González, J. (2014). SNSequate: Standard and nonstandard statistical models and methods for test equating.
#' \emph{Journal of Statistical Software, 59}, 1-30.
#'
#' Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of item response theory}. Newbury Park, CA: Sage.
#'
#' Han, K. T. (2016). Maximum likelihood score estimation method with fences for short-length tests and computerized adaptive tests.
#' \emph{Applied psychological measurement, 40}(4), 289-301.
#'
#' Howard, J. P. (2017). \emph{Computational methods for numerical analysis with R}. New York:
#' Chapman and Hall/CRC.
#'
#' Kolen, M. J. & Brennan, R. L. (2004). \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
#' Springer
#'
#' Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT proficiency estimates.
#' \emph{Educational Measurement: Issues and Practice, 29}(3), 8-14.
#'
#' Lim, H., Davey, T., & Wells, C. S. (2020). A recursion-based analytical approach to evaluate the performance of MST.
#' \emph{Journal of Educational Measurement, 58}(2), 154-178.
#'
#' Magis, D., & Barrada, J. R. (2017). Computerized adaptive testing with R: Recent updates of the package catR.
#' \emph{Journal of Statistical Software, 76}, 1-19.
#'
#' Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
#' \emph{Journal of Educational and Behavioral Statistics, 21}(4), 365-389.
#'
#' Thissen, D. & Orlando, M. (2001). Item response theory for items scored in two categories. In D. Thissen & H. Wainer (Eds.),
#' \emph{Test scoring} (pp.73-140). Mahwah, NJ: Lawrence Erlbaum.
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item Response Theory
#' for Scores on Tests Including Polytomous Items with Ordered Responses. \emph{Applied Psychological
#' Measurement, 19}(1), 39-49.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory. \emph{Psychometrika, 54}(3),
#' 427-450.
#'
#' @examples
#' ## the use of a "-prm.txt" file obtained from a flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # read item parameters and transform them to item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # generate examinees abilities
#' set.seed(12)
#' theta <- rnorm(10)
#'
#' # simulate the item response data
#' data <- simdat(x, theta, D = 1)
#'
#' \donttest{
#' # estimate the abilities using ML
#' est_score(x, data, D = 1, method = "ML", range = c(-4, 4), se = TRUE)
#'
#' # estimate the abilities using WL
#' est_score(x, data, D = 1, method = "WL", range = c(-4, 4), se = TRUE)
#'
#' # estimate the abilities using MLF with default fences of item difficulty parameters
#' est_score(x, data, D = 1, method = "MLF", fence.a = 3.0, fence.b = NULL, se = TRUE)
#'
#' # estimate the abilities using MLF with different fences of item difficulty parameters
#' est_score(x, data, D = 1, method = "MLF", fence.a = 3.0, fence.b = c(-7, 7), se = TRUE)
#'
#' # estimate the abilities using MAP
#' est_score(x, data, D = 1, method = "MAP", norm.prior = c(0, 1), nquad = 30, se = TRUE)
#'
#' # estimate the abilities using EAP
#' est_score(x, data, D = 1, method = "EAP", norm.prior = c(0, 1), nquad = 30, se = TRUE)
#'
#' # estimate the abilities using EAP summed scoring
#' est_score(x, data, D = 1, method = "EAP.SUM", norm.prior = c(0, 1), nquad = 30)
#'
#' # estimate the abilities using inverse TCC scoring
#' est_score(x, data, D = 1, method = "INV.TCC", intpol = TRUE, range.tcc = c(-7, 7))
#' }
#'
#' @export
est_score <- function(x, ...) UseMethod("est_score")

#' @describeIn est_score Default method to estimate examinees' latent ability parameters using a data frame \code{x} containing the item metadata.
#' @importFrom reshape2 melt
#' @import dplyr
#' @export
est_score.default <- function(x, data, D = 1, method = "ML", range = c(-5, 5), norm.prior = c(0, 1),
                              nquad = 41, weights = NULL, fence.a = 3.0, fence.b = NULL,
                              tol = 1e-4, max.iter = 100, se = TRUE, stval.opt = 1,
                              intpol = TRUE, range.tcc = c(-7, 7),
                              missing = NA, ncore = 1, ...) {
  # check if the data set is a vector of an examinee
  if (is.vector(data)) {
    data <- rbind(data)
  }

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of examinees
  nstd <- nrow(data)

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # scoring of ML, WL, MLF, MAP, and EAP
  if (method %in% c("ML", "MAP", "WL", "EAP", "MLF")) {
    # check if the method is WL
    # if TRUE, ji = TRUE
    ji <- ifelse(method == "WL", TRUE, FALSE)

    # add two more items and data responses when MLF is used
    if (method == "MLF") {
      # when fence.b = NULL, use the range argument as the fence.b argument
      if (is.null(fence.b)) {
        fence.b <- range
      }

      # add two more response columns for the two fence items
      data <- cbind(data, f.lower = 1, f.upper = 0)

      # create a new item metadata for the two fence items
      x.fence <- shape_df(
        par.drm = list(a = rep(fence.a, 2), b = fence.b, g = rep(0, 2)),
        item.id = c("fence.lower", "fence.upper"), cats = 2,
        model = "3PLM"
      )

      # create the new item metadata by adding two fence items
      x <- dplyr::bind_rows(x, x.fence)
    }

    # check the maximum score category across all items
    max.cats <- max(x$cats)

    # count the number of columns of the item metadata
    max.col <- ncol(x)

    # create an empty elm_item object
    elm_item <- list(pars = NULL, model = NULL, cats = NULL)

    # check the number of CPU cores
    if (ncore < 1) {
      stop("The number of logical CPU cores must not be less than 1.", call. = FALSE)
    }

    # estimation
    if (ncore == 1L) {
      # reshape the data
      if (stval.opt == 2) {
        data <-
          data.frame(x, t(data), check.names = FALSE) %>%
          reshape2::melt(
            id.vars = 1:max.col,
            variable.name = "std",
            value.name = "resp.num",
            na.rm = TRUE,
            factorsAsStrings = FALSE
          )
        data$resp <- factor(data$resp.num, levels = (seq_len(max.cats) - 1))
      } else {
        data <-
          data.frame(x, t(data), check.names = FALSE) %>%
          reshape2::melt(
            id.vars = 1:max.col,
            variable.name = "std",
            value.name = "resp",
            na.rm = TRUE,
            factorsAsStrings = FALSE
          )
        data$resp <- factor(data$resp, levels = (seq_len(max.cats) - 1))
      }
      data$std <- as.numeric(data$std)

      # create a score table to contain the scoring results
      rst <- data.frame(std = 1:nstd)

      # scoring
      est <-
        dplyr::group_by(data, .data$std) %>%
        dplyr::summarise(
          est_score_indiv(
            exam_dat = dplyr::pick(dplyr::everything()), elm_item = elm_item,
            max.col = max.col, D = D, method = method,
            range = range, norm.prior = norm.prior, nquad = nquad,
            weights = weights, tol = tol, max.iter = max.iter, se = se,
            stval.opt = stval.opt, ji = ji
          )
        )

      # merge the scoring results
      rst <- merge(x = rst, y = est, by = "std")
      rst <- rst[, -1]
    } else {
      # create a parallel processing cluster
      # cl <- parallel::makeCluster(ncore)
      cl <- parallel::makeCluster(ncore, ...)

      # divide data into
      quotient <- nstd %/% ncore
      remain <- nstd %% ncore
      data_list <- vector("list", ncore)
      for (k in 1:ncore) {
        if (k == ncore & remain != 0) {
          data_list[[k]] <- data[((k - 1) * quotient + 1):(quotient * k + remain), ]
        } else {
          data_list[[k]] <- data[((k - 1) * quotient + 1):(quotient * k), ]
        }
      }

      # delete 'data' object
      rm(data, envir = environment(), inherits = FALSE)

      # load some specific variable names into processing cluster
      parallel::clusterExport(cl, c(
        "x", "elm_item", "D", "method",
        "max.cats", "max.col", "range", "norm.prior", "nquad",
        "weights", "tol", "max.iter", "se", "stval.opt", "ji",
        "est_score_1core", "est_score_indiv", "idxfinder",
        "ll_score", "drm", "prm", "gpcm", "grm",
        "logprior_deriv", "esprior_norm",
        "info_score", "info_drm", "info_prm",
        "gen.weight"
      ), envir = environment())
      parallel::clusterEvalQ(cl, library(dplyr))
      # parallel::clusterEvalQ(cl, library(reshape2))
      # parallel::clusterEvalQ(cl, library(Rfast))


      # set a function for scoring
      fsm <- function(subdat) {
        est_score_1core(
          x = x, elm_item = elm_item, data = subdat, D = D,
          method = method, max.cats = max.cats, max.col = max.col,
          range = range, norm.prior = norm.prior, nquad = nquad,
          weights = weights, tol = tol, max.iter = max.iter,
          se = se, stval.opt = stval.opt, ji = ji
        )
      }

      # parallel scoring
      est <- parallel::parLapply(cl = cl, X = data_list, fun = fsm)

      # finish
      parallel::stopCluster(cl)

      # combine the results
      rst <- do.call(what = "rbind", args = est)
    }

    # return a warning message when some examinees have all missing responses
    loc.na <- which(is.na(rst$est.theta))
    if (length(loc.na) > 0) {
      memo <- paste("NA values were reterned for examinees with all missing responses.")
      warning(memo, call. = FALSE)
    }
  }

  if (method == "EAP.SUM") {
    rst <- eap_sum(
      x = x, data = data, norm.prior = norm.prior,
      nquad = nquad, weights = weights, D = D
    )
  }

  if (method == "INV.TCC") {
    rst <- inv_tcc(x, data,
      D = D, intpol = intpol, range.tcc = range.tcc,
      tol = tol, max.it = 500
    )
  }

  # return results
  rst
}


#' @describeIn est_score An object created by the function \code{\link{est_irt}}.
#' @import dplyr
#' @export
est_score.est_irt <- function(x, method = "ML", range = c(-5, 5), norm.prior = c(0, 1),
                              nquad = 41, weights = NULL, fence.a = 3.0, fence.b = NULL,
                              tol = 1e-4, max.iter = 100, se = TRUE, stval.opt = 1,
                              intpol = TRUE, range.tcc = c(-7, 7),
                              missing = NA, ncore = 1, ...) {
  # extract information from an object
  data <- x$data
  D <- x$scale.D
  x <- x$par.est

  # check if the data set is a vector of an examinee
  if (is.vector(data)) {
    data <- rbind(data)
  }

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of examinees
  nstd <- nrow(data)

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # scoring of ML, WL, MLF, MAP, and EAP
  if (method %in% c("ML", "MAP", "WL", "EAP", "MLF")) {
    # check if the method is WL
    # if TRUE, ji = TRUE
    ji <- ifelse(method == "WL", TRUE, FALSE)

    # add two more items and data responses when "ML" with Fences method is used
    if (method == "MLF") {
      # when fence.b = NULL, use the range argument as the fence.b argument
      if (is.null(fence.b)) {
        fence.b <- range
      }

      # add two more response columns for the two fence items
      data <- cbind(data, f.lower = 1, f.upper = 0)

      # create a new item metadata for the two fence items
      x.fence <- shape_df(
        par.drm = list(a = rep(fence.a, 2), b = fence.b, g = rep(0, 2)),
        item.id = c("fence.lower", "fence.upper"), cats = 2,
        model = "3PLM"
      )

      # create the new item metadata by adding two fence items
      x <- dplyr::bind_rows(x, x.fence)
    }

    # check the maximum score category across all items
    max.cats <- max(x$cats)

    # count the number of columns of the item metadata
    max.col <- ncol(x)

    # create an empty elm_item object
    elm_item <- list(pars = NULL, model = NULL, cats = NULL)

    # check the number of CPU cores
    if (ncore < 1) {
      stop("The number of logical CPU cores must not be less than 1.", call. = FALSE)
    }

    # estimation
    if (ncore == 1L) {
      # reshape the data
      if (stval.opt == 2) {
        data <-
          data.frame(x, t(data), check.names = FALSE) %>%
          reshape2::melt(
            id.vars = 1:max.col,
            variable.name = "std",
            value.name = "resp.num",
            na.rm = TRUE,
            factorsAsStrings = FALSE
          )
        data$resp <- factor(data$resp.num, levels = (seq_len(max.cats) - 1))
      } else {
        data <-
          data.frame(x, t(data), check.names = FALSE) %>%
          reshape2::melt(
            id.vars = 1:max.col,
            variable.name = "std",
            value.name = "resp",
            na.rm = TRUE,
            factorsAsStrings = FALSE
          )
        data$resp <- factor(data$resp, levels = (seq_len(max.cats) - 1))
      }
      data$std <- as.numeric(data$std)

      # create a score table to contain the scoring results
      rst <- data.frame(std = 1:nstd)

      # scoring
      est <-
        dplyr::group_by(data, .data$std) %>%
        dplyr::summarise(
          est_score_indiv(
            exam_dat = dplyr::pick(dplyr::everything()), elm_item = elm_item,
            max.col = max.col, D = D, method = method,
            range = range, norm.prior = norm.prior, nquad = nquad,
            weights = weights, tol = tol, max.iter = max.iter, se = se,
            stval.opt = stval.opt, ji = ji
          )
        )

      # merge the scoring results
      rst <- merge(x = rst, y = est, by = "std")
      rst <- rst[, -1]
    } else {
      # create a parallel processing cluster
      cl <- parallel::makeCluster(ncore, ...)

      # divide data into
      quotient <- nstd %/% ncore
      remain <- nstd %% ncore
      data_list <- vector("list", ncore)
      for (k in 1:ncore) {
        if (k == ncore & remain != 0) {
          data_list[[k]] <- data[((k - 1) * quotient + 1):(quotient * k + remain), ]
        } else {
          data_list[[k]] <- data[((k - 1) * quotient + 1):(quotient * k), ]
        }
      }

      # delete 'data' object
      rm(data, envir = environment(), inherits = FALSE)

      # load some specific variable names into processing cluster
      parallel::clusterExport(cl, c(
        "x", "elm_item", "D", "method",
        "max.cats", "max.col", "range", "norm.prior", "nquad",
        "weights", "tol", "max.iter", "se", "stval.opt", "ji",
        "est_score_1core", "est_score_indiv", "idxfinder",
        "ll_score", "drm", "prm", "gpcm", "grm",
        "logprior_deriv", "esprior_norm",
        "info_score", "info_drm", "info_prm",
        "gen.weight"
      ), envir = environment())
      parallel::clusterEvalQ(cl, library(dplyr))
      # parallel::clusterEvalQ(cl, library(reshape2))
      # parallel::clusterEvalQ(cl, library(Rfast))


      # set a function for scoring
      fsm <- function(subdat) {
        est_score_1core(
          x = x, elm_item = elm_item, data = subdat, D = D,
          method = method, max.cats = max.cats, max.col = max.col,
          range = range, norm.prior = norm.prior, nquad = nquad,
          weights = weights, tol = tol, max.iter = max.iter,
          se = se, stval.opt = stval.opt, ji = ji
        )
      }

      # parallel scoring
      est <- parallel::parLapply(cl = cl, X = data_list, fun = fsm)

      # finish
      parallel::stopCluster(cl)

      # combine the results
      rst <- do.call(what = "rbind", args = est)
    }

    # return a warning message when some examinees have all missing responses
    loc.na <- which(is.na(rst$est.theta))
    if (length(loc.na) > 0) {
      memo <- paste("NA values were reterned for examinees with all missing responses.")
      warning(memo, call. = FALSE)
    }
  }

  if (method == "EAP.SUM") {
    rst <- eap_sum(
      x = x, data = data, norm.prior = norm.prior,
      nquad = nquad, weights = weights, D = D
    )
  }

  if (method == "INV.TCC") {
    rst <- inv_tcc(x, data,
      D = D, intpol = intpol, range.tcc = range.tcc,
      tol = tol, max.it = 500
    )
  }

  # return results
  rst
}


# This function is used for each single core computation
#' @import dplyr
est_score_1core <- function(x, elm_item, data, D = 1, method = "ML",
                            max.cats, max.col, range = c(-4, 4), norm.prior = c(0, 1), nquad = 41,
                            weights = NULL, tol = 1e-4, max.iter = 30, se = TRUE,
                            stval.opt = 1, ji = FALSE) {
  # check the number of examinees
  nstd <- nrow(data)

  # reshape the data
  if (stval.opt == 2) {
    data <-
      data.frame(x, t(data), check.names = FALSE) %>%
      reshape2::melt(
        id.vars = 1:max.col,
        variable.name = "std",
        value.name = "resp.num",
        na.rm = TRUE,
        factorsAsStrings = FALSE
      )
    data$resp <- factor(data$resp.num, levels = (seq_len(max.cats) - 1))
  } else {
    data <-
      data.frame(x, t(data), check.names = FALSE) %>%
      reshape2::melt(
        id.vars = 1:max.col,
        variable.name = "std",
        value.name = "resp",
        na.rm = TRUE,
        factorsAsStrings = FALSE
      )
    data$resp <- factor(data$resp, levels = (seq_len(max.cats) - 1))
  }
  data$std <- as.numeric(data$std)

  # create a score table to contain the scoring results
  rst <- data.frame(std = 1:nstd)

  # scoring
  est <-
    dplyr::group_by(data, .data$std) %>%
    dplyr::summarise(
      est_score_indiv(
        exam_dat = dplyr::pick(dplyr::everything()), elm_item = elm_item,
        max.col = max.col, D = D, method = method,
        range = range, norm.prior = norm.prior, nquad = nquad,
        weights = weights, tol = tol, max.iter = max.iter, se = se,
        stval.opt = stval.opt, ji = ji
      )
    )

  # merge the scoring results
  rst <- merge(x = rst, y = est, by = "std")
  rst <- rst[, -1]

  # return results
  rst
}

# This function computes an ability estimate for an examinee (ML, MLF, MAP, EAP)
#' @importFrom stats xtabs na.pass
est_score_indiv <- function(exam_dat, elm_item, max.col, D = 1, method = "ML",
                            range = c(-4, 4), norm.prior = c(0, 1), nquad = 41, weights = NULL,
                            tol = 1e-4, max.iter = 30, se = TRUE, stval.opt = 1, ji = FALSE) {
  # extract the required objects from the individual exam data
  elm_item$pars <- data.matrix(exam_dat[, 4:max.col])
  elm_item$model <- exam_dat$model
  elm_item$cats <- exam_dat$cats
  resp <- exam_dat$resp
  n.resp <- length(resp)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # calculate the score categories
  tmp.id <- 1:n.resp
  freq.cat <-
    matrix(
      stats::xtabs(~ tmp.id + resp,
        na.action = stats::na.pass, addNA = FALSE
      ),
      nrow = n.resp
    )

  ## ----------------------------------------------------
  ## ML, MLF and MAP
  if (method %in% c("ML", "WL", "MLF", "MAP")) {
    # set a starting value
    if (stval.opt == 1) {
      # use a brute force + smart starting value method
      # prepare the discrete theta values
      theta.nodes <- seq(from = range[1], to = range[2], by = 0.1)

      # compute the negative log-likelihood values for all the discrete theta values
      ll_tmp <- ll_score(
        theta = theta.nodes, elm_item = elm_item, freq.cat = freq.cat,
        method = method, idx.drm = idx.drm, idx.prm = idx.prm, D = D,
        norm.prior = norm.prior, logL = TRUE
      )

      # find the locations of thetas where the sign of slope changes
      # in the negative loglikelihood function
      loc_change <- which(diff(sign(diff(ll_tmp))) > 0L) + 1

      # select a theta value that has the minimum of negative log-likelihood value
      stval_tmp1 <- theta.nodes[loc_change][which.min(ll_tmp[loc_change])]

      # if there is no selected starting value (this means that the negative ll function is
      # monotonically increasing or descreasing), use the 0
      theta <- ifelse(length(stval_tmp1) > 0L, stval_tmp1, 0)
    } else if (stval.opt == 2) {
      # compute a perfect NC score
      total.nc <- sum(elm_item$cats - 1)

      # compute a starting value based on the observed sum score
      obs.sum <- sum(exam_dat$resp.num)
      if (obs.sum == 0) {
        theta <- log(1 / total.nc)
      } else if (obs.sum == total.nc) {
        theta <- log(total.nc / 1)
      } else {
        theta <- log(obs.sum / (total.nc - obs.sum))
      }
    } else if (stval.opt == 3) {
      # use 0 as a starting value
      theta <- 0
    }

    # estimate an ability using Newton-Raphson
    # set the iteration number to 0
    i <- 0
    abs_delta <- 1
    while (abs_delta >= tol) {
      # update the iteration number
      i <- i + 1

      # compute the gradient (gradient of negative log-likelihood)
      # and the fisher information (negative expectation of second derivative of log-likelihood)
      gr_fi <-
        info_score(
          theta = theta, elm_item = elm_item, freq.cat = freq.cat,
          idx.drm = idx.drm, idx.prm = idx.prm, method = method, D = D,
          norm.prior = norm.prior, grad = TRUE, ji = ji
        )
      grad <- gr_fi$grad
      finfo <- gr_fi$finfo

      # protect the fisher information having value close to 0
      finfo[finfo < 1e-5] <- 1e-5

      # compute the theta correction factor (delta)
      delta <- grad / finfo

      # if abs(delta) > 1, assign assign 1 or -1 to the delta
      abs_delta <- abs(delta)
      delta[abs_delta > 1] <- sign(delta)

      # update the theta value
      theta <- theta - delta

      if (i == max.iter) break
    }

    # assign boundary score when the theta estimate is beyond the boundary
    theta[theta <= range[1]] <- range[1]
    theta[theta >= range[2]] <- range[2]
    est.theta <- theta

    # compute the standard error
    if (se) {
      if (est.theta %in% range) {
        se.theta <- 99.9999
      } else {
        finfo <-
          info_score(
            theta = est.theta, elm_item = elm_item, freq.cat = freq.cat,
            idx.drm = idx.drm, idx.prm = idx.prm, method = method, D = D,
            norm.prior = norm.prior, grad = FALSE, ji = FALSE
          )$finfo
        se.theta <- 1 / sqrt(finfo)
      }
    } else {
      se.theta <- NA
    }
  }

  ## ----------------------------------------------------
  ## EAP scoring
  if (method == "EAP") {
    # generate quadrature points and weights
    if (is.null(weights)) {
      popdist <- gen.weight(n = nquad, dist = "norm", mu = norm.prior[1], sigma = norm.prior[2])
    } else {
      popdist <- data.frame(weights)
    }

    # compute the posterior distribution
    posterior <-
      ll_score(
        theta = popdist[, 1], elm_item = elm_item, freq.cat = freq.cat,
        idx.drm = idx.drm, idx.prm = idx.prm, D = D, logL = FALSE
      ) * popdist[, 2]

    # Expected A Posterior
    posterior <- posterior / sum(posterior)
    est.theta <- sum(popdist[, 1] * posterior)

    if (se) {
      # calculating standard error
      ex2 <- sum(popdist[, 1]^2 * posterior)
      var <- ex2 - (est.theta)^2
      se.theta <- sqrt(var)
    } else {
      se.theta <- NA
    }
  }

  # combine theta and se into a list
  rst <- data.frame(est.theta = est.theta, se.theta = se.theta)

  # return results
  rst
}
