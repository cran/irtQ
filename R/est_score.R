#' Estimate examinees' ability (proficiency) parameters
#'
#' This function estimates examinees' latent ability parameters. Available
#' scoring methods include maximum likelihood estimation (ML), maximum
#' likelihood estimation with fences (MLF; Han, 2016), weighted likelihood
#' estimation (WL; Warm, 1989), maximum a posteriori estimation (MAP; Hambleton
#' et al., 1991), expected a posteriori estimation (EAP; Bock & Mislevy, 1982),
#' EAP summed scoring (Thissen et al., 1995; Thissen & Orlando, 2001), and
#' inverse test characteristic curve (TCC) scoring (e.g., Kolen & Brennan, 2004;
#' Kolen & Tong, 2010; Stocking, 1996).
#'
#' @inheritParams est_irt
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.); or an object of class `est_irt`
#'   obtained from [irtQ::est_irt()], or `est_item` from [irtQ::est_item()].
#'
#'   See [irtQ::est_irt()] or [irtQ::simdat()] for more details about the item
#'   metadata. This data frame can be easily created using the
#'   [irtQ::shape_df()] function.
#' @param method A character string indicating the scoring method to use.
#'   Available options are:
#'   - `"ML"`: Maximum likelihood estimation
#'   - `"MLF"`: Maximum likelihood estimation with fences (Han, 2016)
#'   - `"WL"`: Weighted likelihood estimation (Warm, 1989)
#'   - `"MAP"`: Maximum a posteriori estimation (Hambleton et al., 1991)
#'   - `"EAP"`: Expected a posteriori estimation (Bock & Mislevy, 1982)
#'   - `"EAP.SUM"`: Expected a posteriori summed scoring (Thissen et al., 1995;
#'   Thissen & Orlando, 2001)
#'   - `"INV.TCC"`: Inverse test characteristic curve scoring
#'   (e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996)
#'
#'   Default is `"ML"`.
#' @param range A numeric vector of length two specifying the lower and upper
#'   bounds of the ability scale. This is used for the following scoring
#'   methods: `"ML"`, `"MLF"`, `"WL"`, and `"MAP"`. Default is `c(-5, 5)`.
#' @param norm.prior A numeric vector of length two specifying the mean and
#'   standard deviation of the normal prior distribution. These values are used
#'   to generate the Gaussian quadrature points and weights. Ignored if `method`
#'   is `"ML"`, `"MLF"`, `"WL"`, or `"INV.TCC"`. Default is `c(0, 1)`.
#' @param nquad An integer indicating the number of Gaussian quadrature points
#'   to be generated from the normal prior distribution. Used only when `method`
#'   is `"EAP"` or `"EAP.SUM"`. Ignored for `"ML"`, `"MLF"`, `"WL"`, `"MAP"`,
#'   and `"INV.TCC"`. Default is `41`.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and their corresponding weights (in the second
#'   column) for the latent variable prior distribution. The weights and points
#'   can be conveniently generated using the function [irtQ::gen.weight()].
#'
#'   If `NULL` and `method` is either `"EAP"` or `"EAP.SUM"`, default quadrature
#'   values are generated based on the `norm.prior` and `nquad` arguments.
#'   Ignored if `method` is `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, or `"INV.TCC"`.
#' @param fence.a A numeric value specifying the item slope parameter (i.e.,
#'   *a*-parameter) for the two imaginary items used in MLF. See **Details** below.
#'   Default is 3.0.
#' @param fence.b A numeric vector of length two specifying the lower and upper
#'   bounds of the item difficulty parameters (i.e., *b*-parameters) for the two
#'   imaginary items in MLF. If `fence.b = NULL`, the values specified in the
#'   `range` argument are used instead. Default is NULL.
#' @param tol A numeric value specifying the convergence tolerance for the ML,
#'   MLF, WL, MAP, and inverse TCC scoring methods. Newton-Raphson optimization
#'   is used for ML, MLF, WL, and MAP, while the bisection method is used for
#'   inverse TCC. Default is 1e-4.
#' @param max.iter A positive integer specifying the maximum number of
#'   iterations allowed for the Newton-Raphson optimization. Default is 100.
#' @param stval.opt A positive integer specifying the starting value option for
#'   the ML, MLF, WL, and MAP scoring methods. Available options are:
#'   - 1: Brute-force search (default)
#'   - 2: Based on observed sum scores
#'   - 3: Fixed at 0
#'
#'   See **Details** below for more information.
#' @param se Logical. If `TRUE`, standard errors of ability estimates are
#'   computed. If `method` is "EAP.SUM" or "INV.TCC", standard errors are always
#'   returned regardless of this setting. Default is `TRUE`.
#' @param intpol Logical. If `TRUE` and `method = "INV.TCC"`, linear
#'   interpolation is applied to approximate ability estimates for sum scores
#'   that cannot be directly mapped using the TCC (e.g., when the observed sum
#'   score is less than the total of item guessing parameters). Default is
#'   `TRUE`. See **Details** below.
#' @param range.tcc A numeric vector of length two specifying the lower and
#'   upper bounds of ability estimates when `method = "INV.TCC"`. Default is
#'   `c(-7, 7)`.
#' @param missing A value indicating missing responses in the data set. Default
#'   is `NA`. See **Details** below.
#' @param ncore An integer specifying the number of logical CPU cores to use for
#'   parallel processing. Default is 1. See **Details** below.
#' @param ... Additional arguments passed to [parallel::makeCluster()].
#'
#' @details For the MAP scoring method, only a normal prior distribution is
#'   supported for the population distribution.
#'
#'   When there are missing responses in the data set, the missing value must be
#'   explicitly specified using the `missing` argument. Missing data are
#'   properly handled when using the ML, MLF, WL, MAP, or EAP methods. However,
#'   when using the "EAP.SUM" or "INV.TCC" methods, any missing responses are
#'   automatically treated as incorrect (i.e., recoded as 0s).
#'
#'   In the maximum likelihood estimation with fences (MLF; Han, 2016), two
#'   imaginary items based on the 2PL model are introduced. The first imaginary
#'   item functions as the lower fence, and its difficulty parameter (*b*)
#'   should be smaller than any of the difficulty parameters in the test form.
#'   Similarly, the second imaginary item serves as the upper fence, and its *b*
#'   parameter should be greater than any difficulty value in the test form.
#'   Both imaginary items should also have very steep slopes (i.e., high
#'   *a*-parameter values). See Han (2016) for more details. If `fence.b =
#'   NULL`, the function will automatically assign the lower and upper fences
#'   based on the values provided in the `range` argument.
#'
#'   When the "INV.TCC" method is used with the 3PL model, ability estimates
#'   cannot be obtained for observed sum scores that are less than the sum of
#'   the items' guessing parameters. In such cases, linear interpolation can be
#'   applied by setting `intpol = TRUE`.
#'
#'   Let \eqn{\theta_{min}} and \eqn{\theta_{max}} denote the minimum and
#'   maximum ability estimates, respectively, and let \eqn{\theta_{X}} be the
#'   ability estimate corresponding to the smallest observed sum score, X, that
#'   is greater than or equal to the sum of the guessing parameters.When linear
#'   interpolation is applied, the first value in the `range.tcc` argument is
#'   treated as \eqn{\theta_{min}}. A line is then constructed between the
#'   points \eqn{(x = \theta_{min}, y = 0)} and \eqn{(x = \theta_{X}, y = X)}.
#'   The second value in `range.tcc` is interpreted as \eqn{\theta_{max}}, which
#'   corresponds to the ability estimate for the maximum observed sum score.
#'
#'   For the "INV.TCC" method, standard errors of ability estimates are computed
#'   using the approach proposed by Lim et al. (2020). The implementation of
#'   inverse TCC scoring in this function is based on a modified version of the
#'   `SNSequate::irt.eq.tse()` function from the \pkg{SNSequate} package
#'   (González, 2014).
#'
#'   For the ML, MLF, WL, and MAP scoring methods, different strategies can be
#'   used to determine the starting value for ability estimation based on the
#'   `stval.opt` argument:
#'
#'   - When `stval.opt = 1` (default), a brute-force search is performed by
#'   evaluating the log-likelihood at discrete theta values within the range
#'   specified by `range`, using 0.1 increments. The theta value yielding the
#'   highest log-likelihood is chosen as the starting value.
#'
#'   - When `stval.opt = 2`, the starting value is derived from the observed
#'   sum score using a logistic transformation. For example, if the maximum
#'   possible score (`max.score`) is 30 and the examinee’s observed sum score
#'   (`obs.score`) is 20, the starting value is `log(obs.score / (max.score -
#'   obs.score))`.
#'     - If all responses are incorrect (i.e., `obs.score = 0`), the starting
#'   value is `log(1 / max.score)`.
#'     - If all responses are correct (`obs.score = max.score`), the starting
#'   value is `log(max.score / 1)`.
#'
#'   - When `stval.opt = 3`, the starting value is fixed at 0.
#'
#'   To accelerate ability estimation using the ML, MLF, WL, MAP, and EAP
#'   methods, this function supports parallel processing across multiple logical
#'   CPU cores. The number of cores can be specified via the `ncore` argument
#'   (default is 1).
#'
#'   Note that the standard errors of ability estimates are computed based on
#'   the Fisher expected information for the ML, MLF, WL, and MAP methods.
#'
#'   For the implementation of the WL method, the function references the
#'   `catR::Pi()`, `catR::Ji()`, and `catR::Ii()` functions from the \pkg{catR}
#'   package (Magis & Barrada, 2017).
#'
#' @return When `method` is one of `"ML"`, `"MLF"`, `"WL"`, `"MAP"`, or `"EAP"`,
#' a two-column data frame is returned:
#' * Column 1: Ability estimates
#' * Column 2: Standard errors of the ability estimates
#'
#' When `method` is either `"EAP.SUM"` or `"INV.TCC"`, a list with two
#' components is returned:
#' * Object 1: A three-column data frame including:
#'    * Column 1: Observed sum scores
#'    * Column 2: Ability estimates
#'    * Column 3: Standard errors of the ability estimates
#' * Object 2: A score table showing possible raw sum scores and the corresponding
#' ability and standard error estimates
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::simdat()], [irtQ::shape_df()],
#'   [irtQ::gen.weight()]
#'
#' @references Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of
#'   ability in a microcomputer environment. *Psychometrika, 35*, 179-198.
#'
#'   González, J. (2014). SNSequate: Standard and nonstandard statistical models
#'   and methods for test equating.
#' *Journal of Statistical Software, 59*, 1-30.
#'
#'   Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).*Fundamentals of
#'   item response theory*. Newbury Park, CA: Sage.
#'
#'   Han, K. T. (2016). Maximum likelihood score estimation method with fences
#'   for short-length tests and computerized adaptive tests.
#' *Applied psychological measurement, 40*(4), 289-301.
#'
#'   Howard, J. P. (2017). *Computational methods for numerical analysis with
#'   R*. New York: Chapman and Hall/CRC.
#'
#'   Kolen, M. J. & Brennan, R. L. (2004). *Test Equating, Scaling, and Linking*
#'   (2nd ed.). New York: Springer
#'
#'   Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT proficiency
#'   estimates.
#' *Educational Measurement: Issues and Practice, 29*(3), 8-14.
#'
#'   Lim, H., Davey, T., & Wells, C. S. (2020). A recursion-based analytical
#'   approach to evaluate the performance of MST.
#' *Journal of Educational Measurement, 58*(2), 154-178.
#'
#'   Magis, D., & Barrada, J. R. (2017). Computerized adaptive testing with R:
#'   Recent updates of the package catR.
#' *Journal of Statistical Software, 76*, 1-19.
#'
#'   Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
#' *Journal of Educational and Behavioral Statistics, 21*(4), 365-389.
#'
#'   Thissen, D. & Orlando, M. (2001). Item response theory for items scored in
#'   two categories. In D. Thissen & H. Wainer (Eds.),
#' *Test scoring* (pp.73-140). Mahwah, NJ: Lawrence Erlbaum.
#'
#'   Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item
#'   Response Theory for Scores on Tests Including Polytomous Items with Ordered
#'   Responses. *Applied Psychological Measurement, 19*(1), 39-49.
#'
#'   Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#'   response theory. *Psychometrika, 54*(3), 427-450.
#'
#' @examples
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameters and convert them into item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Generate examinee ability values
#' set.seed(12)
#' theta <- rnorm(10)
#'
#' # Simulate item response data based on the item metadata and abilities
#' data <- simdat(x, theta, D = 1)
#'
#' \donttest{
#' # Estimate abilities using maximum likelihood (ML)
#' est_score(x, data, D = 1, method = "ML", range = c(-4, 4), se = TRUE)
#'
#' # Estimate abilities using weighted likelihood (WL)
#' est_score(x, data, D = 1, method = "WL", range = c(-4, 4), se = TRUE)
#'
#' # Estimate abilities using MLF with default fences
#' # based on the `range` argument
#' est_score(x, data, D = 1, method = "MLF",
#'   fence.a = 3.0, fence.b = NULL, se = TRUE)
#'
#' # Estimate abilities using MLF with user-specified fences
#' est_score(x, data, D = 1, method = "MLF", fence.a = 3.0,
#'   fence.b = c(-7, 7), se = TRUE)
#'
#' # Estimate abilities using maximum a posteriori (MAP)
#' est_score(x, data, D = 1, method = "MAP", norm.prior = c(0, 1),
#'   nquad = 30, se = TRUE)
#'
#' # Estimate abilities using expected a posteriori (EAP)
#' est_score(x, data, D = 1, method = "EAP", norm.prior = c(0, 1),
#'   nquad = 30, se = TRUE)
#'
#' # Estimate abilities using EAP summed scoring
#' est_score(x, data, D = 1, method = "EAP.SUM", norm.prior = c(0, 1),
#'   nquad = 30)
#'
#' # Estimate abilities using inverse TCC scoring
#' est_score(x, data, D = 1, method = "INV.TCC", intpol = TRUE,
#'   range.tcc = c(-7, 7))
#' }
#'
#' @export
est_score <- function(x, ...) UseMethod("est_score")

#' @describeIn est_score Default method to estimate examinees' latent ability
#'  parameters using a data frame `x` containing the item metadata.
#' @importFrom reshape2 melt
#' @import dplyr
#' @export
est_score.default <- function(x,
                              data,
                              D = 1,
                              method = "ML",
                              range = c(-5, 5),
                              norm.prior = c(0, 1),
                              nquad = 41,
                              weights = NULL,
                              fence.a = 3.0,
                              fence.b = NULL,
                              tol = 1e-4,
                              max.iter = 100,
                              se = TRUE,
                              stval.opt = 1,
                              intpol = TRUE,
                              range.tcc = c(-7, 7),
                              missing = NA,
                              ncore = 1,
                              ...) {

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


#' @describeIn est_score An object created by the function [irtQ::est_irt()].
#' @import dplyr
#' @export
est_score.est_irt <- function(x,
                              method = "ML",
                              range = c(-5, 5),
                              norm.prior = c(0, 1),
                              nquad = 41,
                              weights = NULL,
                              fence.a = 3.0,
                              fence.b = NULL,
                              tol = 1e-4,
                              max.iter = 100,
                              se = TRUE,
                              stval.opt = 1,
                              intpol = TRUE,
                              range.tcc = c(-7, 7),
                              missing = NA,
                              ncore = 1,
                              ...) {
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
est_score_1core <- function(x,
                            elm_item,
                            data,
                            D = 1,
                            method = "ML",
                            max.cats,
                            max.col,
                            range = c(-4, 4),
                            norm.prior = c(0, 1),
                            nquad = 41,
                            weights = NULL,
                            tol = 1e-4,
                            max.iter = 30,
                            se = TRUE,
                            stval.opt = 1,
                            ji = FALSE) {
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
      # monotonically increasing or decreasing), use the 0
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
