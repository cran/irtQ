#' Fixed ability parameter calibration
#'
#' @description This function performs the fixed ability parameter calibration (FAPC), often called
#' Method A, which is the maximum likelihood estimation of item parameters given the ability
#' estimates (Baker & Kim, 2004; Ban, Hanson, Wang, Yi, & Harris, 2001; Stocking, 1988). Also, this could be
#' considered as a special type of the joint maximum likelihood estimation where only one cycle of parameter
#' estimation is implemented given the ability estimates (Birnbaum, 1968). FAPC is one of potentially useful
#' online item calibration methods for computerized adaptive testing (CAT) to put the parameter estimates of
#' pretest items on the same scale of operational item parameter estimates and recalibrate the operational
#' items to evaluate the parameter drifts of the operational items (Chen & Wang, 2016; Stocking, 1988).
#'
#' @param x A data frame containing the item metadata. This metadata is necessary to obtain the information of
#' each item (i.e., number of score categories and IRT model) to be calibrated. You can easily create an empty
#' item metadata using the function \code{\link{shape_df}}. When \code{use.startval = TRUE}, the item parameters
#' specified in the item metadata are used as the starting values for the item parameter estimation.
#' If \code{x = NULL}, the arguments of \code{model} and \code{cats} must be specified. See \code{\link{irtfit}},
#' \code{\link{info}} or \code{\link{simdat}} for more details about the item metadata. Default is NULL.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param score A vector of examinees' ability estimates. Length of the vector must be the same as the number of rows in the
#' response data set.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param model A vector of character strings indicating what IRT model is used to calibrate each item. Available IRT models are
#' "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items. "GRM" and "GPCM" represent the graded
#' response model and (generalized) partial credit model, respectively. Note that "DRM" is considered as "3PLM" in this function.
#' If a single character of the IRT model is specified, that model will be recycled across all items. The provided information in the \code{model}
#' argument is used only when \code{x = NULL}. Default is NULL.
#' @param cats A numeric vector specifying the number of score categories for each item. For example, a dichotomous
#' item has two score categories. If a single numeric value is specified, that value will be recycled across all items. If \code{cats = NULL}
#' and all specified models in the \code{model} argument are the dichotomous models (i.e., 1PLM, 2PLM, 3PLM, or DRM), it assumes
#' that all items have two score categories. The provided information in the \code{cats} argument is used only when \code{x = NULL}.
#' Default is NULL.
#' @param item.id A character vector of item IDs. If NULL, the item IDs are generated automatically. Default is NULL.
#' @param fix.a.1pl A logical value. If TRUE, the slope parameters of the 1PLM items are fixed to a specific value specified in the argument
#' \code{a.val.1pl}. Otherwise, the slope parameters of all 1PLM items are constrained to be equal and estimated. Default is FALSE.
#' @param fix.a.gpcm A logical value. If TRUE, the GPCM items are calibrated with the partial credit model and the slope parameters of
#' the GPCM items are fixed to a specific value specified in the argument \code{a.val.gpcm}. Otherwise, the slope parameter of each GPCM item
#' is estimated. Default is FALSE.
#' @param fix.g A logical value. If TRUE, the guessing parameters of the 3PLM items are fixed to a specific value specified in the argument
#' \code{g.val}. Otherwise, the guessing parameter of each 3PLM item is estimated. Default is FALSE.
#' @param a.val.1pl A numeric value. This value is used to fixed the slope parameters of the 1PLM items.
#' @param a.val.gpcm A numeric value. This value is used to fixed the slope parameters of the GPCM items.
#' @param g.val A numeric value. This value is used to fixed the guessing parameters of the 3PLM items.
#' @param use.aprior A logical value. If TRUE, a prior distribution for the slope parameters is used for the parameter calibration
#' across all items. Default is FALSE.
#' @param use.bprior A logical value. If TRUE, a prior distribution for the difficulty (or threshold) parameters is used for the parameter calibration
#' across all items. Default is FALSE.
#' @param use.gprior A logical value. If TRUE, a prior distribution for the guessing parameters is used for the parameter calibration
#' across all 3PLM items. Default is TRUE.
#' @param aprior A list containing the information of the prior distribution for item slope parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param bprior A list containing the information of the prior distribution for item difficulty (or threshold) parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param gprior A list containing the information of the prior distribution for item guessing parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param use.startval A logical value. If TRUE, the item parameters provided in the item metadata (i.e., the argument \code{x}) are used as
#' the starting values for the item parameter estimation. Otherwise, internal starting values of this function are used. Default is FALSE.
#' @param control A list of control parameters to be passed to the optimization function of \code{nlminb()} in the \pkg{stats} package. The control parameters
#' set the conditions of the item parameter estimation process such as the maximum number of iterations. See \code{nlminb()} in the \pkg{stats} package for details.
#' @param verbose A logical value. If FALSE, all progress messages are suppressed. Default is TRUE.
#'
#' @details In most cases, the function \code{\link{est_item}} will return successfully converged item parameter estimates using
#' the default internal starting values. However, if there is a convergence problem in the calibration, one possible solution is using
#' different starting values. When the item parameter values are specified in the item metadata (i.e., the argument \code{x}), those values
#' can be used as the starting values for the item parameter calibration by setting \code{use.startval = TRUE}.
#'
#' @return This function returns an object of class \code{\link{est_item}}. Within this object, several internal objects are contained such as:
#' \item{estimates}{A data frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
#' \item{par.est}{A data frame containing the item parameter estimates.}
#' \item{se.est}{A data frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
#' observed information functions.}
#' \item{pos.par}{A data frame containing the position number of item parameters being estimated. The position information is useful
#' when interpreting the variance-covariance matrix of item parameter estimates.}
#' \item{covariance}{A matrix of variance-covariance matrix of item parameter estimates.}
#' \item{loglikelihood}{A sum of the log-likelihood values of the complete data set across all estimated items.}
#' \item{data}{A data frame of the examinees' response data set.}
#' \item{score}{A vector of the examinees' ability values used as the fixed effects.}
#' \item{scale.D}{A scaling factor in IRT models.}
#' \item{convergence}{A string indicating the convergence status of the item parameter estimation.}
#' \item{nitem}{A total number of items included in the response data.}
#' \item{deleted.item}{The items which have no item response data. Those items are excluded from the item parameter estimation.}
#' \item{npar.est}{A total number of the estimated parameters.}
#' \item{n.response}{An integer vector indicating the number of item responses for each item used to estimate the item parameters.}
#' \item{TotalTime}{Time (in seconds) spent for total compuatation.}
#'
#' The internal objects can be easily extracted using the function \code{\link{getirt}}.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{irtfit}}, \code{\link{info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{sx2_fit}},
#' \code{\link{traceline.est_item}}, \code{\link{getirt}}
#'
#' @references
#' Baker, F. B., & Kim, S. H. (2004). \emph{Item response theory: Parameter estimation techniques.} CRC Press.
#'
#' Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A comparative study of on-line pretest item calibration/scaling methods
#' in computerized adaptive testing. \emph{Journal of Educational Measurement, 38}(3), 191-212.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In F. M. Lord & M. R. Novick (Eds.),
#' \emph{Statistical theories of mental test scores} (pp. 397-479). Reading, MA: Addison-Wesley.
#'
#' Chen, P., & Wang, C. (2016). A new online calibration method for multidimensional computerized adaptive testing.
#' \emph{Psychometrika, 81}(3), 674-701.
#'
#' Stocking, M. L. (1988). \emph{Scale drift in on-line calibration} (Research Rep. 88-28). Princeton, NJ: ETS.
#'
#' @examples
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # select the item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # modify the item metadata so that some items follow 1PLM, 2PLM and GPCM
#' x[c(1:3, 5), 3] <- "1PLM"
#' x[c(1:3, 5), 4] <- 1
#' x[c(1:3, 5), 6] <- 0
#' x[c(4, 8:12), 3] <- "2PLM"
#' x[c(4, 8:12), 6] <- 0
#' x[54:55, 3] <- "GPCM"
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(500, mean = 0, sd = 1)
#'
#' # simulate the response data
#' data <- simdat(x = x, theta = score, D = 1)
#'
#' \donttest{
#' # 1) item parameter estimation: constrain the slope parameters of the 1PLM to be equal
#' (mod1 <- est_item(x, data, score,
#'   D = 1, fix.a.1pl = FALSE, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 17)), use.startval = FALSE
#' ))
#' summary(mod1)
#'
#' # extract the item parameter estimates
#' getirt(mod1, what = "par.est")
#'
#' # 2) item parameter estimation: fix the slope parameters of the 1PLM to 1
#' (mod2 <- est_item(x, data, score,
#'   D = 1, fix.a.1pl = TRUE, a.val.1pl = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 17)), use.startval = FALSE
#' ))
#' summary(mod2)
#'
#' # extract the standard error estimates
#' getirt(mod2, what = "se.est")
#'
#' # 3) item parameter estimation: fix the guessing parameters of the 3PLM to 0.2
#' (mod3 <- est_item(x, data, score,
#'   D = 1, fix.a.1pl = TRUE, fix.g = TRUE, a.val.1pl = 1, g.val = .2,
#'   use.startval = FALSE
#' ))
#' summary(mod3)
#'
#' # extract both item parameter and standard error estimates
#' getirt(mod2, what = "estimates")
#' }
#'
#' @importFrom Matrix sparseMatrix
#' @import dplyr
#' @export
#'
est_item <- function(x = NULL, data, score, D = 1, model = NULL, cats = NULL, item.id = NULL, fix.a.1pl = FALSE, fix.a.gpcm = FALSE,
                     fix.g = FALSE, a.val.1pl = 1, a.val.gpcm = 1, g.val = .2, use.aprior = FALSE, use.bprior = FALSE, use.gprior = TRUE,
                     aprior = list(dist = "lnorm", params = c(0, 0.5)), bprior = list(dist = "norm", params = c(0.0, 1.0)),
                     gprior = list(dist = "beta", params = c(5, 17)), missing = NA, use.startval = FALSE,
                     control = list(eval.max = 500, iter.max = 500), verbose = TRUE) {
  # check start time
  start.time <- Sys.time()

  # match.call
  cl <- match.call()

  ## -------------------------------------------------------------------------------------------------------
  ## 1. preparation of item parameter estimation
  if (verbose) {
    cat("Starting...", "\n")
  }

  # check if the starting values are available
  if (use.startval & is.null(x)) {
    stop(paste0(
      "To use starting values for item parameter estimation, \n",
      "the item metadata must be specified in the argument 'x'."
    ), call. = FALSE)
  }

  # transform a data set to data.frame
  if (is.vector(data)) {
    data <- data.frame(X1 = data, stringsAsFactors = FALSE)
  } else {
    data <- data.frame(data, stringsAsFactors = FALSE)
  }

  # extract information about the number of score categories and models
  if (verbose) {
    cat("Parsing input...", "\n")
  }
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # extract score categories for all items
    cats <- x$cats

    # extract model names for all items
    model <- x$model
  } else {
    # make the model names as upper cases
    model <- toupper(model)

    # when a character string scalar is provided in the model argument
    # create a vector of the same model names
    if (length(model) == 1) {
      model <- rep(model, ncol(data))
    }

    # check if the score category information is provided in the cats argument
    if (is.null(cats)) {
      if (all(model %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
        cats <- rep(2, ncol(data))
      } else {
        stop("The number of score categories for the items should be specified in the argument 'cats'.", call. = FALSE)
      }
    }

    # when an integer scalar is provided in the cats argument
    # create a vector of the same category numbers
    if (length(cats) == 1) {
      cats <- rep(cats, ncol(data))
    }

    # create an empty item metadata data frame
    x <- shape_df(item.id = item.id, cats = cats, model = model, default.par = TRUE)

    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # extract score categories for all items
    cats <- x$cats

    # extract model names for all items
    model <- x$model
  }

  # check the total number of item in the response data set
  nitem <- ncol(data)

  # check the total number of examinees
  nstd <- nrow(data)

  # check whether included data are correct
  if (nrow(x) != nitem) stop("The number of items included in 'x' and 'data' must be the same.", call. = FALSE)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of item responses across all items
  n.resp <- colSums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)

  # delete items which have all missing responses from the item metadata set
  if (length(loc_allmiss) > 0L) {
    loc_nomiss <- which(n.resp > 0L)
    x_pre <- x
    x <- x[-loc_allmiss, ]
    model <- x[, 3]
    cats <- x[, 2]
    data <- data[, loc_nomiss]
    memo2 <- paste0(
      paste0("item ", loc_allmiss, collapse = ", "),
      " is/are excluded in the item parameter estimation because the item(s) has/have no item response data."
    )
    warning(memo2, call. = FALSE)
  } else {
    loc_allmiss <- NULL
  }

  # computes the group parameter estimates: mean & var of ability parameters
  group.par <- c(mu = mean(score), sigma = stats::sd(score))

  # transform scores to a vector if the score has other formats
  if (is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # copy scores to theta values
  theta <- score

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of items with other models
  if ("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")

    # count the number of 1PLM items to be constrained
    n.1PLM <- length(loc_1p_const)
  } else {
    loc_1p_const <- NULL
    loc_else <- 1:nrow(x)
    n.1PLM <- NULL
  }

  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  param_loc <- parloc(
    x = x, loc_1p_const = loc_1p_const, loc_else = loc_else,
    fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g
  )

  # factorize the response values
  resp <- purrr::map2(.x = data, .y = cats, .f = function(k, m) factor(k, levels = (seq_len(m) - 1)))

  # create a contingency table of score categories for each item
  # and then, transform the table to a matrix format
  std.id <- 1:nstd
  freq.cat <-
    purrr::map(
      .x = resp,
      .f = function(k) {
        stats::xtabs(~ std.id + k,
          na.action = stats::na.pass, addNA = FALSE
        ) %>%
          # as.numeric() %>%
          matrix(nrow = length(k))
      }
    )

  # delete 'resp' object
  rm(resp, envir = environment(), inherits = FALSE)

  ## -------------------------------------------------------------------------------------------------------
  ## 2. item parameter estimation
  ## ---------------------------------------------------------------
  # Item parameter estimation
  # create empty vectors to contain results
  est_par <- NULL
  est_se <- NULL
  loc_items <- NULL
  convergence <- NULL
  noconv_items <- NULL
  objective <- NULL
  cov_list <- NULL

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # then, divide items into three groups
  # (1) DRM 1PL items with the constrained slop
  # (2) other DRM items
  # (3) PRM items
  if (sum(loc_else) == 0) {
    drm.else <- NULL
  } else {
    drm.else <- loc_else
  }
  idx4est <- list(
    drm.slc = loc_1p_const,
    drm.else = drm.else,
    prm = idx.prm
  )

  # add the rowsum margin to the contingency matrix of the DRM items
  if (!is.null(idx.drm)) {
    freq.cat.drm <-
      purrr::map(
        .x = freq.cat[idx.drm],
        ~ {
          stats::addmargins(A = .x, margin = 2, quiet = TRUE)
        }
      )
    freq.cat[idx.drm] <- freq.cat.drm

    # delete 'freq.cat.drm' object
    rm(freq.cat.drm, envir = environment(), inherits = FALSE)
  }

  # create the lower and upper bounds of the item parameters
  parbd <- lubound(model, cats, n.1PLM, idx4est, fix.a.1pl, fix.g, fix.a.gpcm)

  # start estimation
  if (verbose) {
    cat("Estimating item parameters...", "\n")
  }

  # (1) estimation of DRM items: 1PLM with the constrained slope value
  if (!is.null(loc_1p_const)) {
    # prepare input files to estimate the 1PLM item parameters
    f_i <- r_i <- s_i <- array(0, c(nstd, n.1PLM))
    for (k in 1:n.1PLM) {
      s_i[, k] <- freq.cat[loc_1p_const][[k]][, 1]
      r_i[, k] <- freq.cat[loc_1p_const][[k]][, 2]
      f_i[, k] <- freq.cat[loc_1p_const][[k]][, 3]
    }

    # set the starting values
    if (use.startval) {
      startval <-
        set_startval(
          pars = elm_item$pars, item = loc_1p_const,
          use.startval = TRUE, mod = "1PLM",
          score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
          fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
        )
    } else {
      startval <-
        set_startval(
          pars = elm_item$pars, item = loc_1p_const,
          use.startval = FALSE, mod = "1PLM",
          score.cat = 2, fix.a.1pl = FALSE, fix.g = fix.g,
          fix.a.gpcm = fix.a.gpcm, n.1PLM = n.1PLM
        )
    }

    # bounds of the item parameters
    lower <- parbd$drm.slc$lower
    upper <- parbd$drm.slc$upper

    # parameter estimation
    est <- estimation1(
      f_i = f_i, r_i = r_i, s_i, theta = theta, mod = "1PLM", D = D, nstd = nstd,
      fix.a.1pl = FALSE, n.1PLM = n.1PLM, aprior = aprior, bprior = bprior,
      use.aprior = use.aprior, use.bprior = use.bprior,
      control = control, startval = startval, lower = lower, upper = upper
    )

    # extract the results
    # item parameter estimates
    a <- est$pars[1]
    b <- est$par[-1]
    pars <- purrr::map(1:n.1PLM, .f = function(x) c(a, b[x], NA))
    est_par <- c(est_par, pars)

    # standard errors
    a.se <- est$se[1]
    b.se <- est$se[-1]
    pars.se <- purrr::map(1:n.1PLM, .f = function(x) c(NA, b.se[x], NA))
    pars.se[[1]][1] <- a.se
    est_se <- c(est_se, pars.se)

    # covariance matrix
    cov_list <- c(cov_list, list(est$covariance))

    # convergence indicator
    convergence <- c(convergence, est$convergence)
    if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_1p_const)

    # negative loglikelihood value
    objective <- c(objective, est$objective)
  }

  # all other items
  if (length(loc_else) >= 1) {
    for (i in 1:length(loc_else)) {
      # prepare information to estimate item parameters
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a DRM item
      if (score.cat == 2) {
        s_i <- freq.cat[loc_else][[i]][, 1]
        r_i <- freq.cat[loc_else][[i]][, 2]
        f_i <- freq.cat[loc_else][[i]][, 3]

        # set the starting values
        if (use.startval) {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = TRUE, mod = mod,
              score.cat = score.cat, fix.a.1pl = TRUE, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        } else {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = FALSE, mod = mod,
              score.cat = score.cat, fix.a.1pl = TRUE, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        }

        # bounds of the item parameters
        loc.tmp <- which(idx4est$drm.else == loc_else[i])
        lower <- parbd$drm.else[[loc.tmp]]$lower
        upper <- parbd$drm.else[[loc.tmp]]$upper

        # parameter estimation
        est <- estimation1(
          f_i = f_i, r_i = r_i, s_i = s_i, theta = theta, mod = mod, D = D,
          nstd = nstd, fix.a.1pl = ifelse(mod == "1PLM", TRUE, FALSE),
          fix.g = fix.g, a.val.1pl = a.val.1pl, g.val = g.val, n.1PLM = NULL,
          aprior = aprior, bprior = bprior, gprior = gprior,
          use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
          control = control, startval = startval, lower = lower, upper = upper
        )

        # extract the results
        # item parameter estimates
        a <- ifelse(mod == "1PLM", a.val.1pl, est$pars[1])
        b <- ifelse(mod == "1PLM", est$pars[1], est$pars[2])
        g <- ifelse(mod == "3PLM", ifelse(fix.g, g.val, est$pars[3]), NA)
        pars <- c(a, b, g)
        est_par <- c(est_par, list(pars))

        # standard errors
        a.se <- ifelse(mod == "1PLM", NA, est$se[1])
        b.se <- ifelse(mod == "1PLM", est$se[1], est$se[2])
        g.se <- ifelse(mod == "3PLM", ifelse(fix.g, NA, est$se[3]), NA)
        pars.se <- c(a.se, b.se, g.se)
        est_se <- c(est_se, list(pars.se))

        # covariance matrix
        cov_list <- c(cov_list, list(est$covariance))

        # convergence indicator
        convergence <- c(convergence, est$convergence)
        if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        # negative loglikelihood value
        objective <- c(objective, est$objective)
      }

      # in case of a PRM item
      if (score.cat > 2) {
        r_i <- freq.cat[loc_else][[i]]

        # set the starting values
        if (use.startval) {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = TRUE, mod = mod,
              score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        } else {
          startval <-
            set_startval(
              pars = elm_item$pars, item = loc_else[i],
              use.startval = FALSE, mod = mod,
              score.cat = score.cat, fix.a.1pl = fix.a.1pl, fix.g = fix.g,
              fix.a.gpcm = fix.a.gpcm, n.1PLM = NULL
            )
        }

        # bounds of the item parameters
        loc.tmp <- which(idx4est$prm == loc_else[i])
        lower <- parbd$prm[[loc.tmp]]$lower
        upper <- parbd$prm[[loc.tmp]]$upper

        # parameter estimation
        est <- estimation1(
          r_i = r_i, theta = theta, mod = mod, score.cat = score.cat, D = D,
          nstd = nstd, fix.a.gpcm = ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
          a.val.gpcm = a.val.gpcm, n.1PLM = NULL, aprior = aprior, bprior = bprior,
          use.aprior = use.aprior, use.bprior = use.bprior,
          control = control, startval = startval, lower = lower, upper = upper
        )

        # extract the results
        # item parameter estimates
        a <- ifelse(mod == "GRM", est$pars[1], ifelse(fix.a.gpcm, a.val.gpcm, est$pars[1]))
        if (mod == "GRM") {
          bs <- est$pars[-1]
        } else {
          if (fix.a.gpcm) {
            bs <- est$pars
          } else {
            bs <- est$pars[-1]
          }
        }
        pars <- c(a, bs)
        est_par <- c(est_par, list(pars))

        # standard errors
        a.se <- ifelse(mod == "GRM", est$se[1], ifelse(fix.a.gpcm, NA, est$se[1]))
        if (mod == "GRM") {
          bs.se <- est$se[-1]
        } else {
          if (fix.a.gpcm) {
            bs.se <- est$se
          } else {
            bs.se <- est$se[-1]
          }
        }
        pars.se <- c(a.se, bs.se)
        est_se <- c(est_se, list(pars.se))

        # covariance matrix
        cov_list <- c(cov_list, list(est$covariance))

        # convergence indicator
        convergence <- c(convergence, est$convergence)
        if (est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        # negative loglikelihood value
        objective <- c(objective, est$objective)
      }
    }
  }

  ## ---------------------------------------------------------------
  # check the convergence of parameter estimation
  if (sum(convergence) == 0L) {
    note <- "All item parameters were successfully converged."
  } else {
    note <- paste0(paste0("Item ", sort(noconv_items), collapse = ", "), " was(were) not successfully converged.")
  }

  ## ---------------------------------------------------------------
  # compute the sum of loglikelihood values
  llike <- -sum(objective)

  ## ---------------------------------------------------------------
  # arrange the estimated item parameters and standard errors
  par_df <- data.frame(bind.fill(est_par, type = "rbind"))
  par_df$loc <- c(loc_1p_const, loc_else)
  par_df <-
    par_df %>%
    dplyr::arrange("loc") %>%
    dplyr::select(-"loc")
  se_df <- data.frame(bind.fill(est_se, type = "rbind"))
  se_df$loc <- c(loc_1p_const, loc_else)
  se_df <-
    se_df %>%
    dplyr::arrange("loc") %>%
    dplyr::select(-"loc")

  # combine all covariance matrices
  cov_mat <- as.matrix(Matrix::bdiag(cov_list))

  # relocate the matrix elements according to the parameter location
  pos_df <- param_loc$loc.par
  reloc.par <- param_loc$reloc.par
  cov_mat2 <-
    Matrix::sparseMatrix(
      i = c(row(cov_mat)[reloc.par, ]), j = c(col(cov_mat)[, reloc.par]), x = c(cov_mat)
    ) %>%
    as.matrix()

  # arrange the estimated item parameters and standard errors
  # when there is any item which have all missing response data
  if (length(loc_allmiss) > 0L) {
    tmp <-
      x_pre[, -c(1:3)] %>%
      dplyr::mutate_all(.funs = function(x) NA) %>%
      purrr::map(.x = 1:nrow(x_pre), .f = function(i, .data) as.numeric(.data[i, ]))
    par_tmp <- se_tmp <- pos_tmp <- tmp
    par_tmp[loc_nomiss] <-
      par_df %>%
      purrr::map(.x = 1:nrow(par_df), .f = function(i, .data) as.numeric(.data[i, ]))
    se_tmp[loc_nomiss] <-
      se_df %>%
      purrr::map(.x = 1:nrow(se_df), .f = function(i, .data) as.numeric(.data[i, ]))
    pos_tmp[loc_nomiss] <-
      pos_df %>%
      purrr::map(.x = 1:nrow(pos_df), .f = function(i, .data) as.numeric(.data[i, ]))
    par_df <- data.frame(bind.fill(par_tmp, type = "rbind"))
    se_df <- data.frame(bind.fill(se_tmp, type = "rbind"))
    pos_df <- data.frame(bind.fill(pos_tmp, type = "rbind"))
  }

  # create a full data.frame for the item parameter estimates
  if (length(loc_allmiss) > 0L) {
    x <- x_pre
  }
  full_par_df <-
    data.frame(x[, 1:3], par_df) %>%
    confirm_df(g2na = TRUE)

  # create a full data.frame for the standard error estimates
  full_se_df <-
    data.frame(x[, 1:3], se_df) %>%
    confirm_df(g2na = TRUE)

  # create a full data.frame containing the position of item parameter estimates
  # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
  full_pos_df <-
    data.frame(full_par_df[, 1:3], pos_df) %>%
    confirm_df(g2na = TRUE)

  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow = nrow(x), ncol = 2 * ncol(par_df)))
  all_df[, seq(1, 2 * ncol(par_df), 2)] <- par_df
  all_df[, seq(2, 2 * ncol(par_df), 2)] <- se_df
  col.names <- rep(NA, 2 * ncol(par_df))
  col.names[seq(1, 2 * ncol(par_df), 2)] <- paste0("par.", 1:ncol(par_df))
  col.names[seq(2, 2 * ncol(par_df), 2)] <- paste0("se.", 1:ncol(se_df))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(full_par_df[, 1:3], all_df)

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  # return results
  rst <- structure(
    list(
      estimates = full_all_df, par.est = full_par_df, se.est = full_se_df, pos.par = full_pos_df, covariance = cov_mat2,
      loglikelihood = llike, group.par = group.par, data = data, score = score, scale.D = D, convergence = note, nitem = nitem,
      deleted.item = as.numeric(loc_allmiss), npar.est = length(param_loc$reloc.par),
      n.response = as.numeric(n.resp), TotalTime = est_time
    ),
    class = "est_item"
  )
  rst$call <- cl

  if (verbose) {
    cat("Estimation is finished.", "\n")
  }
  return(rst)
}
