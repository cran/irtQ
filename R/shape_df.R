#' Create a data frame of item metadata
#'
#' @description This function creates a data frame which includes item meta (e.g., item parameter, categories, models ...) to be
#' used for the IRT model-data fit analysis as well as other analyses.
#'
#' @param par.drm A list containing three vectors of dichotomous item parameters. Namely, the item discrimination (a), item difficulty (b),
#' and item guessing parameters.
#' @param par.prm A list containing a vector of polytomous item discrimination (or slope) parameters and a list of polytomous item threshold
#' (or step) parameters. In this list, the argument \code{a} should have a vector of slope parameters and the argument \code{d} should include
#' a list of threshold (or step) parameters. See below for more details.
#' @param item.id A character vector of item IDs. If NULL, an ID is automatically given to each item.
#' @param cats A vector containing the number of score categories for items.
#' @param model A character vector of IRT models corresponding to items. The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for
#' dichotomous items, and "GRM" and "GPCM" for polytomous items. Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and
#' "3PLM") and "GRM" and "GPCM" represent the graded response model and (generalized) partial credit model, respectively.
#' @param default.par A logical value to create an item meta with default item parameters. If TRUE, the number of score categories
#' and corresponding IRT models should be specified in the arguments of \code{cats} and \code{model}, respectively. In the default
#' item meta, the item slope parameter has a fixed value of 1, the item difficulty (or threshold) parameter(s) has(have) a fixed value of 0,
#' and the item guessing parameter has a fixed value of .2. Default is FALSE.
#'
#' @details For any item where "1PLM" or "2PLM" is specified in \code{model}, the item guessing parameter will be NA. If \code{model} is
#' a vector of \eqn{length = 1}, the specified model is replicated across all items. As in the function \code{\link{simdat}}, it is important
#' to clearly specify \code{cats} according to the order of items in the test form when a data frame for a mixed-format test needs to be created.
#' See \code{\link{simdat}} for more details about how to specify \code{cats}.
#'
#' When specifying item parameters in \code{par.drm} and/or \code{par.prm}, keep the order of item parameter types. For example,
#' in the \code{par.drm} argument, the first argument \code{a} should contain the slope parameter vector, the second argument \code{b}
#' should contain the difficulty vector, and the third argument \code{g} should contain the guessing parameter vector.
#' In the \code{par.drm} argument, the first argument \code{a} should contain the slope parameter vector and the second argument \code{d}
#' should contain a list including vectors of item threshold (or step) parameters for polytomous response IRT models. Note that when an item follows
#' the (generalized) partial credit model, the item step parameters are the overall item difficulty (or location) parameter subtracted by
#' the difficulty (or threshold) parameter for each category. Thus, the number of step parameters for item with m categories is m-1 because
#' a step parameter for the first category does not affect the category probabilities.
#'
#' @return This function returns a data frame.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{info}}
#'
#' @examples
#' ## a mixed-item format test form
#' ## with five dichotomous and two polytomous items
#' # create a list containing the dichotomous item parameters
#' par.drm <- list(
#'   a = c(1.1, 1.2, 0.9, 1.8, 1.4),
#'   b = c(0.1, -1.6, -0.2, 1.0, 1.2),
#'   g = rep(0.2, 5)
#' )
#'
#' # create a list containing the polytomous item parameters
#' par.prm <- list(
#'   a = c(1.4, 0.6),
#'   d = list(
#'     c(0.0, -1.9, 1.2),
#'     c(0.4, -1.1, 1.5, 0.2)
#'   )
#' )
#'
#' # create a numeric vector of score categories for the items
#' cats <- c(2, 4, 2, 2, 5, 2, 2)
#'
#' # create a character vector of IRT models for the items
#' model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")
#'
#' # create an item meta set
#' shape_df(par.drm = par.drm, par.prm = par.prm, cats = cats, model = model)
#'
#' ## an empty item meta with five dichotomous and two polytomous items
#' # create a numeric vector of score categories for the items
#' cats <- c(2, 4, 3, 2, 5, 2, 2)
#'
#' # create a character vector of IRT models for the items
#' model <- c("1PLM", "GRM", "GRM", "2PLM", "GPCM", "DRM", "3PLM")
#'
#' # create an empty item meta set
#' shape_df(cats = cats, model = model, default.par = TRUE)
#'
#' ## an item meta for a single-item format test form with five dichotomous
#' shape_df(par.drm = par.drm, cats = rep(2, 5), model = "DRM")
#'
#' @export
#'
shape_df <- function(par.drm = list(a = NULL, b = NULL, g = NULL), par.prm = list(a = NULL, d = NULL),
                     item.id = NULL, cats, model, default.par = FALSE) {
  # ensure that the model names are all upper cases
  model <- toupper(model)

  # check model names
  if (!all(model %in% c("1PLM", "2PLM", "3PLM", "DRM", "GRM", "GPCM"))) {
    stop(paste0(
      "At least, one model name is mis-specified in the model argument. \n",
      "Available model names are 1PLM, 2PLM, 3PLM, DRM, GRM, and GPCM"
    ), call. = FALSE)
  }

  # only to create an empty item meta
  if (default.par) {
    if (missing(cats) | missing(model)) {
      stop("The number of score categories and IRT models must be specified.", call. = FALSE)
    }

    # find the index of drm items
    idx.drm <- which(cats == 2)
    if (sum(idx.drm) == 0) idx.drm <- NULL

    # find the index of prm items
    idx.prm <- which(cats > 2)
    if (sum(idx.prm) == 0) idx.prm <- NULL

    # assign default values to the item parameters
    if (!is.null(idx.drm)) {
      par.drm <- list(
        a = rep(1, length(idx.drm)),
        b = rep(0, length(idx.drm)), g = rep(NA_real_, length(idx.drm))
      )
      par.drm$g[model[idx.drm] == "3PLM" | model[idx.drm] == "DRM"] <- 0.2
    } else {
      par.drm <- list(a = NULL, b = NULL, g = NULL)
    }
    if (!is.null(idx.prm)) {
      par.prm <- list(a = rep(1, length(idx.prm)), d = vector("list", length(idx.prm)))
      for (i in 1:length(idx.prm)) {
        par.prm$d[[i]] <- rep(0, cats[idx.prm[i]] - 1)
      }
    } else {
      par.prm <- list(a = NULL, d = NULL)
    }
  }

  # number of the items
  nitem <- length(par.drm$b) + length(par.prm$d)

  # max score categories
  max.cat <- max(cats)

  # create a vector of item ids when item.id = NULL
  if (is.null(item.id)) item.id <- paste0("V", 1:nitem)

  # create a vector of score categories when length(cats) = 1
  if (length(cats) == 1) cats <- rep(cats, nitem)

  # create a vector of model names when length(model) = 1
  if (length(model) == 1) model <- rep(model, nitem)

  # find the index of DRM and PLM items when default.par = FALSE
  if (!default.par) {
    # find the index of drm items
    idx.drm <- which(cats == 2)
    if (sum(idx.drm) == 0) idx.drm <- NULL

    # find the index of prm items
    idx.prm <- which(cats > 2)
    if (sum(idx.prm) == 0) idx.prm <- NULL
  }

  # create an empty matrix to contain item parameters
  if (is.null(idx.prm)) {
    par_mat <- array(NA, c(nitem, 3))
  } else {
    par_mat <- array(NA, c(nitem, max.cat))
  }

  # if drm items exist
  if (!is.null(idx.drm)) {
    if (is.null(par.drm[[3]])) par.drm[[3]] <- rep(0, length(idx.drm))
    par_mat[idx.drm, 1:3] <- bind.fill(par.drm, type = "cbind")
  }

  # if prm items exist
  if (!is.null(idx.prm)) {
    par_mat[idx.prm, 1] <- par.prm[[1]]
    par_mat[idx.prm, 2:max.cat] <- bind.fill(par.prm[[2]], type = "rbind")
  }

  # create an item metadata
  x <- data.frame(id = item.id, cats = cats, model = model, par_mat, stringsAsFactors = FALSE)

  # re-assign column names
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # assign NAs to the par.3 column for the 1PLM, 2PLM, items
  x[x$model %in% c("1PLM", "2PLM"), "par.3"] <- NA_real_

  # assign 0s to the par.3 column for the 3PLM when par.3 = NA
  x[x$model == "3PLM" & is.na(x$par.3), "par.3"] <- 0

  # last check
  if (any(x[x$cats == 2, 3] %in% c("GRM", "GPCM"))) {
    stop("Dichotomous items must have models among '1PLM', '2PLM', '3PLM', and 'DRM'.", call. = FALSE)
  }
  if (any(x[x$cats > 2, 3] %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
    stop("Polytomous items must have models among 'GRM' and 'GPCM'.", call. = FALSE)
  }

  # return the results
  x
}


# This function creates an item meta containing the starting values
startval_df <- function(cats, model, item.id = NULL) {
  # ensure that the model names are all upper cases
  model <- toupper(model)

  # check model names
  if (!all(model %in% c("1PLM", "2PLM", "3PLM", "DRM", "GRM", "GPCM"))) {
    stop(paste0(
      "At least, one model name is mis-specified in the model argument. \n",
      "Available model names are 1PLM, 2PLM, 3PLM, DRM, GRM, and GPCM"
    ), call. = FALSE)
  }

  # find the index of drm items
  idx.drm <- which(cats == 2)
  if (sum(idx.drm) == 0) idx.drm <- NULL

  # find the index of prm items
  idx.prm <- which(cats > 2)
  if (sum(idx.prm) == 0) idx.prm <- NULL

  # assign default values to the item parameters
  if (!is.null(idx.drm)) {
    par.drm <- list(
      a = rep(1, length(idx.drm)),
      b = rep(0, length(idx.drm)), g = rep(0, length(idx.drm))
    )
    par.drm$g[model[idx.drm] == "3PLM" | model[idx.drm] == "DRM"] <- 0.2
  } else {
    par.drm <- list(a = NULL, b = NULL, g = NULL)
  }
  if (!is.null(idx.prm)) {
    par.prm <- list(a = rep(1, length(idx.prm)), d = vector("list", length(idx.prm)))
    for (i in 1:length(idx.prm)) {
      par.prm$d[[i]] <- seq(-1.0, 1.0, length.out = (cats[idx.prm[i]] - 1))
    }
  } else {
    par.prm <- list(a = NULL, d = NULL)
  }

  # number of the items
  nitem <- length(par.drm$b) + length(par.prm$d)

  # max score categories
  max.cat <- max(cats)

  # create a vector of item ids when item.id = NULL
  if (is.null(item.id)) item.id <- paste0("V", 1:nitem)

  # create a vector of score categories when length(cats) = 1
  if (length(cats) == 1) cats <- rep(cats, nitem)

  # create a vector of model names when length(model) = 1
  if (length(model) == 1) model <- rep(model, nitem)

  # create an empty matrix to contain item parameters
  if (is.null(idx.prm)) {
    par_mat <- array(NA, c(nitem, 3))
  } else {
    par_mat <- array(NA, c(nitem, max.cat))
  }

  # if drm items exist
  if (!is.null(idx.drm)) {
    if (is.null(par.drm[[3]])) par.drm[[3]] <- rep(0, length(idx.drm))
    par_mat[idx.drm, 1:3] <- bind.fill(par.drm, type = "cbind")
  }

  # if prm items exist
  if (!is.null(idx.prm)) {
    par_mat[idx.prm, 1] <- par.prm[[1]]
    par_mat[idx.prm, 2:max.cat] <- bind.fill(par.prm[[2]], type = "rbind")
  }

  # create an item metadata
  x <- data.frame(id = item.id, cats = cats, model = model, par_mat, stringsAsFactors = FALSE)

  # re-assign column names
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # assign NAs to the par.3 column for the 1PLM, 2PLM, items
  # x[x$model %in% c("1PLM", "2PLM"), "par.3"] <- NA_real_

  # assign 0s to the par.3 column for the 3PLM when par.3 = NA
  x[x$model == "3PLM" & is.na(x$par.3), "par.3"] <- 0

  # last check
  if (any(x[x$cats == 2, 3] %in% c("GRM", "GPCM"))) {
    stop("Dichotomous items must have models among '1PLM', '2PLM', '3PLM', and 'DRM'.", call. = FALSE)
  }
  if (any(x[x$cats > 2, 3] %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
    stop("Polytomous items must have models among 'GRM' and 'GPCM'.", call. = FALSE)
  }

  # return the results
  x
}
