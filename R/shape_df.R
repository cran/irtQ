#' Create a Data Frame of Item Metadata
#'
#' This function creates a data frame of item metadata—including item
#' parameters, the number of score categories, and IRT model specifications—to
#' be used in various IRT-related analyses within the \pkg{irtQ} package.
#'
#' @param par.drm A list containing three numeric vectors for dichotomous item
#'   parameters: item discrimination (`a`), item difficulty (`b`), and guessing
#'   parameters (`g`).
#' @param par.prm A list containing polytomous item parameters. The list must
#'   include a numeric vector `a` for item discrimination (slope) parameters,
#'   and a list `d` of numeric vectors specifying difficulty (or threshold)
#'   parameters for each item. See the **Details** section for more information.
#' @param item.id A character vector of item IDs. If `NULL`, default IDs (e.g.,
#'   "V1", "V2", ...) are assigned automatically.
#' @param cats A numeric vector indicating the number of score categories for
#'   each item.
#' @param model A character vector specifying the IRT model for each item.
#'   Available options are `"1PLM"`, `"2PLM"`, `"3PLM"`, and `"DRM"` for
#'   dichotomous items, and `"GRM"` and `"GPCM"` for polytomous items. The label
#'   `"DRM"` serves as a general category that encompasses all dichotomous
#'   models (`"1PLM"`, `"2PLM"`, and `"3PLM"`), while `"GRM"` and `"GPCM"` refer
#'   to the graded response model and (generalized) partial credit model,
#'   respectively.
#' @param default.par Logical. If `TRUE`, default item parameters are generated
#'   based on the specified `cats` and `model`. In this case, the slope
#'   parameter is set to 1, all difficulty (or threshold) parameters are set to
#'   0, and the guessing parameter is set to 0.2 for `"3PLM"` or `"DRM"` items.
#'   The default is `FALSE`.
#'
#' @details For any item where `"1PLM"` or `"2PLM"` is specified in `model`, the
#'   guessing parameter will be set to `NA`. If `model` is a vector of length 1,
#'   the specified model will be replicated across all items.
#'
#'   As in the [irtQ::simdat()] function, when constructing a mixed-format test
#'   form, it is important to specify the `cats` argument to reflect the correct
#'   number of score categories for each item, in the exact order that the items
#'   appear. See [irtQ::simdat()] for further guidance on how to specify `cats`.
#'
#'   When specifying item parameters using `par.drm` and/or `par.prm`, the
#'   internal structure and ordering of elements must be followed.
#'   - `par.drm` should be a list with three components:
#'     - `a`: a numeric vector of slope parameters
#'     - `b`: a numeric vector of difficulty parameters
#'     - `g`: a numeric vector of guessing parameters
#'   - `par.prm` should be a list with two components:
#'     - `a`: a numeric vector of slope parameters for polytomous items
#'     - `d`: a list of numeric vectors specifying threshold (or step) parameters
#'       for each polytomous item
#'
#'
#' For items following the (generalized) partial credit model (`"GPCM"`), the
#' threshold (or step) parameters are computed as the overall item difficulty
#' (location) minus the category-specific thresholds. Therefore, for an item
#' with `m` score categories, `m - 1` step parameters must be provided, since
#' the first category threshold is fixed and does not contribute to category
#' probabilities.
#'
#' @return A data frame containing item metadata, including item IDs, number of
#'   score categories, IRT model types, and associated item parameters. This
#'   data frame can be used as input for other functions in the \pkg{irtQ}
#'   package, such as [irtQ::est_irt()] or [irtQ::simdat()].
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::simdat()], [irtQ::shape_df_fipc()]
#'
#' @examples
#' ## A mixed-format test form
#' ## containing five dichotomous items and two polytomous items
#' # Create a list of dichotomous item parameters
#' par.drm <- list(
#'   a = c(1.1, 1.2, 0.9, 1.8, 1.4),
#'   b = c(0.1, -1.6, -0.2, 1.0, 1.2),
#'   g = rep(0.2, 5)
#' )
#'
#' # Create a list of polytomous item parameters
#' par.prm <- list(
#'   a = c(1.4, 0.6),
#'   d = list(
#'     c(0.0, -1.9, 1.2),
#'     c(0.4, -1.1, 1.5, 0.2)
#'   )
#' )
#'
#' # Create a numeric vector indicating the number of score categories for each item
#' cats <- c(2, 4, 2, 2, 5, 2, 2)
#'
#' # Create a character vector specifying the IRT model for each item
#' model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")
#'
#' # Generate an item metadata set using the specified parameters
#' shape_df(par.drm = par.drm, par.prm = par.prm, cats = cats, model = model)
#'
#' ## An empty item metadata frame with five dichotomous items and two polytomous items
#' # Create a numeric vector indicating the number of score categories for each item
#' cats <- c(2, 4, 3, 2, 5, 2, 2)
#'
#' # Create a character vector specifying the IRT model for each item
#' model <- c("1PLM", "GRM", "GRM", "2PLM", "GPCM", "DRM", "3PLM")
#'
#' # Generate an item metadata frame with default parameters
#' shape_df(cats = cats, model = model, default.par = TRUE)
#'
#' ## A single-format test form consisting of five dichotomous items
#' # Generate the item metadata
#' shape_df(par.drm = par.drm, cats = rep(2, 5), model = "DRM")
#'
#' @export
#'
shape_df <- function(par.drm = list(a = NULL, b = NULL, g = NULL),
                     par.prm = list(a = NULL, d = NULL),
                     item.id = NULL,
                     cats,
                     model,
                     default.par = FALSE) {

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
