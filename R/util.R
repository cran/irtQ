#' Bind Fill
#'
#' @description This function creates a cbind matrix or rbind matrix using a list containing different length
#' of numeric vectors.
#' @param List A list containing different length of numeric vectors
#' @param type A character string specifying whether rbind is used or cbind is used.
#' @param fill The value used to fill in missing data when aligning datasets. For \code{type = "cbind"},
#' this fills missing rows in shorter columns. For \code{type = "rbind"}, this fills missing columns in shorter rows.
#' Accepts any R object (e.g., numeric, character, logical). Defaults to NA.
#'
#' @return A matrix.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' # sample list
#' score_list <- list(item1=c(0:3), item2=c(0:2), item3=c(0:5), item3=c(0:4))
#'
#' # examples
#' # 1) create a rbind with the sample score list
#' bind.fill(score_list, type="rbind")
#'
#' # 2) create a cbind with the sample score list
#' bind.fill(score_list, type="cbind")
#'
#' # 3) create a cbind with the sample score list,
#' #    and fill missing data with 0s.
#' bind.fill(score_list, type="cbind", fill = 0L)
#'
#' @export
#' @import dplyr
bind.fill <- function(List, type=c("rbind", "cbind"), fill = NA){
  type <- tolower(type)
  type <- match.arg(type)
  nm <- List
  nm <- purrr::map(nm, as.matrix)
  names(nm) <- 1:length(nm)
  n <- max(purrr::map_dbl(nm, nrow))
  df <-
    purrr::map_dfc(nm, function(x) {rbind(x, matrix(fill, n - nrow(x), ncol(x)))}) %>%
    as.matrix()
  switch(type,
         cbind = unname(df),
         rbind = unname(t(df))
  )

}

# this function finds the index of the DRM items and PLM items given the object of
# "simdat()" function
idxfinder <- function(x) {

  # find the index of drm items
  idx.drm <- which(x$cats == 2)
  if(sum(idx.drm) == 0) idx.drm <- NULL

  # find the index of prm items
  idx.prm <- which(x$cats > 2)
  if(sum(idx.prm) == 0) idx.prm <- NULL

  # return the results
  list(idx.drm = idx.drm, idx.prm = idx.prm)

}

# a function to calculate a mean and variance at each theta point
cal_moment <- function(node, weight) {
  mu <- sum(node * weight)
  sigma2 <- sum(node^2 * weight) - mu^2
  rst <- c(mu=mu, sigma2=sigma2)
  rst
}


# This function divides the item response data sets into the two DRM responses (correct and incorrect)
# and one PRM item parts.
#' @importFrom Matrix Matrix
divide_data <- function(data, idx.item, freq.cat) {

  # divide the response data set into DRM and PRM parts
  if(!is.null(idx.item$idx.drm)) {
    data_drm_p <- Matrix::Matrix(data[, idx.item$idx.drm], sparse = TRUE)
    data_drm_q <- Matrix::Matrix(1 - data_drm_p, sparse = TRUE)
    data_drm_p[is.na(data_drm_p)] <- 0
    data_drm_q[is.na(data_drm_q)] <- 0
  } else {
    data_drm_p <- NULL
    data_drm_q <- NULL
  }
  if(!is.null(idx.item$idx.prm)) {
    data_prm <-
      Matrix::Matrix(do.call(what='cbind',
                             freq.cat[idx.item$idx.prm]), sparse = TRUE)
  } else {
    data_prm <- NULL
  }

  # create a response data matrix including all incorrect + correct responses
  data_all <- Matrix::Matrix(do.call(what = "cbind", freq.cat), sparse = TRUE)

  # return the results
  list(data_drm_p = data_drm_p, data_drm_q = data_drm_q,
       data_prm = data_prm, data_all = data_all)

}

# This function returns the column numbers of the frequency
# matrix corresponding to all items
cols4item <- function(nitem, cats, loc_1p_const=NULL) {

  cols.all <- vector('list', nitem)
  for(i in 1:nitem) {
    if(i == 1) {
      cat.st <- 1
      cat.ed <- cats[i]
    } else {
      cat.st <- cat.ed + 1
      cat.ed <- cat.st + (cats[i] - 1)
    }
    cols.tmp <- cat.st:cat.ed
    cols.all[[i]] <- cols.tmp
  }

  if(!is.null(loc_1p_const)) {
    cols.1pl <- unlist(cols.all[loc_1p_const])
  } else {
    cols.1pl <- NULL
  }

  # return
  rst <- list(cols.all = cols.all, cols.1pl=cols.1pl)
  rst

}


#' @export
coef.est_irt <- function(object, ...) {
  object$estimates
}

#' @export
coef.est_mg <- function(object, ...) {
  object$estimates
}

#' @export
coef.est_item <- function(object, ...) {
  object$estimates
}

#' @export
logLik.est_irt <- function(object, ...) {
  object$loglikelihood
}

#' @export
logLik.est_mg <- function(object, ...) {
  object$loglikelihood
}

#' @export
logLik.est_item <- function(object, ...) {
  object$loglikelihood
}

#' @export
vcov.est_irt <- function(object, ...) {
  object$covariance
}

#' @export
vcov.est_mg <- function(object, ...) {
  object$covariance
}

#' @export
vcov.est_item <- function(object, ...) {
  object$covariance
}

