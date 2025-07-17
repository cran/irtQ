#' Write a "-prm.txt" File for flexMIRT
#'
#' @description This function writes a flexMIRT-compatible "-prm.txt" file (Cai,
#'   2017). It currently supports only unidimensional IRT models. This function
#'   was developed by modifying `read.flexmirt()` from Pritikin & Falk (2020).
#'
#' @param x A data frame of item metadata (e.g., item parameters, number of
#'   categories, model types) for a single group, or a list of such data frames
#'   for multiple groups. See [irtQ::est_irt()] or [irtQ::simdat()] for item
#'   metadata format. You can also create metadata using [irtQ::shape_df()].
#' @param file A character string specifying the destination file path (with a
#'   ".txt" extension).
#' @param norm.pop A numeric vector of length two specifying the mean and
#'   standard deviation of the normal population ability distribution for a
#'   single group, or a list of such vectors for multiple groups. When a list is
#'   provided, each internal vector must contain the mean and standard deviation
#'   for a group's ability distribution (e.g., `norm.pop = list(c(0, 1), c(0,
#'   0.8), c(0.5, 1.2))` for three groups). If `mgroup = TRUE` and a single
#'   vector is provided (e.g., `norm.pop = c(0, 1)`), it will be recycled across
#'   all groups. The default is `c(0, 1)`.
#' @param rePar A logical value indicating whether the item parameters are
#'   reparameterized. If `TRUE`, item intercepts and logits of guessing
#'   parameters are assumed. If `FALSE`, item difficulty and guessing parameters
#'   are assumed.
#' @param mgroup A logical value indicating whether the file includes multiple
#'   groups. Default is `FALSE`.
#' @param group.name A character vector of group names. If `NULL`, group names
#'   are automatically generated (e.g., "Group1", "Group2", ...).
#'
#' @return This function creates a flexMIRT-style "-prm.txt" file at the
#'   specified path.
#'
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional
#'   item analysis and test scoring (Computer Software). Chapel Hill, NC: Vector
#'   Psychometric Group.
#'
#'   Pritikin, J. N., & Falk, C. F. (2020). OpenMx: A modular research
#'   environment for item response theory method development. *Applied
#'   Psychological Measurement, 44*(7-8), 561-562.
#'
#' @examples
#' \donttest{
#' ## 1. Create a "-prm.txt" file for a single group
#' ##    using the simulated CAT data
#' # 1-(1) Extract the item metadata
#' x <- simCAT_MX$item.prm
#'
#' # 1-(2) Set the name of the "-prm.txt" file
#' temp_prm <- file.path(tempdir(), "single_group_temp-prm.txt")
#'
#' # 1-(3) Write the "-prm.txt" file
#' write.flexmirt(x, file = temp_prm, norm.pop = c(0, 1), rePar = FALSE)
#'
#' ## 2. Create a "-prm.txt" file for multiple groups
#' ##    using simulated multi-group data
#' # 2-(1) Extract the item metadata
#' x <- simMG$item.prm
#'
#' # Set the name of the "-prm.txt" file
#' temp_prm <- file.path(tempdir(), "mg_group_temp-prm1.txt")
#'
#' # Write the "-prm.txt" file
#' write.flexmirt(x,
#'   file = temp_prm, norm.pop = list(c(0, 1), c(0.5, 0.8), c(-0.3, 1.3)),
#'   rePar = FALSE, mgroup = TRUE, group.name = c("GR1", "GR2", "GR3")
#' )
#'
#' # Or write the "-prm.txt" file so that
#' # all groups share the same ability distribution
#' # and group names are generated automatically
#' temp_prm <- file.path(tempdir(), "mg_group_temp-prm2.txt")
#' write.flexmirt(x,
#'   file = temp_prm, norm.pop = c(0, 1),
#'   rePar = FALSE, mgroup = TRUE, group.name = NULL
#' )
#' }
#'
#' @export
write.flexmirt <- function(x,
                           file = NULL,
                           norm.pop = c(0, 1),
                           rePar = TRUE,
                           mgroup = FALSE,
                           group.name = NULL) {

  if (!mgroup) {
    ngroup <- 1
    x_gr <- list(x)
    norm.pop <- list(norm.pop)
    if (is.null(group.name)) {
      group.name <- "Group1"
    }
  } else {
    x_gr <- x
    ngroup <- length(x)
    if (!is.list(norm.pop)) {
      norm.pop <- purrr::map(.x = 1:ngroup, ~ {
        norm.pop
      })
    }
    if (is.null(group.name)) {
      group.name <- paste0("Group", 1:ngroup)
    }
  }

  # open a new file
  prm_file <- file(file, open = "w")

  for (g in 1:ngroup) {
    # give column names
    x <- data.frame(x_gr[[g]])
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

    # warning message
    if (is.null(file)) {
      stop("You must specify the destination of the file.")
    }

    # select parameter columns
    param <- dplyr::select(x, dplyr::starts_with("par."))

    # write prm file one item by one item
    for (i in 1:nrow(x)) {
      # check metadata for each item
      id <- as.character(x$id[i])
      cats <- x$cats[i]
      model <- as.character(x$model[i])
      group <- g
      n.factor <- 1

      # (1) Dichotomous item
      if (model %in% c("3PLM", "DRM")) {
        a <- param[i, 1]
        if (!rePar) {
          c <- -a * param[i, 2]
          if (param[i, 3] == 0L) {
            param[i, 3] <- 1e-100
          }
          logitg <- log(param[i, 3] / (1 - param[i, 3]))
        } else {
          c <- param[i, 2]
          logitg <- param[i, 3]
        }
        a <- format(a, nsmall = 7)
        c <- format(c, nsmall = 7)
        logitg <- format(logitg, nsmall = 7)
        cat(c(1, id, group, n.factor, 1, cats, logitg, c, a),
          sep = "\t",
          file = prm_file, fill = TRUE
        )
      }

      if (model %in% c("1PLM", "2PLM")) {
        a <- param[i, 1]
        if (!rePar) {
          c <- -a * param[i, 2]
        } else {
          c <- param[i, 2]
        }
        a <- format(a, nsmall = 7)
        c <- format(c, nsmall = 7)
        cat(c(1, id, group, n.factor, 2, cats, c, a),
          sep = "\t",
          file = prm_file, fill = TRUE
        )
      }

      # (2) Polytomous item: GRM
      if (model %in% c("GRM")) {
        a <- param[i, 1]
        if (!rePar) {
          cs <- -a * as.numeric(param[i, 2:cats])
        } else {
          cs <- as.numeric(param[i, 2:cats])
        }
        a <- format(a, nsmall = 7)
        cs <- format(cs, nsmall = 7)
        cat(c(1, id, group, n.factor, 2, cats, cs, a),
          sep = "\t",
          file = prm_file, fill = TRUE
        )
      }

      # (3) Polytomous item: PCM and GPCM
      if (model %in% c("PCM", "GPCM")) {
        a <- param[i, 1]
        if (!rePar) {
          alpha <- c(1, rep(0, cats - 2))
          ds.new <- param[i, 2:cats]
          b <- sum(ds.new) / 4
          ds <- b - ds.new

          Tmat <- matrix(0, nrow = cats, ncol = (cats - 1))
          Tmat[, 1] <- 0:(cats - 1)
          for (k in 2:(cats - 1)) {
            for (j in 2:(cats - 1)) {
              Tmat[k, j] <- sin(pi * (j - 1) * (k - 1) / (cats - 1))
            }
          }

          c.vec <- rep(0, cats)
          c.vec[cats] <- -a * (cats - 1) * b
          for (k in length(ds):2) {
            c.vec[k] <- c.vec[k + 1] - a * (as.numeric(ds[k]) - b)
          }

          gam <- as.numeric(solve(Tmat[-1, ]) %*% c.vec[-1])
        } else {
          alpha <- as.numeric(param[i, 2:(cats - 1)])
          gam <- as.numeric(param[i, (1 + cats - 2 + 1 + 1):(1 + cats - 2 + 1 + 1 + cats - 2)])
        }
        a <- format(a, nsmall = 7)
        alpha <- format(alpha, nsmall = 7)
        gam <- format(gam, nsmall = 7)
        cat(c(1, id, group, n.factor, 3, cats, 0, alpha, a, 0, gam),
          sep = "\t",
          file = prm_file, fill = TRUE
        )
      }
    }

    # group information
    cat(c(0, group.name[g], g, n.factor, 0, format(norm.pop[[g]][1], nsmall = 7), format(norm.pop[[g]][2], nsmall = 7)),
      sep = "\t", file = prm_file, fill = TRUE
    )
  }

  # close the file
  close(prm_file)
}
