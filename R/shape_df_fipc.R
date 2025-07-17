#' Combine fixed and new item metadata for fixed-item parameter calibration
#' (FIPC)
#'
#' This function merges existing fixed-item metadata with automatically
#' generated metadata for new items, producing a single data frame ordered by
#' specified test positions, to facilitate fixed item parameter calibration
#' using [irtQ::est_irt()].
#'
#' @param x A data.frame of metadata for items whose parameters remain fixed
#'   (e.g., output from [irtQ::shape_df()]).
#' @param fix.loc An integer vector specifying the row positions in the final
#'   output where fixed items should be placed.
#' @param item.id A character vector of IDs for new items whose parameters will
#'   be estimated.If `NULL`, default IDs (e.g., "V1", "V2", ...) are assigned
#'   automatically.
#' @param cats An integer vector indicating the number of response categories
#'   for each new item; order must match `item.id`.
#' @param model A character vector of IRT model names for each new item. Valid
#'   options for dichotomous items: "1PLM", "2PLM", "3PLM", "DRM"; for
#'   polytomous items: "GRM", "GPCM".
#'
#' @details To use this function, first prepare a metadata frame `x` containing
#'   only fixed items—either created by [irtQ::shape_df()] or imported from
#'   external software (e.g., via [irtQ::bring.flexmirt()]), which must include
#'   columns `id`, `cats`, `model`, and all relevant parameter columns (`par.1`,
#'   `par.2`, etc.). The `fix.loc` argument should then specify the exact row
#'   positions in the final test form where these fixed items should remain. The
#'   length of `fix.loc` must match the number of rows in `x`, and the order of
#'   positions in `fix.loc` determines where each fixed-item row is placed.
#'
#'   Next, provide information for the new items whose parameters will be
#'   estimated. Supply vectors for `item.id`, `cats`, and `model` matching the
#'   number of new items (equal to total form length minus length of `fix.loc`).
#'   If `item.id` is `NULL`, unique IDs are generated automatically.
#'
#'
#' @return A data.frame containing combined metadata for all items (fixed and
#'   new), ordered by test position.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::shape_df()]
#'
#' @examples
#' ## Import the flexMIRT parameter output file
#' prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#' x_fixed <- bring.flexmirt(file = prm_file, "par")$Group1$full_df
#'
#' ## Define positions of fixed items in the test form
#' fixed_pos <- c(1:40, 43:57)
#'
#' ## Specify IDs, models, and category counts for new items
#' new_ids <- paste0("NI", 1:6)
#' new_models <- c("3PLM", "1PLM", "2PLM", "GRM", "GRM", "GPCM")
#' new_cats <- c(2, 2, 2, 4, 5, 6)
#'
#' ## Generate combined metadata for FIPC
#' shape_df_fipc(x = x_fixed, fix.loc = fixed_pos, item.id = new_ids,
#'   cats = new_cats, model = new_models)
#'
#' @export
shape_df_fipc <- function(x, fix.loc = NULL, item.id = NULL, cats, model) {

  # Validate and standardize the fixed-item metadata
  x_fix <- confirm_df(x)

  # Generate default metadata for the new items
  x_new <- shape_df(item.id = item.id, cats = cats, model = model, default.par = TRUE)

  # Merge fixed and new metadata into a single data frame
  x_all <- dplyr::bind_rows(x_fix, x_new)

  # Determine the total number of items
  nitem <- nrow(x_all)

  # Identify row positions reserved for new items
  nfix.loc <- setdiff(seq_len(nitem), fix.loc)

  # Extract the original fixed-item rows
  x_fix2 <- x_all[seq_len(nrow(x_fix)), ]

  # Place fixed items and new items at their specified positions
  x_all[fix.loc, ] <- x_fix2
  x_all[nfix.loc, ] <- x_new

  # Return the combined and ordered metadata
  x_all

}

