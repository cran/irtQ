# "collapse_ftable_prm" function
# Vectorised wrapper around collapse_ftable() for the PRM category-collapsing
# step in sx2_fit().  Processes all summed-score groups for one polytomous item.
#
# Algorithm:
#   1. Rfast::rowmins() detects - in one vectorised pass - which score groups
#      require category collapsing (min expected freq < min.collapse).
#   2. Rows not needing collapse are extracted in bulk with asplit(), avoiding
#      per-row data.frame construction and function-call overhead.
#   3. Only rows that actually need collapsing invoke the original
#      collapse_ftable(), preserving bit-for-bit identical numerical output.
#
# exp_mat     : data.frame or matrix (nrows x K), expected frequencies
# obs_mat     : data.frame or matrix (nrows x K), observed frequencies
# min.collapse: minimum expected frequency per category (same as collapse_ftable)
# Returns: list(exp_table, obs_table), each a length-nrows list of numeric vectors
#' @importFrom Rfast rowMins
collapse_ftable_prm <- function(exp_mat, obs_mat, min.collapse = 1) {
  # convert to plain matrices for consistent row access and Rfast compatibility
  exp_m <- as.matrix(exp_mat)
  obs_m <- as.matrix(obs_mat)
  nrows <- nrow(exp_m)

  # vectorised scan: one pass to find the minimum expected freq per score group
  # value = TRUE returns the actual minimum value (default FALSE returns the index)
  row_min_exp <- Rfast::rowMins(exp_m, value = TRUE)
  needs_idx   <- which(row_min_exp < min.collapse)   # groups needing collapse
  pass_idx    <- which(row_min_exp >= min.collapse)  # groups with no collapse

  exp_table <- vector("list", nrows)
  obs_table <- vector("list", nrows)

  # fast path: bulk extraction for rows not needing any category collapsing;
  # asplit() splits the submatrix by row without per-row function-call overhead
  if (length(pass_idx) > 0L) {
    exp_table[pass_idx] <- lapply(asplit(exp_m[pass_idx, , drop = FALSE], 1L), unname)
    obs_table[pass_idx] <- lapply(asplit(obs_m[pass_idx, , drop = FALSE], 1L), unname)
  }

  # collapse path: rows with at least one category below threshold;
  # delegates to the original collapse_ftable() to guarantee identical output
  for (j in needs_idx) {
    x   <- data.frame(exp = exp_m[j, ], obs = obs_m[j, ])
    tmp <- collapse_ftable(x = x, col = 1L, min.collapse = min.collapse)
    exp_table[[j]] <- tmp$exp
    obs_table[[j]] <- tmp$obs
  }

  list(exp_table = exp_table, obs_table = obs_table)
}


# "collapse_ftable" function
# This function collapses the cells of a contingency table according to a column
collapse_ftable <- function(x, col, min.collapse = 1) {
  tmp <- x

  # check the locations of the cells in the selected column that have frequencies less than the minimum criterion
  loc_less <- which(tmp[, col] < min.collapse)

  # collapse cells
  while (length(loc_less) > 0) {
    # check the last row number in the contingency table
    last.num <- nrow(tmp)

    # check the center point
    center <- round(mean(1:last.num) + 0.01, 0)

    # relocate the location of rows that have frequencies less than minimum criterion
    # "loc_less_low" is the rows in which location is less than or equal to the center point
    loc_less_low <- loc_less[loc_less <= center]

    # "loc_less_high" is the rows in which location is greater than the center point
    loc_less_high <- rev(loc_less[loc_less > center])

    # relocation of the rows
    loc_less <- c(loc_less_low, loc_less_high)

    # check the location of the first selected cell
    start.num <- loc_less[1]

    # when the location of the first selected cell is 1
    if (start.num == 1) {
      row.1 <- tmp[start.num, ]
      row.2 <- tmp[start.num + 1, ]
      row.sum <- row.1 + row.2
      row.other <- tmp[-c(start.num, start.num + 1), ]
      tmp <- rbind(row.sum, row.other)
    }

    # when the location of the first selected cell is between 1 and last row number of the contingency table
    if (start.num >= 2 & start.num < last.num) {
      dif.num <- start.num - center
      sel <- ifelse(dif.num > 0, 1, 2)
      row.1 <- tmp[start.num, ]
      row.2 <- rbind(tmp[start.num - 1, ], tmp[start.num + 1, ])[sel, ]
      row.sum <- row.1 + row.2
      collapsed.nums <- sort(c(start.num, c(start.num - 1, start.num + 1)[sel]))
      row.other1 <- tmp[-c(collapsed.nums[1]:last.num), ]
      row.other2 <- tmp[-c(1:collapsed.nums[2]), ]
      tmp <- rbind(row.other1, row.sum, row.other2)
    }

    # when the location of the first selected cell is the last row number of the contingency table
    if (start.num == last.num) {
      row.1 <- tmp[start.num, ]
      row.2 <- tmp[start.num - 1, ]
      row.sum <- row.1 + row.2
      row.other <- tmp[-c(start.num, start.num - 1), ]
      tmp <- rbind(row.other, row.sum)
    }

    # after collapsing cells,
    # recheck the locations of the cells in the selected column that have frequencies
    # less than the minimum criterion
    loc_less <- which(round(tmp[, col], 15) < min.collapse)

    # count the number of remaining cells after collapsing
    # if the number of remaining cells is 1 and and the expected cell frequency is still less than
    # the specified value, then stop collapsing cell
    if (length(loc_less) == 1 & nrow(tmp) == 1) {
      break
    }
  }

  tmp
}
