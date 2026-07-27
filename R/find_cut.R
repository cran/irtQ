# find_cut.R
# ---------------------------------------------------------------------------
# find_cut() -- Identify TIF-crossing cut scores for MST routing
#
# For each pair of adjacent modules (in the order provided by the user,
# i.e., ascending module index) at each stage transition, this function finds
# the theta value where the right module's TIF curve crosses the left module's
# TIF curve with a positive slope (a "proper" crossing).  The resulting cut
# scores can be passed directly to run_mst() as the cut_score argument.
# ---------------------------------------------------------------------------


#' Find TIF-Crossing Cut Scores for MST Routing
#'
#' @description
#' Computes routing cut scores for a multistage test (MST) by identifying the
#' theta values at which adjacent modules' Test Information Functions (TIFs)
#' cross in the correct direction. The resulting cut scores can be passed
#' directly to \code{\link{run_mst}} as its \code{cut_score} argument
#' (with \code{route_method = NULL}), providing a principled and psychometrically
#' grounded alternative to heuristic cut-score selection.
#'
#' @param x A data frame of item metadata in the standard \pkg{irtQ} format
#'   (columns: \code{id}, \code{cats}, \code{model}, \code{par.1},
#'   \code{par.2}, \ldots). See \code{\link{shape_df}} for details.
#' @param module A binary integer matrix of dimensions \emph{J} \eqn{\times}
#'   \emph{M}, where \emph{J} is the number of items and \emph{M} is the
#'   number of modules. A value of 1 at row \emph{j}, column \emph{m}
#'   indicates that item \emph{j} belongs to module \emph{m}. This is the
#'   same \code{module} argument used in \code{\link{run_mst}}.
#' @param route_map A binary square matrix of dimensions \emph{M} \eqn{\times}
#'   \emph{M} defining the MST transition structure. A value of 1 at row
#'   \emph{i}, column \emph{j} indicates that examinees can be routed from
#'   module \emph{i} to module \emph{j}. This is the same \code{route_map}
#'   argument used in \code{\link{run_mst}}.
#' @param D A numeric scaling constant for the IRT model. Default is
#'   \code{1.702}.
#' @param theta_range A numeric vector of length 2 specifying the theta range
#'   over which to search for TIF crossings. Default is \code{c(-6, 6)}.
#' @param n_grid An integer specifying the number of equally spaced theta
#'   points used for the initial sign-change scan. A finer grid reduces the
#'   chance of missing a crossing but increases computation time. Default is
#'   \code{2001}.
#' @param ref_theta A single numeric value used to break ties when multiple
#'   proper crossings are found for an adjacent module pair. The crossing
#'   closest to \code{ref_theta} is selected. Default is \code{0}.
#'
#' @return An object of class \code{"find_cut"}, which is a named list with
#'   three elements:
#' \describe{
#'   \item{\code{cut_score}}{A list of length \emph{S} - 1 (where \emph{S}
#'     is the number of stages), compatible with the \code{cut_score} argument
#'     of \code{\link{run_mst}}. Each element is a numeric vector of cut
#'     scores for the transition into that stage, sorted in ascending order.
#'     Stages with only one module have an empty numeric vector
#'     (\code{numeric(0)}).}
#'   \item{\code{details}}{A named list of per-stage diagnostic information,
#'     including module indices, mean item locations, difficulty-sorted module
#'     order, and per-pair crossing details (proper crossings, anomalous
#'     crossings, and the selected cut score).}
#'   \item{\code{tif_data}}{A \code{tibble} with columns \code{theta},
#'     \code{module}, \code{stage}, and \code{tif}, providing TIF curves for
#'     all modules at all stages (including stage 1). Used by
#'     \code{\link{plot.find_cut}} to visualise TIF curves and crossing points.}
#' }
#'
#' @details
#' \subsection{Motivation: the MFI routing problem}{
#' When \code{\link{run_mst}} uses Maximum Fisher Information (MFI) routing
#' (\code{route_method = "mfi"}), it selects the next-stage module that
#' maximizes TIF at the examinee's current ability estimate. In principle, a
#' harder module should always have a higher TIF for high-ability examinees,
#' and a lower TIF for low-ability examinees, making MFI routing work as
#' intended. However, in practice, a harder module's TIF curve can peak at
#' a theta value that is unexpectedly low, so that the harder module
#' \emph{also} has the highest TIF near the low end of the ability scale.
#' In this case, a low-ability examinee is incorrectly routed into the harder
#' module -- the opposite of what the MST is designed to do, which is called an
#' \emph{anomalous routing} or \emph{path reversal}.
#' }
#'
#' \subsection{Solution: TIF-crossing cut scores}{
#' \code{find_cut()} resolves this problem by computing cut scores from the
#' TIF curves themselves, rather than relying on pure MFI at the time of
#' routing. Specifically, for each pair of adjacent modules (ordered by mean
#' item difficulty), the function finds the theta value at which the harder
#' module's TIF first overtakes the easier module's TIF \emph{from below}
#' as theta increases. This crossing point is the natural boundary between
#' the two modules: below it, the easier module is more informative; above it,
#' the harder module is more informative. Using this point as a fixed cut
#' score pre-empts the path reversal problem entirely.
#' }
#'
#' \subsection{Proper vs. anomalous crossings}{
#' Two TIF curves may cross more than once. \code{find_cut()} distinguishes
#' two types of crossing:
#' \itemize{
#'   \item \strong{Proper crossing}: the difference
#'     \eqn{\text{TIF}_{\text{harder}}(\theta) - \text{TIF}_{\text{easier}}(\theta)}
#'     changes sign from negative to positive (positive slope at the crossing).
#'     This is the correct, expected crossing in the middle of the theta scale.
#'   \item \strong{Anomalous crossing}: the difference changes sign from
#'     positive to negative (negative slope at the crossing). This corresponds
#'     to the pathological situation in which the harder module temporarily
#'     dominates at the low end of the theta scale -- exactly the pattern that
#'     causes path reversals under MFI routing.
#' }
#' Only proper crossings are used as cut scores. Anomalous crossings are
#' reported in \code{$details} for diagnostic purposes but excluded from the
#' returned \code{cut_score}. If no proper crossing is found for a module pair,
#' an informative error is thrown. If multiple proper crossings are found
#' (unusual for well-designed modules), a warning is issued and the one
#' closest to \code{ref_theta} is selected.
#' }
#'
#' \subsection{Algorithm}{
#' The function proceeds stage by stage (from stage 2 onward). At each stage:
#' \enumerate{
#'   \item For each adjacent pair of modules in the order supplied by the
#'     user (ascending module index), the difference
#'     \eqn{\text{TIF}_{\text{right}}(\theta) -
#'     \text{TIF}_{\text{left}}(\theta)} is evaluated on a fine theta grid
#'     of \code{n_grid} points, where "left" and "right" refer to the lower-
#'     and higher-index module in the pair. Sign changes on the grid are
#'     detected and then refined with \code{\link[stats]{uniroot}}.
#'   \item The slope at each refined root is estimated by a finite difference
#'     (step size 0.01) to classify it as proper or anomalous.
#'   \item The selected proper crossing is stored as the cut score for that pair.
#' }
#' }
#'
#' \subsection{Module index vs. difficulty order}{
#' \code{\link{run_mst}} internally orders candidate modules by ascending
#' module index and maps routing rank 1 to the lowest-index module. For
#' cut-score routing to assign low-theta examinees to the easiest module, the
#' modules at each stage must be numbered in ascending difficulty order
#' (easiest module = lowest index). \code{find_cut()} checks this assumption and issues a warning if the
#' difficulty order does not match the index order. Regardless of this warning,
#' cut scores are always computed between consecutive module pairs in the
#' supplied index order (e.g., modules 2 & 3 and 3 & 4, not 2 & 4 and 4 & 3).
#' The returned cut scores within each stage are sorted in ascending order,
#' as required by \code{\link[base]{cut}}.
#' }
#'
#' @seealso \code{\link{run_mst}}, \code{\link{panel_info}},
#'   \code{\link{reval_mst}}, \code{\link{info}}
#'
#' @importFrom stats setNames
#'
#' @examples
#' ## -- Setup: use the built-in simMST 1-3-3 panel --------------------------
#' ## simMST is a 7-module, 3-stage MST panel (8 dichotomous 3PLM items each).
#' ## Modules 1 (routing), 2-4 (stage 2: easy/medium/hard),
#' ##         5-7 (stage 3: easy/medium/hard).
#' x         <- simMST$item_bank
#' module    <- simMST$module
#' route_map <- simMST$route_map
#'
#' ## -- Find TIF-crossing cut scores -----------------------------------------
#' ## For each adjacent module pair at stages 2 and 3, find_cut() identifies
#' ## the theta at which the harder module's TIF first exceeds the easier
#' ## module's TIF (proper crossing), and returns it as a cut score.
#' cut_result <- find_cut(x = x, module = module, route_map = route_map)
#'
#' ## Print a summary: crossing points found, anomalous crossings excluded,
#' ## and the final selected cut scores per stage transition.
#' print(cut_result)
#'
#' ## -- Visualise TIF curves and cut scores ---------------------------------
#' ## plot() shows TIF curves for every module faceted by stage.
#' ## Selected cut scores appear as solid black vertical lines with theta
#' ## labels at the crossing point. Anomalous crossings (if any) appear
#' ## as dashed red lines.
#' plot(cut_result)                         # stages stacked vertically (default)
#' plot(cut_result, layout = "horizontal")  # stages side by side
#' plot(cut_result, show_anomalous = FALSE) # hide anomalous crossing markers
#'
#' ## Inspect the cut_score element -- a list directly compatible with run_mst()
#' ## cut_score[[1]]: two cut scores for the stage-1 -> stage-2 transition
#' ## cut_score[[2]]: cut scores for the stage-2 -> stage-3 transition
#' cut_result$cut_score
#'
#' ## Compare with the manually specified cut scores stored in simMST
#' simMST$cut_score
#'
#' ## -- Use the cut scores in run_mst() --------------------------------------
#' ## Pass cut_result$cut_score directly to run_mst() with route_method = NULL.
#' ## This replaces pure MFI routing with TIF-crossing-based fixed cut scores,
#' ## preventing anomalous path reversals at the extremes of the theta scale.
#' \dontrun{
#' set.seed(1)
#' theta_true <- rnorm(500)
#' result <- run_mst(
#'   x            = x,
#'   route_map    = route_map,
#'   module       = module,
#'   theta        = theta_true,
#'   route_method = NULL,
#'   cut_score    = cut_result$cut_score
#' )
#' }
#'
#' @export
find_cut <- function(x,
                     module,
                     route_map,
                     D           = 1.702,
                     theta_range = c(-6, 6),
                     n_grid      = 2001L,
                     ref_theta   = 0) {

  # -- 0. Input validation ----------------------------------------------------
  if (!is.data.frame(x))
    stop("'x' must be a data frame of item metadata in irtQ format.")
  if (!is.matrix(module) || !is.numeric(module))
    stop("'module' must be a numeric binary matrix (J x M).")
  if (nrow(x) != nrow(module))
    stop("Number of rows in 'x' must equal number of rows in 'module' (J items).")
  # route_map may be a matrix or a data.frame (panel_info() accepts both)
  if (!(is.matrix(route_map) || is.data.frame(route_map)) ||
      !all(vapply(as.data.frame(route_map), is.numeric, logical(1L))))
    stop("'route_map' must be a numeric binary square matrix (M x M).")
  if (ncol(module) != nrow(route_map))
    stop("Number of columns in 'module' (M modules) must equal dimensions of 'route_map'.")
  if (!is.numeric(D) || length(D) != 1L || !is.finite(D))
    stop("'D' must be a single finite numeric scalar.")
  if (!is.numeric(theta_range) || length(theta_range) != 2L ||
      theta_range[1L] >= theta_range[2L])
    stop("'theta_range' must be a numeric vector of length 2 with theta_range[1] < theta_range[2].")
  n_grid <- as.integer(n_grid)                    # coerce to integer
  if (n_grid < 100L)
    stop("'n_grid' must be an integer >= 100.")
  if (!is.numeric(ref_theta) || length(ref_theta) != 1L || !is.finite(ref_theta))
    stop("'ref_theta' must be a single finite numeric scalar.")

  # -- 1. Panel structure -----------------------------------------------------
  # Use panel_info() to identify which modules belong to each stage
  pinfo  <- panel_info(route_map)
  config <- pinfo$config    # named list: config[[s]] = integer vector of module indices
  n_stg  <- pinfo$n.stage   # total number of stages (integer)
  n_mods <- ncol(module)    # total number of modules M

  # Build a list of item metadata subsets, one per module
  # item_mod[[m]]: rows of x assigned to module m (module[,m] == 1)
  item_mod <- lapply(seq_len(n_mods), function(m) {
    x[module[, m] == 1L, , drop = FALSE]
  })

  # Fine theta grid for sign-change scan (root detection)
  theta_grid <- seq(theta_range[1L], theta_range[2L], length.out = n_grid)

  # -- 2. Allocate output containers ------------------------------------------
  # cut_list[[k]]: cut score vector for transition from stage k to stage k+1
  # Indexed 1 ... n_stg-1, matching the cut_score argument of run_mst()
  cut_list     <- vector("list", n_stg - 1L)
  names(cut_list) <- paste0("stage.", seq(2L, n_stg))

  # details[[k]]: per-stage diagnostic information
  details      <- vector("list", n_stg - 1L)
  names(details) <- names(cut_list)

  # tif_rows: accumulator for tif_data tibble rows
  tif_rows     <- list()

  # -- 2b. Collect TIF data for stage 1 (routing stage) ----------------------
  # Stage 1 is not processed in the main loop (no routing decision needed),
  # but its TIF is included in tif_data so plot.find_cut() can display it
  # as the top facet for context.
  for (m in config[[1L]]) {
    tif_vals <- info(x = item_mod[[m]], theta = theta_grid,
                     D = D, tif = TRUE)$tif
    tif_rows[[length(tif_rows) + 1L]] <- tibble::tibble(
      theta  = theta_grid,
      module = m,
      stage  = 1L,
      tif    = tif_vals
    )
  }

  # -- 3. Main loop: process each stage from stage 2 onward ------------------
  for (s in seq(2L, n_stg)) {

    mod_ids <- config[[s]]    # module indices at stage s (integer vector)
    n_mod_s <- length(mod_ids) # number of modules at this stage

    # -- 3a. Single module: no routing decision at this stage -----------------
    if (n_mod_s == 1L) {
      cut_list[[s - 1L]] <- numeric(0L)   # empty: no cut score needed
      details[[s - 1L]]  <- list(
        modules      = mod_ids,
        mean_locs    = NULL,
        sorted_order = NULL,
        pairs        = NULL,
        message      = "Single module at this stage; no cut score required."
      )
      # Still collect TIF data for the single module (for tif_data tibble)
      tif_vals <- info(x = item_mod[[mod_ids]], theta = theta_grid,
                       D = D, tif = TRUE)$tif
      tif_rows[[length(tif_rows) + 1L]] <- tibble::tibble(
        theta  = theta_grid,
        module = mod_ids,
        stage  = s,
        tif    = tif_vals
      )
      next   # move to next stage
    }

    # -- 3b. Precompute TIF matrix for all modules at stage s -----------------
    # tif_mat: n_grid x n_mod_s matrix; column k = TIF of mod_ids[k]
    tif_mat <- vapply(
      X         = mod_ids,
      FUN       = function(m) info(x = item_mod[[m]], theta = theta_grid,
                                   D = D, tif = TRUE)$tif,
      FUN.VALUE = numeric(n_grid)
    )

    # Accumulate TIF rows for visualization tibble
    for (ki in seq_len(n_mod_s)) {
      tif_rows[[length(tif_rows) + 1L]] <- tibble::tibble(
        theta  = theta_grid,
        module = mod_ids[ki],   # module index
        stage  = s,
        tif    = tif_mat[, ki]  # TIF values for this module
      )
    }

    # -- 3c. Sort modules by ascending mean item location (difficulty) ---------
    # mean_loc() returns a vector of per-item mean location parameters
    mean_locs <- vapply(
      X         = mod_ids,
      FUN       = function(m) mean(mean_loc(item_mod[[m]])),
      FUN.VALUE = numeric(1L)
    )

    diff_ord    <- order(mean_locs)       # ascending difficulty -> index in mod_ids
    sorted_mods <- mod_ids[diff_ord]      # module ids sorted by ascending difficulty
    sorted_locs <- mean_locs[diff_ord]    # corresponding mean locations

    # -- 3d. Check: index order should match difficulty order ------------------
    # run_mst() maps routing rank 1 -> next_possible[1] (lowest module index).
    # For cut-score routing to assign low theta -> easy module, the modules at
    # each stage must have ascending indices in ascending difficulty order.
    if (!identical(sorted_mods, sort(mod_ids))) {
      warning(sprintf(paste0(
        "Stage %d: difficulty order of modules does not match their index order.\n",
        "  Ascending index order    : modules %s\n",
        "  Ascending difficulty order: modules %s (mean locations: %s)\n",
        "  run_mst() routes by ascending module index (rank 1 -> lowest index).\n",
        "  If index order differs from difficulty order, cut-score routing may\n",
        "  not produce the intended module assignments.\n",
        "  Consider redesigning the panel so that easier modules have lower indices."
      ),
      s,
      paste(sort(mod_ids),         collapse = ", "),
      paste(sorted_mods,           collapse = ", "),
      paste(round(sorted_locs, 3L), collapse = ", ")
      ))
    }

    # -- 3e. Find cut score for each adjacent pair in module index order --------
    # IMPORTANT: pairs always follow the input module index order (e.g., 2<->3,
    # 3<->4), NOT the difficulty-sorted order.  This matches how run_mst() maps
    # routing rank to module index via next_possible (ascending index sort).
    n_pairs   <- n_mod_s - 1L                  # number of adjacent pairs
    pair_cuts <- numeric(n_pairs)              # one cut score per pair
    pair_info <- vector("list", n_pairs)       # diagnostic info per pair

    for (k in seq_len(n_pairs)) {

      left_m  <- mod_ids[k]          # left (lower-index) module
      right_m <- mod_ids[k + 1L]    # right (higher-index) module

      # Column positions in tif_mat match the order of mod_ids directly
      col_l    <- k          # left module is the k-th column in tif_mat
      col_r    <- k + 1L    # right module is the (k+1)-th column

      # Difference of TIF values on the grid: positive means right > left
      diff_grid <- tif_mat[, col_r] - tif_mat[, col_l]

      # Inline evaluation function for uniroot() and derivative estimation
      # Captures left_m and right_m from the current loop iteration
      l_m <- left_m    # local copy to avoid closure pitfalls
      r_m <- right_m
      diff_tif_fn <- function(theta_val) {
        # Evaluate TIF difference (right minus left) at an arbitrary theta
        tif_r <- info(x = item_mod[[r_m]], theta = theta_val, D = D, tif = TRUE)$tif
        tif_l <- info(x = item_mod[[l_m]], theta = theta_val, D = D, tif = TRUE)$tif
        as.numeric(tif_r - tif_l)
      }

      # Detect sign-change positions on the fine grid
      # diff() of sign() gives +2, -2, or 0; nonzero means a sign change occurred
      chng_idx <- which(diff(sign(diff_grid)) != 0L)

      # Error: no sign changes -> no crossing in the search range
      if (length(chng_idx) == 0L) {
        stop(sprintf(paste0(
          "Stage %d, pair (module %d vs module %d): ",
          "no TIF crossing detected in theta range [%.2f, %.2f].\n",
          "The TIF curves do not intersect within this range.\n",
          "Consider widening 'theta_range' or reviewing module item composition."
        ), s, left_m, right_m, theta_range[1L], theta_range[2L]))
      }

      # Refine each sign-change interval with uniroot() and classify the root
      all_roots  <- numeric(0L)      # all refined roots (proper + anomalous)
      root_types <- character(0L)    # "proper" or "anomalous" for each root

      for (ci in chng_idx) {
        lo <- theta_grid[ci]        # left endpoint of sign-change interval
        hi <- theta_grid[ci + 1L]  # right endpoint

        # Re-evaluate at endpoints (guard against numerical ties on grid)
        val_lo <- diff_tif_fn(lo)
        val_hi <- diff_tif_fn(hi)
        if (sign(val_lo) == sign(val_hi)) next   # not a genuine crossing

        # Refine the root position with uniroot()
        root_i <- tryCatch(
          stats::uniroot(
            f        = diff_tif_fn,
            interval = c(lo, hi),
            tol      = .Machine$double.eps^0.5   # tight tolerance
          )$root,
          error = function(e) NA_real_
        )
        if (is.na(root_i)) next    # uniroot failed -- skip this candidate

        # Classify the crossing by the sign of the slope (finite difference, h = 0.01)
        slope_i <- (diff_tif_fn(root_i + 0.01) - diff_tif_fn(root_i - 0.01)) / 0.02
        type_i  <- if (slope_i > 0) "proper" else "anomalous"

        all_roots  <- c(all_roots,  root_i)
        root_types <- c(root_types, type_i)
      }

      proper_roots    <- all_roots[root_types == "proper"]    # valid cut candidates
      anomalous_roots <- all_roots[root_types == "anomalous"] # pathological crossings

      # -- Error: no proper crossings ------------------------------------------
      if (length(proper_roots) == 0L) {
        if (length(anomalous_roots) > 0L) {
          # Crossings exist but ALL are anomalous (right module dominates at low theta)
          stop(sprintf(paste0(
            "Stage %d, pair (module %d vs module %d): ",
            "no valid cut score found.\n",
            "All %d crossing(s) detected have a negative slope (anomalous pattern):\n",
            "  theta = [%s]\n",
            "This indicates that the right (higher-index) module has higher TIF at\n",
            "low theta values, which may cause misrouting of low-ability examinees.\n",
            "Review module item composition or consult a test design specialist."
          ),
          s, left_m, right_m,
          length(anomalous_roots),
          paste(round(anomalous_roots, 4L), collapse = ", ")))
        } else {
          # Sign changes existed on grid but uniroot() found nothing -- numerical issue
          stop(sprintf(paste0(
            "Stage %d, pair (module %d vs module %d): ",
            "sign changes were detected on the grid but no root could be refined.\n",
            "Try increasing 'n_grid' or widening 'theta_range'."
          ), s, left_m, right_m))
        }
      }

      # -- Warning: multiple proper crossings -> select closest to ref_theta ---
      if (length(proper_roots) > 1L) {
        selected_root <- proper_roots[which.min(abs(proper_roots - ref_theta))]
        warning(sprintf(paste0(
          "Stage %d, pair (module %d vs module %d): ",
          "%d proper crossings found at theta = [%s].\n",
          "Selecting the crossing closest to ref_theta = %.4f: theta = %.4f.\n",
          "To override the selection, specify a different 'ref_theta'."
        ),
        s, left_m, right_m,
        length(proper_roots),
        paste(round(proper_roots, 4L), collapse = ", "),
        ref_theta, selected_root))
      } else {
        selected_root <- proper_roots[1L]    # exactly one proper crossing
      }

      pair_cuts[k] <- selected_root    # store selected cut score for this pair

      # Record per-pair diagnostic information
      pair_info[[k]] <- list(
        left_module         = left_m,            # lower-index module (left in pair)
        right_module        = right_m,           # higher-index module (right in pair)
        proper_crossings    = proper_roots,      # all valid crossing points
        anomalous_crossings = anomalous_roots,   # excluded pathological crossings
        selected_cut        = selected_root      # the cut score used
      )
      names(pair_info)[k] <- sprintf("mod%d_vs_mod%d", left_m, right_m)

    }  # end pair loop (k)

    # -- 3f. Store stage-level results -----------------------------------------
    # Cut scores must be sorted ascending: give_path() uses them as breaks in
    # cut(x, breaks = c(-Inf, cut_sc, Inf)), so ascending order is required.
    cut_list[[s - 1L]] <- sort(pair_cuts)

    # Collect diagnostic information for this stage
    details[[s - 1L]] <- list(
      modules      = mod_ids,                              # module indices at this stage
      mean_locs    = setNames(mean_locs, paste0("mod", mod_ids)),  # named mean locations
      sorted_order = sorted_mods,                         # difficulty-sorted module order
      pairs        = pair_info                            # per-pair crossing details
    )

  }  # end stage loop (s)

  # -- 4. Build tif_data tibble -----------------------------------------------
  # Tidy data frame with one row per (theta, module) combination
  # Suitable for ggplot2 visualisation of TIF curves and crossing points
  tif_data <- dplyr::bind_rows(tif_rows)

  # -- 5. Assemble and return the result object -------------------------------
  rst <- list(
    cut_score = cut_list,    # list: pass directly to run_mst(cut_score = ...)
    details   = details,     # list: per-stage diagnostic info
    tif_data  = tif_data     # tibble: TIF curves for all stages (1 and above)
  )
  class(rst) <- "find_cut"   # S3 class for print dispatch
  rst
}


# ---------------------------------------------------------------------------
# plot.find_cut() -- S3 plot method for find_cut objects
# ---------------------------------------------------------------------------

#' Plot TIF Curves and Cut Scores from a \code{find_cut} Result
#'
#' @description
#' Produces a \pkg{ggplot2} visualisation of the Test Information Function
#' (TIF) curves for each module, faceted by stage, with vertical lines marking
#' the selected cut scores (proper TIF crossings) and, optionally, any
#' anomalous crossings detected by \code{\link{find_cut}}.
#'
#' Each stage occupies one panel. Stage 1 (routing module) appears first.
#' Stages with routing decisions show TIF curves for all candidate modules and
#' mark the selected cut score(s) with a solid black vertical line and a theta
#' value label at the crossing point.
#'
#' @param x An object of class \code{"find_cut"} returned by
#'   \code{\link{find_cut}}.
#' @param theta_range A numeric vector of length 2 specifying the theta range
#'   shown on the x-axis. Defaults to the full range stored in
#'   \code{x$tif_data}.
#' @param show_anomalous Logical. If \code{TRUE} (default), anomalous TIF
#'   crossings (negative-slope crossings excluded from the cut scores) are
#'   shown as dashed red vertical lines.
#' @param label_cuts Logical. If \code{TRUE} (default), each selected cut
#'   score is annotated with its theta value above the crossing point.
#' @param layout Character string controlling the panel arrangement.
#'   \code{"vertical"} (default) stacks each stage in its own row (one
#'   column of panels, stage 1 at the top). \code{"horizontal"} places each
#'   stage in its own column (one row of panels, stage 1 at the left).
#' @param ... Currently unused. Reserved for future arguments.
#'
#' @return A \code{ggplot} object. The object is printed automatically when
#'   called interactively. It can be further customised with standard
#'   \pkg{ggplot2} functions such as \code{ggplot2::theme()},
#'   \code{ggplot2::scale_color_manual()}, etc.
#'
#' @seealso \code{\link{find_cut}}, \code{\link{run_mst}}
#'
#' @importFrom rlang .data
#'
#' @examples
#' ## Use the built-in simMST 1-3-3 panel
#' cut_result <- find_cut(
#'   x         = simMST$item_bank,
#'   module    = simMST$module,
#'   route_map = simMST$route_map
#' )
#'
#' ## Default: stages stacked vertically (stage 1 at the top)
#' plot(cut_result)
#'
#' ## Horizontal layout: stages side by side
#' plot(cut_result, layout = "horizontal")
#'
#' ## Hide anomalous crossing markers
#' plot(cut_result, show_anomalous = FALSE)
#'
#' @export
plot.find_cut <- function(x,
                          theta_range    = NULL,
                          show_anomalous = TRUE,
                          label_cuts     = TRUE,
                          layout         = c("vertical", "horizontal"),
                          ...) {

  # -- 0. Argument validation -------------------------------------------------
  if (!inherits(x, "find_cut"))
    stop("'x' must be an object of class \"find_cut\".")
  layout <- match.arg(layout)    # "vertical" or "horizontal"
  if (!is.logical(show_anomalous) || length(show_anomalous) != 1L)
    stop("'show_anomalous' must be a single logical value (TRUE or FALSE).")
  if (!is.logical(label_cuts) || length(label_cuts) != 1L)
    stop("'label_cuts' must be a single logical value (TRUE or FALSE).")

  # -- 1. Prepare TIF tibble for plotting ------------------------------------
  tif_plot <- x$tif_data    # columns: theta, module, stage, tif

  # Apply optional theta range filter
  if (!is.null(theta_range)) {
    if (!is.numeric(theta_range) || length(theta_range) != 2L ||
        theta_range[1L] >= theta_range[2L])
      stop("'theta_range' must be a numeric vector of length 2 with ",
           "theta_range[1] < theta_range[2].")
    tif_plot <- tif_plot[tif_plot$theta >= theta_range[1L] &
                         tif_plot$theta <= theta_range[2L], ]
  }

  # Create ordered stage factor (stage 1 first = top / leftmost panel)
  all_stages   <- sort(unique(tif_plot$stage))
  stage_levels <- paste0("Stage ", all_stages)

  tif_plot$stage_label  <- factor(paste0("Stage ", tif_plot$stage),
                                  levels = stage_levels)
  tif_plot$module_label <- paste0("Module ", tif_plot$module)

  # -- 2. Build annotation data for vertical lines and crossing dots ----------
  vline_rows <- list()   # rows for vline data frame
  point_rows <- list()   # rows for crossing-point dot data frame

  has_unselected <- FALSE   # any unselected proper crossings present?
  has_anomalous  <- FALSE   # any anomalous crossings present?

  for (stage_name in names(x$details)) {

    # Parse stage index from name string (e.g. "stage.2" -> 2L)
    stage_num <- as.integer(sub("stage\\.", "", stage_name))
    stg_info  <- x$details[[stage_name]]

    if (is.null(stg_info$pairs)) next    # single-module stage: no crossings

    for (pair in stg_info$pairs) {

      sel       <- pair$selected_cut         # scalar: the chosen cut score
      proper    <- pair$proper_crossings     # numeric vector: all proper roots
      anomalous <- pair$anomalous_crossings  # numeric vector: all anomalous roots

      # -- Selected cut score ------------------------------------------------
      if (length(sel) == 1L && !is.na(sel)) {
        vline_rows[[length(vline_rows) + 1L]] <- data.frame(
          stage         = stage_num,
          xintercept    = sel,
          crossing_type = "selected",
          stringsAsFactors = FALSE
        )
        # Interpolate the TIF value at the crossing using the left module curve
        left_sub <- tif_plot[tif_plot$module == pair$left_module &
                              tif_plot$stage  == stage_num, ]
        if (nrow(left_sub) >= 2L) {
          y_cross <- stats::approx(x    = left_sub$theta,
                                   y    = left_sub$tif,
                                   xout = sel)$y
          if (!is.na(y_cross)) {
            point_rows[[length(point_rows) + 1L]] <- data.frame(
              stage  = stage_num,
              theta  = sel,
              tif    = y_cross,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      # -- Unselected proper crossings (when multiple exist) ------------------
      other_proper <- proper[!is.na(proper) & proper != sel]
      if (length(other_proper) > 0L) {
        has_unselected <- TRUE
        for (op in other_proper) {
          vline_rows[[length(vline_rows) + 1L]] <- data.frame(
            stage         = stage_num,
            xintercept    = op,
            crossing_type = "unselected",
            stringsAsFactors = FALSE
          )
        }
      }

      # -- Anomalous crossings (excluded from cut scores) --------------------
      if (show_anomalous &&
          length(anomalous) > 0L && !all(is.na(anomalous))) {
        has_anomalous <- TRUE
        for (ac in anomalous[!is.na(anomalous)]) {
          vline_rows[[length(vline_rows) + 1L]] <- data.frame(
            stage         = stage_num,
            xintercept    = ac,
            crossing_type = "anomalous",
            stringsAsFactors = FALSE
          )
        }
      }

    }  # end pair loop
  }  # end stage loop

  # Combine annotation rows into data frames (NULL if empty)
  if (length(vline_rows) > 0L) {
    vline_df <- do.call(rbind, vline_rows)
    vline_df$stage_label <- factor(paste0("Stage ", vline_df$stage),
                                   levels = stage_levels)
  } else {
    vline_df <- NULL
  }

  if (length(point_rows) > 0L) {
    point_df <- do.call(rbind, point_rows)
    point_df$stage_label <- factor(paste0("Stage ", point_df$stage),
                                   levels = stage_levels)
  } else {
    point_df <- NULL
  }

  # -- 3. Build plot caption describing vertical line types -------------------
  caption_parts <- "Vertical lines: solid black = selected cut score"
  if (has_unselected)
    caption_parts <- paste0(caption_parts,
                            " - grey dashed = unselected proper crossing")
  if (show_anomalous && has_anomalous)
    caption_parts <- paste0(caption_parts,
                            " - red dashed = anomalous crossing")

  # -- 4. Choose facet layout -------------------------------------------------
  # "vertical"   -> ncol = 1, one row per stage, stage 1 at the top
  # "horizontal" -> nrow = 1, one column per stage, stage 1 at the left
  facet_spec <- if (layout == "vertical") {
    ggplot2::facet_wrap(ggplot2::vars(.data$stage_label),
                        ncol   = 1L,
                        scales = "free_y")
  } else {
    ggplot2::facet_wrap(ggplot2::vars(.data$stage_label),
                        nrow   = 1L,
                        scales = "free_y")
  }

  # -- 5. Assemble the base ggplot --------------------------------------------
  p <- ggplot2::ggplot(
    data    = tif_plot,
    mapping = ggplot2::aes(x     = .data$theta,
                           y     = .data$tif,
                           color = .data$module_label,
                           group = .data$module_label)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    facet_spec +
    ggplot2::labs(
      x       = expression(theta),
      y       = "Test Information (TIF)",
      color   = "Module",
      title   = "TIF Curves and TIF-Crossing Cut Scores by Stage",
      caption = caption_parts
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "right",
      strip.background = ggplot2::element_rect(fill  = "grey92",
                                               color = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold", size = 11),
      plot.caption     = ggplot2::element_text(hjust  = 0,
                                               size   = 8,
                                               color  = "grey40")
    )

  # -- 6. Overlay vertical lines for each crossing type ----------------------
  if (!is.null(vline_df)) {

    # Selected cut scores -- solid black, linewidth 1.0
    sel_df <- vline_df[vline_df$crossing_type == "selected", ]
    if (nrow(sel_df) > 0L) {
      p <- p + ggplot2::geom_vline(
        data        = sel_df,
        mapping     = ggplot2::aes(xintercept = .data$xintercept),
        color       = "black",
        linetype    = "solid",
        linewidth   = 1.0,
        alpha       = 0.8,
        inherit.aes = FALSE
      )
    }

    # Unselected proper crossings -- dashed grey50
    unsel_df <- vline_df[vline_df$crossing_type == "unselected", ]
    if (nrow(unsel_df) > 0L) {
      p <- p + ggplot2::geom_vline(
        data        = unsel_df,
        mapping     = ggplot2::aes(xintercept = .data$xintercept),
        color       = "grey50",
        linetype    = "dashed",
        linewidth   = 0.7,
        alpha       = 0.8,
        inherit.aes = FALSE
      )
    }

    # Anomalous crossings -- dashed firebrick
    anom_df <- vline_df[vline_df$crossing_type == "anomalous", ]
    if (nrow(anom_df) > 0L) {
      p <- p + ggplot2::geom_vline(
        data        = anom_df,
        mapping     = ggplot2::aes(xintercept = .data$xintercept),
        color       = "firebrick",
        linetype    = "dashed",
        linewidth   = 0.7,
        alpha       = 0.8,
        inherit.aes = FALSE
      )
    }
  }

  # -- 7. Add open-circle dot at each selected crossing point -----------------
  # The dot sits at (theta_cut, TIF_left(theta_cut)) -- the exact intersection
  # of the two adjacent module TIF curves at the cut score theta.
  if (!is.null(point_df)) {
    p <- p + ggplot2::geom_point(
      data        = point_df,
      mapping     = ggplot2::aes(x = .data$theta, y = .data$tif),
      color       = "black",
      fill        = "white",
      shape       = 21L,       # open circle with border
      size        = 2.5,
      stroke      = 1.2,
      inherit.aes = FALSE
    )
  }

  # -- 8. Annotate each selected cut score with its theta value --------------
  # Label is placed above the crossing dot; if no dot exists, it is suppressed.
  if (label_cuts && !is.null(point_df)) {
    p <- p + ggplot2::geom_text(
      data        = point_df,
      mapping     = ggplot2::aes(x     = .data$theta,
                                 y     = .data$tif,
                                 label = sprintf("theta=%.3f", .data$theta)),
      vjust       = -0.9,     # slightly above the dot
      hjust       = 0.5,
      size        = 3.0,
      color       = "black",
      fontface    = "bold",
      inherit.aes = FALSE
    )
  }

  p    # return ggplot object (auto-printed when called interactively)
}
