#' Multistage-Adaptive Test (MST) Simulation
#'
#' @description
#' Simulates an MST administration for a panel of examinees. At each stage,
#' examinees are routed to a module based on their current ability estimate
#' (or a set of cut scores), and a final ability estimate is computed from the
#' accumulated responses across all administered stages. The function supports
#' all IRT models available in \pkg{irtQ} (1PLM, 2PLM, 3PLM, GRM, GPCM) and
#' all scoring methods available in \code{\link{est_score}}.
#'
#' The overall structure - separate \code{route_score} and \code{final_score}
#' list arguments, and the \code{route_map}/\code{module} matrix conventions -
#' was inspired by the \code{randomMST()} function in the \pkg{mstR} package
#' (Magis et al., 2017).
#'
#' @param x A data frame of item bank metadata following the irtQ format
#'   (columns: \code{id}, \code{cats}, \code{model}, \code{par.1},
#'   \code{par.2}, ...). Use \code{\link{shape_df}} to create this frame or
#'   \code{\link{bring.flexmirt}} / \code{\link{bring.bilog}} /
#'   \code{\link{bring.parscale}} / \code{\link{bring.mirt}} to import from
#'   IRT software output. The function calls \code{confirm_df(x)} internally
#'   to validate and normalise the metadata.
#'
#' @param route_map A binary square matrix defining the MST transition
#'   structure. A 1 at row \emph{i}, column \emph{j} means that a test taker
#'   can be routed from module \emph{i} to module \emph{j}. Equivalent to the
#'   \code{transMatrix} argument in \code{randomMST()} from \pkg{mstR}
#'   (Magis et al., 2017). See \code{\link{reval_mst}} for details.
#'
#' @param module A binary matrix mapping items (rows) to modules (columns).
#'   \code{module[j, m] = 1} means item \emph{j} belongs to module \emph{m}.
#'   Equivalent to the \code{modules} argument in \code{randomMST()} from
#'   \pkg{mstR} (Magis et al., 2017). Must have \code{nrow(module) == nrow(x)}
#'   and \code{ncol(module) == total number of modules in the panel}.
#'
#' @param theta A numeric vector of true ability values, one per examinee
#'   (length \emph{N}). Used to simulate item responses when \code{response}
#'   is \code{NULL}, and stored in the result for bias/RMSE computation.
#'   Either \code{theta} or \code{response} (or both) must be provided.
#'   If \code{theta} is omitted and only \code{response} is given,
#'   \code{true.theta} in the result will be \code{NULL} and bias/RMSE
#'   cannot be computed.
#'
#' @param response An optional \emph{N} x \emph{J} matrix of observed item
#'   responses, where \emph{N} is the number of examinees and \emph{J} is the
#'   total number of items in the item bank (\code{nrow(x)}). Responses for
#'   items not administered at a given stage may be \code{NA}. If \code{NULL}
#'   and \code{theta} is provided, responses are simulated internally via
#'   \code{\link{simdat}}. If both \code{theta} and \code{response} are
#'   provided, \code{response} is used as-is and \code{theta} is retained for
#'   evaluation (bias/RMSE computation) only. Providing a pre-generated
#'   \code{response} matrix (e.g., from \code{simdat(x, theta, D)}) is
#'   useful when you want to compare different routing methods or scoring
#'   options on the \emph{same} set of item responses.
#'
#' @param D A numeric scaling constant for the IRT model. Use \code{D = 1}
#'   (default) for the logistic metric or \code{D = 1.702} to approximate the
#'   normal ogive.
#'
#' @param ini_mod Either \code{NULL} (default) or a single positive integer
#'   index into the stage-1 module list
#'   (\code{panel_info(route_map)$config[[1]]}). When \code{NULL}, each
#'   examinee is independently and randomly assigned to one of the stage-1
#'   modules with equal probability - appropriate for MST panels that use
#'   spiralling or booklet-based stage-1 assignment, or for panels with only
#'   a single stage-1 module (where random assignment is equivalent to a
#'   fixed assignment). When an integer is supplied, all examinees begin with
#'   that stage-1 module; use \code{ini_mod = 1} to force all examinees to
#'   the first (and typically only) stage-1 module.
#'
#' @param route_method A character string specifying the adaptive routing
#'   method used to select the next module at each stage transition.
#'   \describe{
#'     \item{\code{"bmat"}}{B-matching: selects the module whose average item
#'       location (difficulty) is closest to the current ability estimate. For
#'       dichotomous items the \emph{b} parameter is used; for polytomous items
#'       (GRM, GPCM) the mean of the threshold / step parameters is used.}
#'     \item{\code{"mfi"}}{Maximum Fisher Information: selects the module that
#'       provides the highest test information function (TIF) value at the
#'       current ability estimate.}
#'     \item{\code{NULL}}{Cut-score routing: uses the cut scores supplied in
#'       \code{cut_score} to route examinees to the next module. Each
#'       examinee's current routing score is compared against the cut scores to
#'       determine which of the reachable next-stage modules to administer.}
#'   }
#'   Default is \code{"bmat"}.
#'
#' @param cut_score A list of numeric vectors specifying the routing cut scores
#'   used when \code{route_method = NULL}. The list should have
#'   \code{n.stage - 1} elements, one per stage transition. Element \emph{s}
#'   contains the cut scores for routing from stage \emph{s} to stage
#'   \emph{s}+1. For example, in a 1-3-3 MST, \code{cut_score = list(c(-0.5,
#'   0.5), c(-0.6, 0.6))} routes examinees whose stage-1 score is below
#'   \eqn{-0.5} to the easiest stage-2 module, between \eqn{-0.5} and
#'   \eqn{0.5} to the medium module, and above \eqn{0.5} to the hardest
#'   module. Ignored when \code{route_method} is \code{"bmat"} or
#'   \code{"mfi"}. Default is \code{NULL}.
#'
#' @param route_score A named list specifying the scoring method and options
#'   used to obtain intermediate ability estimates for routing decisions
#'   (applied after every stage except the last). Fields:
#'   \describe{
#'     \item{\code{method}}{Character. One of \code{"ML"} (maximum likelihood),
#'       \code{"WL"} (weighted likelihood; Warm, 1989), \code{"MLF"}
#'       (maximum likelihood with fences; Han, 2016), \code{"MAP"} (maximum
#'       a posteriori), \code{"EAP"} (expected a posteriori; Bock & Mislevy,
#'       1982), \code{"EAP.SUM"} (EAP summed scoring; Thissen et al., 1995),
#'       or \code{"INV.TCC"} (inverse test characteristic curve scoring; Lim
#'       et al., 2021). For \code{"EAP.SUM"} and \code{"INV.TCC"}, a
#'       sum-score-to-theta lookup table is pre-computed once per module before
#'       the simulation loop; routing theta is then obtained by a single
#'       named-vector lookup, making the approach efficient for large \eqn{N}.
#'       Default: \code{"ML"}.}
#'     \item{\code{range}}{Numeric vector of length 2: lower and upper bounds
#'       for the ability scale. Default: \code{c(-5, 5)}.}
#'     \item{\code{norm.prior}}{Numeric vector of length 2: mean and SD of the
#'       normal prior for MAP/EAP. Default: \code{c(0, 1)}.}
#'     \item{\code{nquad}}{Integer: number of quadrature points for EAP.
#'       Default: \code{41}.}
#'     \item{\code{tol}}{Numeric: convergence tolerance for iterative methods.
#'       Default: \code{1e-4}.}
#'     \item{\code{max.iter}}{Integer: maximum iterations. Default: \code{100}.}
#'     \item{\code{fence.a}}{Numeric: discrimination parameter of fence items
#'       for MLF. Default: \code{3.0}.}
#'     \item{\code{fence.b}}{Numeric vector of length 2 or \code{NULL}:
#'       difficulty parameters of the lower and upper fence items. If
#'       \code{NULL}, \code{range} is used. Default: \code{NULL}.}
#'     \item{\code{intpol}}{Logical: enable linear interpolation for
#'       \code{"INV.TCC"}. Default: \code{TRUE}.}
#'     \item{\code{range.tcc}}{Numeric vector of length 2: theta search range
#'       for \code{"INV.TCC"}. Default: \code{c(-7, 7)}.}
#'     \item{\code{max.it}}{Integer: maximum bisection iterations for
#'       \code{"INV.TCC"}. Default: \code{500L}.}
#'   }
#'   Unspecified fields take their defaults; the full default is
#'   \code{list(method = "ML", range = c(-5, 5), norm.prior = c(0, 1),}
#'   \code{nquad = 41L, tol = 1e-4, max.iter = 100L, fence.a = 3.0,}
#'   \code{fence.b = NULL, intpol = TRUE, range.tcc = c(-7, 7), max.it = 500L)}.
#'   Example: \code{route_score = list(method = "INV.TCC", range.tcc = c(-5, 5))}.
#'
#' @param final_score A named list specifying the scoring method and options
#'   for the final ability estimate (applied to all items accumulated across
#'   the entire administered pathway). Supports all fields in
#'   \code{route_score} plus:
#'   \describe{
#'     \item{\code{method}}{Character. Supports all \code{route_score} methods,
#'       plus \code{"EAP.SUM"} (EAP summed scoring; Thissen et al., 1995) and
#'       \code{"INV.TCC"} (inverse test characteristic curve scoring; Lim et
#'       al., 2021). Default: \code{"ML"}.}
#'     \item{\code{intpol}}{Logical: enable linear interpolation for
#'       \code{"INV.TCC"}. Default: \code{TRUE}.}
#'     \item{\code{range.tcc}}{Numeric vector of length 2: theta search range
#'       for \code{"INV.TCC"}. Default: \code{c(-7, 7)}.}
#'     \item{\code{max.it}}{Integer: maximum bisection iterations for
#'       \code{"INV.TCC"}. Default: \code{500}.}
#'   }
#'   The full default is equivalent to
#'   \code{list(method = "ML", range = c(-5, 5), norm.prior = c(0, 1),}
#'   \code{nquad = 41L, tol = 1e-4, max.iter = 100L, fence.a = 3.0, fence.b = NULL,}
#'   \code{intpol = TRUE, range.tcc = c(-7, 7), max.it = 500L)}.
#'
#' @param se Logical. Whether to compute and return the standard error of the
#'   final ability estimate. SE is computed for all \code{final_score} methods:
#'   for \code{"ML"}, \code{"WL"}, \code{"MLF"}, \code{"MAP"}, and \code{"EAP"}
#'   via the Fisher information / posterior variance; for \code{"EAP.SUM"} and
#'   \code{"INV.TCC"} via the posterior standard deviation stored in the
#'   pre-computed lookup table. Returns \code{NA} when \code{se = FALSE}, when
#'   all responses are missing, or when the sum score falls outside the
#'   estimable range (e.g. a perfect or zero score beyond \code{range.tcc}).
#'   Default is \code{TRUE}.
#'
#' @param missing A scalar specifying how missing (not-administered) responses
#'   are represented in \code{response}. Default is \code{NA}.
#'
#' @param verbose Logical. If \code{TRUE}, prints progress messages to the
#'   console during the simulation (approximately every 10\% of examinees).
#'   Default is \code{TRUE}.
#'
#' @param return_full_resp Logical. If \code{TRUE}, the \code{full.resp}
#'   element of the returned list is an \emph{N} x \emph{J} integer matrix,
#'   where \emph{N} is the number of examinees and \emph{J} = \code{nrow(x)}
#'   (all items in the item bank). Cell \code{[i, j]} contains examinee
#'   \emph{i}'s response to item \emph{j} if that item was administered, or
#'   \code{NA} if it was not. Row names are \code{"examinee.1"},
#'   \code{"examinee.2"}, etc.; column names are the item IDs from
#'   \code{x$id}. If \code{FALSE} (default), \code{full.resp} is \code{NULL}.
#'   Because this matrix can be large for large panels and many examinees, it
#'   is not populated by default.
#'
#' @details
#' \strong{Panel structure}: The function first calls \code{\link{panel_info}}
#' to identify stages, module membership, and valid pathways from
#' \code{route_map}. Module \emph{m}'s items are those where
#' \code{module[, m] == 1}.
#'
#' \strong{Response simulation}: If \code{response} is \code{NULL}, item
#' responses are simulated from \code{theta} using \code{\link{simdat}}.
#'
#' \strong{Stage-1 module assignment}: When \code{ini_mod = NULL} (default),
#' each examinee is independently assigned to a stage-1 module by uniform
#' random sampling from \code{panel_info(route_map)$config[[1]]}. For panels
#' with a single stage-1 module (the most common design), this is equivalent
#' to a fixed assignment. For panels with multiple stage-1 modules (e.g., a
#' 2-3-3 or 3-2-3 design), random assignment approximates the spiralling or
#' booklet-based allocation used in operational testing. When \code{ini_mod}
#' is an integer, all examinees start at the specified stage-1 module.
#'
#' \strong{Routing (stages 1 to n.stage - 1)}: After each non-final stage,
#' an intermediate ability estimate is obtained from the current stage's
#' responses using the \code{route_score} method. This estimate is used to
#' select the next-stage module:
#' \itemize{
#'   \item \code{"bmat"}: the reachable module (via \code{route_map}) whose
#'     mean item location is closest to the routing estimate is selected.
#'     For dichotomous items, the location is the \emph{b} parameter; for
#'     GRM/GPCM items, it is the mean of the threshold / step parameters.
#'   \item \code{"mfi"}: the reachable module with the highest test
#'     information at the routing estimate is selected.
#'   \item \code{NULL}: the routing estimate is compared against the cut
#'     scores in \code{cut_score[[s]]} (for the transition from stage
#'     \emph{s} to stage \emph{s}+1) to assign a rank, and the rank-th
#'     reachable module (ordered by module index) is administered. If the
#'     rank exceeds the number of reachable modules, the last module is used.
#' }
#'
#' \strong{Final scoring}: Responses from all administered stages are
#' concatenated, and the final ability estimate is computed using the
#' \code{final_score} method.
#'
#' @return An object of class \code{"run_mst"}, which is a named list
#'   containing:
#' \describe{
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{est.theta}}{A numeric vector of length \emph{N} containing
#'     the final ability estimate for each examinee.}
#'   \item{\code{se.theta}}{A numeric vector of length \emph{N} containing
#'     the standard error of the final ability estimate. All seven
#'     \code{final_score} methods return SE; for \code{"EAP.SUM"} and
#'     \code{"INV.TCC"}, SE is the posterior standard deviation from the
#'     pre-computed lookup table. \code{NA} when \code{se = FALSE}, when all
#'     responses are missing, or when the sum score falls outside the
#'     estimable range.}
#'   \item{\code{theta.route}}{An \emph{N} x \code{n.stage} numeric matrix
#'     of ability estimates. Columns 1 to \code{n.stage - 1} are the
#'     intermediate routing estimates; column \code{n.stage} is the final
#'     estimate (equal to \code{est.theta}).}
#'   \item{\code{path}}{An \emph{N} x \code{n.stage} integer matrix of
#'     administered module indices.}
#'   \item{\code{true.theta}}{The \code{theta} argument (true abilities), or
#'     \code{NULL} if not provided.}
#'   \item{\code{panel}}{The output of \code{\link{panel_info}(route_map)}.}
#'   \item{\code{D}}{The scaling constant used.}
#'   \item{\code{route.method}}{The \code{route_method} argument.}
#'   \item{\code{route.score}}{The fully-merged \code{route_score} list
#'     (with all defaults applied).}
#'   \item{\code{final.score}}{The fully-merged \code{final_score} list
#'     (with all defaults applied).}
#'   \item{\code{N}}{Number of examinees.}
#'   \item{\code{full.resp}}{Always present in the returned list. When
#'     \code{return_full_resp = TRUE}, an \emph{N} x \emph{J} integer matrix
#'     of all item responses, where \emph{N} is the number of examinees and
#'     \emph{J} = \code{nrow(x)} (total items in the item bank). Row names are
#'     \code{"examinee.1"}, \code{"examinee.2"}, etc.; column names are the
#'     item IDs from \code{x$id}. \code{NA} for items not administered to a
#'     given examinee. When \code{return_full_resp = FALSE} (default), this
#'     element is \code{NULL}.}
#' }
#'
#' @references
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
#' in a microcomputer environment. \emph{Applied Psychological Measurement,
#' 6}(4), 431-444. \doi{10.1177/014662168200600405}
#'
#' Han, K. T. (2016). Maximum likelihood score estimation method with fences
#' for short-length tests and computerized adaptive tests.
#' \emph{Applied Psychological Measurement, 40}(4), 289-301.
#' \doi{10.1177/0146621616631317}
#'
#' Lim, H., Davey, T., & Wells, C. S. (2021). A recursion-based analytical
#' approach to evaluate the performance of MST.
#' \emph{Journal of Educational Measurement, 58}(2), 154-178.
#' \doi{10.1111/jedm.12276}
#'
#' Magis, D., Yan, D., & von Davier, A. A. (2017). \emph{Computerized adaptive
#' and multistage testing with R: Using packages catR and mstR}. Springer.
#' \doi{10.1007/978-3-319-69218-0}
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995).
#' Item response theory for scores on tests including polytomous items with
#' ordered responses. \emph{Applied Psychological Measurement, 19}(1), 39-49.
#' \doi{10.1177/014662169501900105}
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. \emph{Psychometrika, 54}(3), 427-450.
#' \doi{10.1007/BF02294627}
#'
#' @seealso \code{\link{panel_info}}, \code{\link{reval_mst}},
#'   \code{\link{est_score}}, \code{\link{simdat}}, \code{\link{simMST}}
#'
#' @examples
#' \donttest{
#' ## ---------------------------------------------------------
#' ## Example 1: bmat routing with ML scoring (default)
#' ## Using the simMST 1-3-3 panel (7 modules, 8 items each)
#' ## ---------------------------------------------------------
#' # Load item bank metadata
#' x        <- simMST$item_bank
#' module   <- simMST$module
#' route_map <- simMST$route_map
#'
#' # True abilities for 500 examinees
#' set.seed(42)
#' theta_true <- rnorm(500, mean = 0, sd = 1)
#'
#' # Run MST simulation with bmat routing and ML scoring
#' result_bmat <- run_mst(
#'   x          = x,
#'   route_map  = route_map,
#'   module     = module,
#'   theta      = theta_true,
#'   D          = 1.702,
#'   route_method = "bmat",
#'   route_score  = list(method = "ML", range = c(-4, 4)),
#'   final_score  = list(method = "ML", range = c(-4, 4)),
#'   se = TRUE
#' )
#'
#' # Print simulation summary
#' print(result_bmat)
#'
#' # Final ability estimates
#' head(result_bmat$est.theta)
#'
#' # Module pathway taken by each examinee
#' head(result_bmat$path)
#'
#' ## ---------------------------------------------------------
#' ## Example 2: Cut-score routing with EAP routing and ML final scoring
#' ## EAP keeps the routing estimate finite even from a short Stage 1
#' ## module (e.g., a perfect or zero score), while ML avoids the shrinkage
#' ## that EAP's prior would otherwise introduce in the reported final score.
#' ## ---------------------------------------------------------
#' cut_score <- simMST$cut_score   # list(c(-0.45, 0.46), c(-0.42, 0.45))
#'
#' result_cut <- run_mst(
#'   x            = x,
#'   route_map    = route_map,
#'   module       = module,
#'   theta        = theta_true,
#'   D            = 1.702,
#'   route_method = NULL,
#'   cut_score    = cut_score,
#'   route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
#'   final_score  = list(method = "ML", range = c(-4, 4)),
#'   se = TRUE
#' )
#' print(result_cut)
#'
#' ## ---------------------------------------------------------
#' ## Example 3: Using a pre-generated response matrix
#' ## (useful when comparing multiple routing methods on the
#' ##  same simulated responses, or when using real observed data)
#' ## ---------------------------------------------------------
#' # Generate responses once from theta using simdat()
#' set.seed(123)
#' theta_true2 <- rnorm(500, mean = 0, sd = 1)
#' resp_matrix <- simdat(x = x, theta = theta_true2, D = 1.702)
#' # resp_matrix is an N x J matrix (N examinees x J total items in bank)
#'
#' # Now run run_mst() twice with the SAME responses but different routing methods
#' # Method A: bmat routing
#' result_A <- run_mst(
#'   x            = x,
#'   route_map    = route_map,
#'   module       = module,
#'   theta        = theta_true2,   # kept for bias/RMSE evaluation
#'   response     = resp_matrix,   # pre-generated N x J response matrix
#'   D            = 1.702,
#'   route_method = "bmat",
#'   route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
#'   final_score  = list(method = "ML", range = c(-4, 4)),
#'   se = FALSE
#' )
#'
#' # Method B: cut-score routing using the same responses
#' result_B <- run_mst(
#'   x            = x,
#'   route_map    = route_map,
#'   module       = module,
#'   theta        = theta_true2,
#'   response     = resp_matrix,   # identical response matrix -> fair comparison
#'   D            = 1.702,
#'   route_method = NULL,
#'   cut_score    = simMST$cut_score,
#'   route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
#'   final_score  = list(method = "ML", range = c(-4, 4)),
#'   se = FALSE
#' )
#'
#' # Compare RMSE of the two routing strategies
#' rmse_A <- sqrt(mean((result_A$est.theta - theta_true2)^2))
#' rmse_B <- sqrt(mean((result_B$est.theta - theta_true2)^2))
#' cat(sprintf("RMSE (bmat): %.4f\n", rmse_A))
#' cat(sprintf("RMSE (cut-score): %.4f\n", rmse_B))
#'
#' ## ---------------------------------------------------------
#' ## Example 4: TIF-based cut scores via find_cut()
#' ## find_cut() computes principled cut scores from TIF crossings.
#' ## This prevents anomalous routing under MFI-like heuristics.
#' ## ---------------------------------------------------------
#' # Derive cut scores from TIF crossings between adjacent modules
#' cut_result <- find_cut(
#'   x         = x,
#'   module    = module,
#'   route_map = route_map
#' )
#'
#' # Inspect the derived cut scores
#' print(cut_result)
#'
#' # Use the derived cut scores for routing (route_method = NULL)
#' result_findcut <- run_mst(
#'   x            = x,
#'   route_map    = route_map,
#'   module       = module,
#'   theta        = theta_true,
#'   D            = 1.702,
#'   route_method = NULL,
#'   cut_score    = cut_result$cut_score,
#'   route_score  = list(method = "EAP", norm.prior = c(0, 1), nquad = 41),
#'   final_score  = list(method = "ML", range = c(-4, 4)),
#'   se = TRUE
#' )
#' print(result_findcut)
#'
#' ## Visualise the TIF-based cut scores used for routing
#' plot(cut_result)
#'
#' ## ---------------------------------------------------------
#' ## Example 5: INV.TCC routing and final scoring
#' ## (pre-computed lookup tables; efficient for large N)
#' ## ---------------------------------------------------------
#' # INV.TCC can be used as both the routing and final scoring method.
#' # Before the examinee loop, run_mst() pre-computes a sum_score -> theta
#' # lookup table for every module (routing) and for every unique complete
#' # pathway (final scoring), consistent with the reval_mst() approach.
#' # SE is now returned correctly for INV.TCC (posterior SD from the table).
#'
#' result_inv <- run_mst(
#'   x            = x,
#'   route_map    = route_map,
#'   module       = module,
#'   theta        = theta_true,
#'   D            = 1.702,
#'   route_method = "bmat",
#'   route_score  = list(method = "INV.TCC", range.tcc = c(-7, 7)),
#'   final_score  = list(method = "INV.TCC", range.tcc = c(-7, 7)),
#'   se           = TRUE
#' )
#' print(result_inv)
#'
#' # SE is now populated (posterior SD from the pre-computed table)
#' head(result_inv$se.theta)
#'
#' # Compare RMSE: INV.TCC routing vs. ML routing (result_bmat from Example 1)
#' rmse_inv <- sqrt(mean((result_inv$est.theta  - theta_true)^2, na.rm = TRUE))
#' rmse_ml  <- sqrt(mean((result_bmat$est.theta - theta_true)^2, na.rm = TRUE))
#' cat(sprintf("RMSE (INV.TCC routing + scoring): %.4f\n", rmse_inv))
#' cat(sprintf("RMSE (ML    routing + scoring):   %.4f\n", rmse_ml))
#' }
#'
#' @export
run_mst <- function(x,
                    route_map,
                    module,
                    theta            = NULL,
                    response         = NULL,
                    D                = 1,
                    ini_mod          = NULL,
                    route_method     = "bmat",
                    cut_score        = NULL,
                    route_score      = list(method     = "ML",
                                           range      = c(-5, 5),
                                           norm.prior = c(0, 1),
                                           nquad      = 41L,
                                           tol        = 1e-4,
                                           max.iter   = 100L,
                                           fence.a    = 3.0,
                                           fence.b    = NULL,
                                           intpol     = TRUE,
                                           range.tcc  = c(-7, 7),
                                           max.it     = 500L),
                    final_score      = list(method     = "ML",
                                           range      = c(-5, 5),
                                           norm.prior = c(0, 1),
                                           nquad      = 41L,
                                           tol        = 1e-4,
                                           max.iter   = 100L,
                                           fence.a    = 3.0,
                                           fence.b    = NULL,
                                           intpol     = TRUE,
                                           range.tcc  = c(-7, 7),
                                           max.it     = 500L),
                    se               = TRUE,
                    missing          = NA,
                    verbose          = TRUE,
                    return_full_resp = FALSE) {

  # Capture the function call for inclusion in the result object
  call <- match.call()

  # --- Input validation ---

  # (A) Either theta or response must be provided
  if (is.null(theta) && is.null(response)) {
    stop("At least one of 'theta' or 'response' must be provided.",
         call. = FALSE)
  }

  # (B) Validate and normalise item metadata via confirm_df()
  x <- confirm_df(x)

  # (C) Check route_method
  valid_route_methods <- c("bmat", "mfi")
  if (!is.null(route_method) && !(route_method %in% valid_route_methods)) {
    stop("'route_method' must be \"bmat\", \"mfi\", or NULL.",
         call. = FALSE)
  }

  # (D) If route_method is NULL, cut_score must be provided
  if (is.null(route_method) && is.null(cut_score)) {
    stop("When 'route_method = NULL', 'cut_score' must be provided.",
         call. = FALSE)
  }

  # (E) Validate route_score method
  # EAP.SUM and INV.TCC are supported via pre-computed module-level lookup tables
  valid_route_score_methods <- c("ML", "WL", "MLF", "MAP", "EAP", "EAP.SUM", "INV.TCC")
  route_score_method <- if (is.null(route_score$method)) "ML" else route_score$method
  if (!(route_score_method %in% valid_route_score_methods)) {
    stop("'route_score$method' must be one of: ",
         paste(valid_route_score_methods, collapse = ", "), ".",
         call. = FALSE)
  }

  # (F) Validate final_score method
  valid_final_score_methods <- c("ML", "WL", "MLF", "MAP", "EAP", "EAP.SUM", "INV.TCC")
  final_score_method <- if (is.null(final_score$method)) "ML" else final_score$method
  if (!(final_score_method %in% valid_final_score_methods)) {
    stop("'final_score$method' must be one of: ",
         paste(valid_final_score_methods, collapse = ", "), ".",
         call. = FALSE)
  }

  # (G) Check module matrix dimensions
  if (nrow(module) != nrow(x)) {
    stop("'module' must have the same number of rows as 'x' (the item bank).",
         call. = FALSE)
  }

  # --- Merge user-supplied lists with defaults using modifyList() ---

  # Default arguments for routing scoring (ML, WL, MLF, MAP, EAP, EAP.SUM, INV.TCC)
  default_route <- list(
    method     = "ML",
    range      = c(-5, 5),
    norm.prior = c(0, 1),
    nquad      = 41L,
    tol        = 1e-4,
    max.iter   = 100L,
    fence.a    = 3.0,
    fence.b    = NULL,
    intpol     = TRUE,      # INV.TCC: use linear interpolation
    range.tcc  = c(-7, 7), # INV.TCC: theta search range for bisection
    max.it     = 500L       # INV.TCC: maximum bisection iterations
  )

  # Default arguments for final scoring (all route methods + EAP.SUM + INV.TCC)
  default_final <- default_route

  # Merge user input with defaults (user values override defaults)
  route_args <- utils::modifyList(default_route, route_score)
  final_args <- utils::modifyList(default_final, final_score)

  # --- Panel structure ---

  # Extract panel info: stage config, valid pathways, module/stage counts
  panel_data <- panel_info(route_map)
  n.stg      <- panel_data$n.stage          # total number of stages
  n.mod      <- panel_data$n.module         # modules per stage
  tn.mod     <- sum(n.mod)                  # total modules across all stages
  pathway    <- panel_data$pathway          # valid pathway matrix

  # Validate ini_mod: NULL (random per-examinee) or a single integer index
  if (!is.null(ini_mod)) {
    if (!is.numeric(ini_mod) || length(ini_mod) != 1L ||
        is.na(ini_mod) || ini_mod < 1L || ini_mod > n.mod[1L]) {
      stop(sprintf(
        "'ini_mod' must be NULL or a single integer between 1 and %d (number of stage-1 modules).",
        n.mod[1L]), call. = FALSE)
    }
    ini_mod <- as.integer(ini_mod)   # coerce to integer for safe indexing
  }

  # Pre-resolve the fixed start module (NULL when random assignment is used)
  # stage1_mods: all module indices that belong to stage 1
  stage1_mods <- panel_data$config[[1L]]          # integer vector, length = n.mod[1]
  fixed_start <- if (!is.null(ini_mod)) stage1_mods[ini_mod] else NULL
  # When fixed_start is NULL, each examinee is assigned a stage-1 module
  # by sample() at the start of the per-examinee loop (see below).

  # Validate cut_score length matches the number of stage transitions
  if (is.null(route_method) && !is.null(cut_score)) {
    if (length(cut_score) != (n.stg - 1)) {
      stop(sprintf(
        "'cut_score' must be a list of length %d (one vector per stage transition).",
        n.stg - 1), call. = FALSE)
    }
  }

  # --- Determine number of examinees ---
  if (!is.null(theta)) {
    N <- length(theta)
  } else {
    N <- nrow(response)
  }

  # --- Build item metadata list for each module ---
  # item_mod[[m]] = item metadata data frame for module m
  item_mod <- purrr::map(seq_len(tn.mod), function(m) {
    # Select items belonging to module m (module matrix column m == 1)
    x[module[, m] == 1, ] %>%
      tibble::remove_rownames()
  })
  names(item_mod) <- paste0("m.", seq_len(tn.mod))

  # Verify every module has at least one item assigned
  mod_item_counts <- colSums(module)
  if (any(mod_item_counts == 0)) {
    empty_mods <- which(mod_item_counts == 0)
    stop(sprintf("Module(s) %s have no items assigned (column sum = 0 in 'module').",
                 paste(empty_mods, collapse = ", ")), call. = FALSE)
  }

  # --- Pre-compute breakdown objects, idx, max.cats for each module ---
  # This pattern mirrors est_score_1core(): move breakdown/idxfinder outside
  # the per-examinee loop.

  elm_mod      <- vector("list", tn.mod)   # breakdown() output per module
  max_cats_mod <- integer(tn.mod)          # max score category per module
  idx_mod      <- vector("list", tn.mod)   # idxfinder() output per module

  for (m in seq_len(tn.mod)) {
    elm_mod[[m]]      <- breakdown(item_mod[[m]])   # convert metadata to list
    max_cats_mod[[m]] <- max(item_mod[[m]]$cats)    # max score category
    idx_mod[[m]]      <- idxfinder(elm_mod[[m]])    # DRM/PRM item indices
  }

  # --- Pre-compute MLF fence item metadata for routing (if MLF is used) ---
  # When route_args$method == "MLF", two synthetic fence items are appended
  # to the module's item metadata before scoring. Pre-compute the fence
  # metadata once, because fence.a and fence.b are fixed across all modules
  # and all examinees.
  elm_mod_mlf_route <- NULL     # will be a list of tn.mod elm_item objects
  max_cats_mlf_route <- NULL
  idx_mlf_route <- NULL

  if (route_args$method == "MLF") {
    # Determine fence.b: if NULL, use route_args$range
    fence_b_route <- if (is.null(route_args$fence.b)) route_args$range
                     else route_args$fence.b

    # Build the fence item metadata data frame (two 3PLM items)
    x_fence_route <- shape_df(
      par.drm  = list(a = rep(route_args$fence.a, 2),
                      b = fence_b_route,
                      g = rep(0, 2)),
      item.id  = c("fence.lower", "fence.upper"),
      cats     = 2,
      model    = "3PLM"
    )

    # Build augmented elm_item for each module (module items + 2 fence items)
    elm_mod_mlf_route  <- vector("list", tn.mod)
    max_cats_mlf_route <- integer(tn.mod)
    idx_mlf_route      <- vector("list", tn.mod)

    for (m in seq_len(tn.mod)) {
      # Coerce id to character so module items (id may be integer/numeric)
      # bind cleanly with the character fence item IDs
      x_aug_m             <- dplyr::bind_rows(
        dplyr::mutate(item_mod[[m]], id = as.character(.data$id)),
        x_fence_route
      )
      elm_mod_mlf_route[[m]]  <- breakdown(x_aug_m)
      max_cats_mlf_route[[m]] <- max(x_aug_m$cats)
      idx_mlf_route[[m]]      <- idxfinder(elm_mod_mlf_route[[m]])
    }
  }

  # --- Pre-compute EAP quadrature weights for routing (if EAP is used) ---
  popdist_route <- NULL
  if (route_args$method == "EAP") {
    popdist_route <- gen.weight(
      n     = route_args$nquad,
      dist  = "norm",
      mu    = route_args$norm.prior[1],
      sigma = route_args$norm.prior[2]
    )
  }

  # --- Pre-compute EAP quadrature weights for final scoring (if EAP final) ---
  popdist_final <- NULL
  if (final_args$method == "EAP") {
    popdist_final <- gen.weight(
      n     = final_args$nquad,
      dist  = "norm",
      mu    = final_args$norm.prior[1],
      sigma = final_args$norm.prior[2]
    )
  }

  # --- Pre-compute sum_score -> theta/SE lookup tables for EAP.SUM/INV.TCC routing ---
  # One table per module (tn.mod tables total); cost is O(tn.mod), independent of N.
  # During the examinee loop, routing theta = named-vector lookup, not a function call.
  #
  # route_tables[[m]]$theta  named numeric: name = "0","1",...; value = theta estimate
  # route_tables[[m]]$se     named numeric: name = "0","1",...; value = SE estimate
  #
  # Tables are built by calling inv_tcc() or eap_sum() with a 1-row dummy response
  # matrix (all zeros). Only $score.table is used - it covers all possible sum scores
  # and is computed independently of the actual response data passed in 'data'.
  route_tables <- NULL   # remains NULL when method is ML/WL/MLF/MAP/EAP

  if (route_args$method %in% c("EAP.SUM", "INV.TCC")) {
    route_tables <- vector("list", tn.mod)
    names(route_tables) <- names(item_mod)   # "m.1", "m.2", ...

    for (m in seq_len(tn.mod)) {
      n_items_m <- nrow(item_mod[[m]])
      # Dummy 1-row all-zero matrix: triggers full score.table computation
      dummy_m <- as.data.frame(matrix(0L, nrow = 1L, ncol = n_items_m))

      if (route_args$method == "INV.TCC") {
        tbl_m <- inv_tcc(
          x         = item_mod[[m]],
          data      = dummy_m,
          D         = D,
          intpol    = route_args$intpol,
          range.tcc = route_args$range.tcc,
          tol       = route_args$tol,
          max.it    = route_args$max.it
        )$score.table

      } else {   # EAP.SUM
        tbl_m <- eap_sum(
          x          = item_mod[[m]],
          data       = dummy_m,
          norm.prior = route_args$norm.prior,
          nquad      = route_args$nquad,
          weights    = NULL,
          D          = D
        )$score.table
      }

      # Build named lookup vectors (name = sum score as string)
      route_tables[[m]] <- list(
        theta = setNames(tbl_m$est.theta, as.character(tbl_m$sum.score)),
        se    = setNames(tbl_m$se.theta,  as.character(tbl_m$sum.score))
      )
    }
  }

  # --- Pre-compute sum_score -> theta/SE lookup tables for EAP.SUM/INV.TCC final scoring ---
  # One table per unique complete pathway; cost is O(n_unique_pathways), independent of N.
  # Fixes the se_theta=NA bug: both theta AND SE are stored and returned from the table.
  #
  # Key format: paste(module_indices_stage1_to_stageN, collapse = "_")  e.g. "1_3_6"
  # final_tables[["1_3_6"]]$theta  named numeric: sum_score -> theta
  # final_tables[["1_3_6"]]$se     named numeric: sum_score -> SE
  final_tables <- NULL   # remains NULL when method is ML/WL/MLF/MAP/EAP

  if (final_args$method %in% c("EAP.SUM", "INV.TCC")) {
    # Enumerate all unique complete pathways (rows of panel_data$pathway)
    uni_path  <- unique(pathway)   # matrix: each row is one complete path
    path_keys <- apply(uni_path, 1L, paste, collapse = "_")
    final_tables <- vector("list", length(path_keys))
    names(final_tables) <- path_keys

    for (p in seq_len(nrow(uni_path))) {
      # Accumulate item metadata for all modules along this complete pathway
      path_mods <- uni_path[p, ]   # module indices (length = n.stg)
      x_path_p  <- dplyr::bind_rows(
        purrr::map(path_mods, ~ item_mod[[.x]])
      ) %>% tibble::remove_rownames()

      n_items_p <- nrow(x_path_p)
      dummy_p   <- as.data.frame(matrix(0L, nrow = 1L, ncol = n_items_p))

      if (final_args$method == "INV.TCC") {
        tbl_p <- inv_tcc(
          x         = x_path_p,
          data      = dummy_p,
          D         = D,
          intpol    = final_args$intpol,
          range.tcc = final_args$range.tcc,
          tol       = final_args$tol,
          max.it    = final_args$max.it
        )$score.table

      } else {   # EAP.SUM
        tbl_p <- eap_sum(
          x          = x_path_p,
          data       = dummy_p,
          norm.prior = final_args$norm.prior,
          nquad      = final_args$nquad,
          weights    = NULL,
          D          = D
        )$score.table
      }

      final_tables[[path_keys[p]]] <- list(
        theta = setNames(tbl_p$est.theta, as.character(tbl_p$sum.score)),
        se    = setNames(tbl_p$se.theta,  as.character(tbl_p$sum.score))
      )
    }
  }

  # --- Generate or prepare the response matrix ---
  if (is.null(response)) {
    # Simulate full item bank responses from true theta using simdat()
    # simdat() returns an N x J response matrix directly (no list wrapper)
    resp_mat <- simdat(x = x, theta = theta, D = D)
  } else {
    resp_mat <- data.matrix(response)
    # Replace any custom missing value with NA
    if (!is.na(missing)) {
      resp_mat[resp_mat == missing] <- NA
    }
  }

  # Confirm dimensions
  if (nrow(resp_mat) != N || ncol(resp_mat) != nrow(x)) {
    stop(sprintf(
      "'response' must be an N x J matrix where N = %d (examinees) and J = %d (items).",
      N, nrow(x)
    ), call. = FALSE)
  }

  # --- Initialise result containers ---

  # Final ability estimates (length N)
  est_theta  <- numeric(N)
  # SE of final estimates (length N; NA if SE not computed)
  se_theta   <- rep(NA_real_, N)
  # N x n.stg routing estimate matrix (last column = final estimate)
  theta_route_mat <- matrix(NA_real_, nrow = N, ncol = n.stg)
  colnames(theta_route_mat) <- paste0("stage.", 1:n.stg)
  # N x n.stg module path matrix
  path_mat   <- matrix(NA_integer_, nrow = N, ncol = n.stg)
  colnames(path_mat) <- paste0("stage.", 1:n.stg)

  # --- Initialise the full response matrix (always allocated) ---
  # Dimensions: N x J  (rows = examinees, columns = items in item bank)
  # full_resp is NULL when return_full_resp = FALSE; matrix otherwise.
  J         <- nrow(x)             # total items in the item bank
  full_resp <- if (return_full_resp) {
    m <- matrix(NA_integer_, nrow = N, ncol = J)   # N examinees x J items
    rownames(m) <- paste0("examinee.", seq_len(N)) # row names: examinee index
    colnames(m) <- x$id                            # column names: item IDs
    m
  } else {
    NULL                                           # not requested; keep as NULL
  }

  # --- Verbose progress reporting ---
  if (verbose) {
    message(sprintf(
      "[run_mst] Panel: %d stages, %d modules (%s design)",
      n.stg, tn.mod,
      paste(n.mod, collapse = "-")
    ))
    if (is.null(response)) {
      message(sprintf("[run_mst] Simulating responses for %d examinees...", N))
    }
    message(sprintf("[run_mst] Scoring %d examinees...", N))
    # Report interval: approximately every 10% of examinees (min 1)
    report_every <- max(1L, N %/% 10L)
  }

  # --- Main loop over examinees ---
  for (i in seq_len(N)) {

    # Print progress at start, every report_every examinees, and at end
    if (verbose && (i == 1L || i %% report_every == 0L || i == N)) {
      message(sprintf("  Processing examinee %d / %d ...", i, N))
    }

    # Determine starting module for this examinee:
    #   fixed_start non-NULL -> same module for every examinee (ini_mod was integer)
    #   fixed_start NULL     -> uniformly sample from all stage-1 modules
    current_mod <- if (!is.null(fixed_start)) fixed_start
                   else sample(stage1_mods, size = 1L)
    # Accumulate item indices and response vector across all administered stages
    items_acc  <- integer(0)       # item row indices in x (item bank)
    resp_acc   <- integer(0)       # corresponding response values

    # --- Iterate over stages ---
    for (s in seq_len(n.stg)) {

      # ---- Select module for this stage ----
      if (s == 1) {
        # Stage 1: module already determined above (current_mod)
        mod_s <- current_mod

      } else {
        # Stages 2..n.stg: select next module based on route_method
        # Identify modules reachable from current_mod via route_map
        next_possible <- which(route_map[current_mod, ] == 1)
        # next_possible is ordered by module index (ascending)

        if (length(next_possible) == 1L) {
          # Only one reachable module - no choice to make
          mod_s <- next_possible[1L]

        } else if (!is.null(route_method) && route_method == "bmat") {
          # b-matching: select module with mean location closest to current theta
          # theta_route_mat[i, s-1] holds the routing estimate from stage s-1
          theta_prev <- theta_route_mat[i, s - 1L]
          mod_locs <- vapply(next_possible, function(m) {
            mean(mean_loc(item_mod[[m]]))   # mean location of module m's items
          }, numeric(1L))
          mod_s <- next_possible[which.min(abs(mod_locs - theta_prev))]

        } else if (!is.null(route_method) && route_method == "mfi") {
          # MFI: select module with highest TIF at current routing estimate
          theta_prev <- theta_route_mat[i, s - 1L]
          tif_vals <- vapply(next_possible, function(m) {
            # info() returns a list; $tif is the test information value
            info(x = item_mod[[m]], theta = theta_prev, D = D, tif = TRUE)$tif
          }, numeric(1L))
          mod_s <- next_possible[which.max(tif_vals)]

        } else {
          # Cut-score routing: use cut_score[[s-1]] for this transition
          # cut_score[[s-1]] has length = (number of categories - 1) for stage s
          theta_prev <- theta_route_mat[i, s - 1L]
          cut_s      <- cut_score[[s - 1L]]
          # give_path() assigns a rank (1, 2, ..., ncats) based on theta_prev
          rank_s <- give_path(score = theta_prev, cut_sc = cut_s)$path
          # Clamp rank to valid range [1, length(next_possible)]
          rank_s <- max(1L, min(rank_s, length(next_possible)))
          mod_s  <- next_possible[rank_s]
        }

        # Update current module tracker
        current_mod <- mod_s
      }

      # ---- Record module in path matrix ----
      path_mat[i, s] <- mod_s

      # ---- Accumulate administered items and responses ----
      items_s   <- which(module[, mod_s] == 1)   # item indices for module mod_s
      resp_s    <- resp_mat[i, items_s]           # this examinee's responses
      items_acc <- c(items_acc, items_s)          # append to accumulated list
      resp_acc  <- c(resp_acc, resp_s)            # append to accumulated responses

      # Record responses into the full N x J matrix (if requested)
      if (return_full_resp) {
        # Row i = examinee i; columns items_s = items of module mod_s
        full_resp[i, items_s] <- as.integer(resp_s)
      }

      # ---- Score responses for routing (all stages except the final) ----
      if (s < n.stg) {
        # Use route_score method on THIS STAGE'S items only
        resp_s_num <- as.integer(resp_s)       # integer 0-based responses

        # Identify non-missing items for this stage
        na_mask_s <- !is.na(resp_s_num)
        if (!any(na_mask_s)) {
          # All responses missing: carry forward last theta or use 0
          theta_s <- if (s == 1L) 0 else theta_route_mat[i, s - 1L]
        } else {
          # Subset to non-missing items for this stage
          resp_sub <- resp_s_num[na_mask_s]
          na_pos   <- which(na_mask_s)

          if (route_args$method %in% c("EAP.SUM", "INV.TCC")) {
            # ----- EAP.SUM / INV.TCC: pre-computed module-level table lookup -----
            # Sum all non-NA responses for this stage; look up theta from the
            # pre-computed table for module mod_s.
            sum_s   <- sum(resp_sub)   # sum of observed responses (NA already excluded)
            theta_s <- route_tables[[mod_s]]$theta[as.character(sum_s)]
            # Fallback when sum_s is outside the estimable range (e.g. all wrong / all right
            # near boundary with extreme range.tcc). Carry forward previous theta or use 0.
            if (is.na(theta_s)) {
              theta_s <- if (s == 1L) 0 else theta_route_mat[i, s - 1L]
            }

          } else {
            # ----- ML / WL / MLF / MAP / EAP: est_score_indiv() -----

            # Choose pre-computed elm_item based on whether MLF is used
            if (route_args$method == "MLF") {
              # Augmented elm_item (module items + 2 fence items)
              elm_s        <- elm_mod_mlf_route[[mod_s]]
              max_cats_s   <- max_cats_mlf_route[[mod_s]]
              idx_full_s   <- idx_mlf_route[[mod_s]]
              # Append fence responses: lower fence = 1 (correct), upper = 0
              resp_sub <- c(resp_sub, 1L, 0L)
              # Extend na_pos to include both fence items (always observed)
              na_pos_fence <- c(na_pos, (length(na_mask_s) + 1L), (length(na_mask_s) + 2L))
              na_pos <- na_pos_fence
            } else {
              elm_s      <- elm_mod[[mod_s]]
              max_cats_s <- max_cats_mod[[mod_s]]
              idx_full_s <- idx_mod[[mod_s]]
            }

            # Map full-item indices to the non-NA subset (mirrors est_score_1core)
            if (all(na_mask_s) && route_args$method != "MLF") {
              idx_drm_s <- idx_full_s$idx.drm
              idx_prm_s <- idx_full_s$idx.prm
            } else {
              idx_drm_s <- if (!is.null(idx_full_s$idx.drm)) {
                matched <- which(na_pos %in% idx_full_s$idx.drm)
                if (length(matched) == 0L) NULL else matched
              } else NULL
              idx_prm_s <- if (!is.null(idx_full_s$idx.prm)) {
                matched <- which(na_pos %in% idx_full_s$idx.prm)
                if (length(matched) == 0L) NULL else matched
              } else NULL
            }

            # Subset elm_item to observed items
            elm_sub_s       <- elm_s
            elm_sub_s$pars  <- elm_s$pars[seq_along(resp_sub), , drop = FALSE]
            elm_sub_s$model <- elm_s$model[seq_along(resp_sub)]
            elm_sub_s$cats  <- elm_s$cats[seq_along(resp_sub)]

            # Estimate routing theta using est_score_indiv()
            route_result <- est_score_indiv(
              resp_vec   = resp_sub,
              elm_item   = elm_sub_s,
              max.cats   = max_cats_s,
              idx.drm    = idx_drm_s,
              idx.prm    = idx_prm_s,
              D          = D,
              method     = route_args$method,
              range      = route_args$range,
              norm.prior = route_args$norm.prior,
              nquad      = route_args$nquad,
              weights    = NULL,
              tol        = route_args$tol,
              max.iter   = route_args$max.iter,
              se         = FALSE,       # SE not needed for routing
              stval.opt  = 1L,
              ji         = (route_args$method == "WL"),
              obs.sum    = NULL,
              popdist    = popdist_route   # pre-computed for EAP; NULL otherwise
            )
            theta_s <- route_result$est.theta
          }   # end EAP.SUM/INV.TCC vs ML/WL/...
        }

        # Store routing estimate for this stage
        theta_route_mat[i, s] <- theta_s

      } else {
        # ---- Final stage: compute final ability estimate ----
        # Use ALL accumulated items and responses

        # Handle missing values
        resp_acc_num <- as.integer(resp_acc)
        na_mask_acc  <- !is.na(resp_acc_num)

        if (!any(na_mask_acc)) {
          # All responses missing across entire path
          est_theta[i]   <- NA_real_
          se_theta[i]    <- NA_real_
          theta_route_mat[i, s] <- NA_real_
          next   # skip to next examinee
        }

        # Build item metadata for the full administered path
        x_path <- x[items_acc, ]  %>% tibble::remove_rownames()

        if (final_args$method %in% c("EAP.SUM", "INV.TCC")) {
          # ----- EAP.SUM or INV.TCC: pre-computed pathway lookup table -----
          # theta and SE are retrieved by a single named-vector lookup.
          # This also fixes the se_theta=NA bug: SE is now returned correctly.
          path_key  <- paste(path_mat[i, ], collapse = "_")   # e.g. "1_3_6"
          sum_total <- sum(resp_acc_num, na.rm = TRUE)         # total accumulated sum score

          # Safety check: pathway must exist in pre-computed tables
          if (is.null(final_tables[[path_key]])) {
            est_theta[i]          <- NA_real_
            se_theta[i]           <- NA_real_
            theta_route_mat[i, s] <- NA_real_
          } else {
            theta_i <- final_tables[[path_key]]$theta[as.character(sum_total)]
            se_i    <- final_tables[[path_key]]$se[as.character(sum_total)]

            # NA can occur at extreme sum scores (outside estimable range)
            est_theta[i]          <- if (!is.na(theta_i)) theta_i else NA_real_
            se_theta[i]           <- if (se && !is.na(se_i)) se_i else NA_real_
            theta_route_mat[i, s] <- est_theta[i]
          }

        } else {
          # ----- ML / WL / MLF / MAP / EAP via est_score_indiv() -----
          resp_sub_acc <- resp_acc_num[na_mask_acc]
          na_pos_acc   <- which(na_mask_acc)

          # Build elm_item for the accumulated path
          if (final_args$method == "MLF") {
            # Determine fence.b
            fence_b_final <- if (is.null(final_args$fence.b)) final_args$range
                             else final_args$fence.b
            # Build fence item metadata
            x_fence_final <- shape_df(
              par.drm  = list(a = rep(final_args$fence.a, 2),
                              b = fence_b_final,
                              g = rep(0, 2)),
              item.id  = c("fence.lower", "fence.upper"),
              cats     = 2L,
              model    = "3PLM"
            )
            # Coerce id to character so path items (id may be integer/numeric)
            # bind cleanly with the character fence item IDs
            x_path_aug <- dplyr::bind_rows(
              dplyr::mutate(x_path, id = as.character(.data$id)),
              x_fence_final
            )
            elm_path   <- breakdown(x_path_aug)
            max_cats_p <- max(x_path_aug$cats)
            idx_path   <- idxfinder(elm_path)
            # Append fence responses (lower = 1, upper = 0)
            resp_sub_acc <- c(resp_sub_acc, 1L, 0L)
            # Update na_pos for fences (always observed)
            n_items_acc <- length(na_mask_acc)
            na_pos_acc  <- c(na_pos_acc, n_items_acc + 1L, n_items_acc + 2L)
          } else {
            elm_path   <- breakdown(x_path)
            max_cats_p <- max(x_path$cats)
            idx_path   <- idxfinder(elm_path)
          }

          # Map indices to non-NA subset
          if (all(na_mask_acc) && final_args$method != "MLF") {
            idx_drm_p <- idx_path$idx.drm
            idx_prm_p <- idx_path$idx.prm
          } else {
            idx_drm_p <- if (!is.null(idx_path$idx.drm)) {
              matched <- which(na_pos_acc %in% idx_path$idx.drm)
              if (length(matched) == 0L) NULL else matched
            } else NULL
            idx_prm_p <- if (!is.null(idx_path$idx.prm)) {
              matched <- which(na_pos_acc %in% idx_path$idx.prm)
              if (length(matched) == 0L) NULL else matched
            } else NULL
          }

          # Subset elm_item to observed items
          elm_sub_p       <- elm_path
          elm_sub_p$pars  <- elm_path$pars[seq_along(resp_sub_acc), , drop = FALSE]
          elm_sub_p$model <- elm_path$model[seq_along(resp_sub_acc)]
          elm_sub_p$cats  <- elm_path$cats[seq_along(resp_sub_acc)]

          final_result <- est_score_indiv(
            resp_vec   = resp_sub_acc,
            elm_item   = elm_sub_p,
            max.cats   = max_cats_p,
            idx.drm    = idx_drm_p,
            idx.prm    = idx_prm_p,
            D          = D,
            method     = final_args$method,
            range      = final_args$range,
            norm.prior = final_args$norm.prior,
            nquad      = final_args$nquad,
            weights    = NULL,
            tol        = final_args$tol,
            max.iter   = final_args$max.iter,
            se         = se,
            stval.opt  = 1L,
            ji         = (final_args$method == "WL"),
            obs.sum    = NULL,
            popdist    = popdist_final
          )

          est_theta[i]          <- final_result$est.theta
          se_theta[i]           <- if (se) final_result$se.theta else NA_real_
          theta_route_mat[i, s] <- final_result$est.theta
        }
      }  # end of stage s
    }  # end of stage loop
  }  # end of examinee loop

  # --- Assemble result list ---
  rst <- list(
    call             = call,
    est.theta        = est_theta,           # final ability estimates (length N)
    se.theta         = se_theta,            # SE of final estimates (length N)
    theta.route      = theta_route_mat,     # N x n.stg routing + final estimates
    path             = path_mat,            # N x n.stg administered module indices
    true.theta       = theta,               # true abilities (NULL if not provided)
    panel            = panel_data,          # panel_info() output
    D                = D,                   # scaling constant
    route.method     = route_method,        # routing method used
    route.score      = route_args,          # merged route scoring arguments
    final.score      = final_args,          # merged final scoring arguments
    N                = N,                   # number of examinees
    full.resp        = full_resp            # N x J matrix or NULL
  )

  # Assign S3 class for print dispatch
  class(rst) <- "run_mst"

  if (verbose) message("[run_mst] Done.")

  rst
}

# ---------------------------------------------------------------------------
# mean_loc() - Internal helper for bmat routing
#
# Computes the mean location parameter for each item in an item metadata
# data frame. Used by the bmat (b-matching) routing method to select the
# module whose average difficulty is closest to the current theta estimate.
#
# For dichotomous items (1PLM, 2PLM, 3PLM, DRM):
#   Returns the b (difficulty) parameter = par.2.
#
# For polytomous items (GRM, GPCM):
#   Returns the mean of all threshold / step parameters
#   (par.2 through par.[cats]), i.e., mean(b_1, b_2, ..., b_{K-1}).
#   This average location captures the overall difficulty of the item.
#
# @param item_meta  A data frame conforming to the irtQ item metadata format
#                   (columns: id, cats, model, par.1, par.2, ...).
# @return           A numeric vector of length nrow(item_meta), one mean
#                   location value per item.
# ---------------------------------------------------------------------------
mean_loc <- function(item_meta) {
  
  # Identify all parameter columns (par.1, par.2, ...)
  par_cols <- grep("^par\\.", names(item_meta), value = TRUE)
  pars     <- item_meta[, par_cols, drop = FALSE]
  
  # Compute mean location for each item
  locs <- purrr::map_dbl(seq_len(nrow(item_meta)), function(i) {
    model_i <- item_meta$model[i]   # IRT model for this item
    cats_i  <- item_meta$cats[i]    # number of score categories
    
    if (model_i %in% c("1PLM", "2PLM", "3PLM", "DRM")) {
      # Dichotomous: location = b parameter (par.2, 2nd column of pars)
      as.numeric(pars[i, 2])
    } else {
      # GRM / GPCM: location = mean of threshold/step parameters
      # par.2 through par.[cats] = columns 2 to cats_i
      thresh_vals <- as.numeric(pars[i, 2:cats_i])
      mean(thresh_vals, na.rm = TRUE)
    }
  })
  
  locs   # return mean location vector
}


