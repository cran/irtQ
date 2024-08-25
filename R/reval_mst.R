#' Recursion-based MST evaluation method
#'
#' @description This function evaluates the measurement precision and bias in
#' Multistage-adaptive Test (MST) panels using a recursion-based evaluation
#' method introduced by Lim et al. (2020). This function computes conditional
#' biases and standard errors of measurement (CSEMs) across a range of IRT
#' ability levels, facilitating efficient and accurate MST panel
#' assessments without extensive simulations.
#'
#' @param x A data frame containing the metadata for the item bank, which includes
#' important information for each item such as the number of score categories and the
#' IRT model applied. This metadata is essential for evaluating the MST panel, with
#' items selected based on the specifications in the \code{module} argument. To construct
#' this item metadata efficiently, the \code{\link{shape_df}} function is recommended.
#' Further details on utilizing item bank metadata along with \code{module} for
#' MST panel evaluation are provided below.
#' @param D A scaling factor in IRT models to make the logistic function as
#' close as possible to the normal ogive function (if set to 1.7). Default is 1.
#' @param route_map A binary square matrix that defines the MST structure, illustrating
#' transitions between modules and stages. This concept and structure are inspired by
#' the \code{transMatrix} argument in the \code{randomMST} function
#' from the \pkg{mstR} package (Magis et al., 2017), which provides a framework for representing MST pathways.
#' For comprehensive understanding and examples of constructing \code{route_map},
#' refer to the \pkg{mstR} package (Magis et al., 2017) documentation. Also see below for details.
#' @param module A binary matrix that maps items from the item bank specified in \code{x}
#' to modules within the MST framework. This parameter's structure is analogous to the
#' \code{modules} argument in the \code{randomMST}  function of the \pkg{mstR}
#' package, enabling precise item-to-module assignments for MST configurations.
#' For detailed instructions on creating and utilizing the \code{module} matrix
#' effectively, consult the documentation of the \pkg{mstR} package (Magis et al., 2017).
#' Also see below for details.
#' @param cut_score A list defining cut scores for routing test takers through MST stages.
#' Each list element is a vector of cut scores for advancing participants to subsequent
#' stage modules. In a 1-3-3 MST configuration, for example, \code{cut_score} might be
#' defined as \code{cut_score = list(c(-0.5, 0.5), c(-0.6, 0.6))}, where c(-0.5, 0.5)
#' are thresholds for routing from the first to the second stage, and c(-0.6, 0.6)
#' for routing from the second to the third stage. This setup facilitates tailored
#' test progression based on performance. Further examples and explanations are
#' available below.
#' @param theta A vector of ability levels (theta) at which the MST panel's performance
#' is assessed. This allows for the evaluation of measurement precision and bias
#' across a continuum of ability levels. The default range is \code{theta = seq(-5, 5, 0.1)}.
#' @param intpol A logical value to enable linear interpolation in the inverse
#' test characteristic curve (TCC) scoring, facilitating ability estimate approximation
#' for observed sum scores not directly obtainable from the TCC, such as those below
#' the sum of item guessing parameters. Default is TRUE, applying interpolation
#' to bridge gaps in the TCC. Refer to \code{\link{est_score}} for more details and
#' consult Lim et al. (2020) for insights into the interpolation technique
#' within inverse TCC scoring.
#' @param range.tcc A vector to define the range of ability estimates for
#' inverse TCC scoring, expressed as the two numeric values for lower and upper
#' bounds. Default is to c(-7, 7).
#' @param tol A numeric value of the convergent tolerance for the  inverse TCC
#' scoring. For the inverse TCC scoring, the bisection method is used for
#' optimization. Default is 1e-4.
#'
#' @details The \code{\link{reval_mst}} function evaluates an MST panel by
#' implementing a recursion-based method to assess measurement precision
#' across IRT ability levels. This approach, detailed in Lim et al. (2020),
#' enables the computation of conditional biases and CSEMs efficiently,
#' bypassing the need for extensive simulations traditionally required
#' for MST evaluation.
#'
#' The \code{module} argument, used in conjunction with the item bank metadata
#' \code{x}, systematically organizes items into modules for MST panel evaluation.
#' Each row of \code{x} corresponds to an item, detailing its characteristics
#' like score categories and IRT model. The \code{module} matrix, structured
#' with the same number of rows as \code{x} and columns representing modules,
#' indicates item assignments with 1s. This precise mapping enables
#' the \code{\link{reval_mst}} function to evaluate the MST panel's
#' performance by analyzing how items within each module contribute
#' to measurement precision and bias, reflecting the tailored progression
#' logic inherent in MST designs.
#'
#' The \code{route_map} argument is essential for defining the MST's structure
#' by indicating possible transitions between modules. Similar to the \code{transMatrix}
#' in the \pkg{mstR} package (Magis et al., 2017), \code{route_map} is a binary matrix that outlines
#' which module transitions are possible within an MST design. Each "1" in the matrix
#' represents a feasible transition from one module (row) to another (column),
#' effectively mapping the flow of test takers through the MST based on their
#' performance. For instance, a "1" at the intersection of row \emph{i} and
#' column \emph{j} indicates the possibility for test takers to progress from
#' the module corresponding to row \emph{i} directly to the module denoted
#' by column \emph{j}. This structure allows \code{\link{reval_mst}} to simulate and evaluate
#' the dynamic routing of test takers through various stages and modules of the MST panel.
#'
#' To further detail the \code{cut_score} argument with an illustration:
#' In a 1-3-3 MST configuration, the list \code{cut_score = list(c(-0.5, 0.5), c(-0.6, 0.6))}
#' operates as a decision guide at each stage. Initially, all test takers start
#' in the first module. Upon completion, their scores determine their next stage module:
#' scores below -0.5 route to the first module of the next stage, between -0.5
#' and 0.5 to the second, and above 0.5 to the third. This pattern allows for
#' dynamic adaptation, tailoring the test path to individual performance levels.
#'
#' @return This function returns a list of seven internal objects. The four objects are:
#'
#' \item{panel.info}{A list of several sub-objects containing detailed information
#' about the MST panel configuration, including:
#'     \describe{
#'       \item{config}{A nested list indicating the arrangement of modules across
#'       stages, showing which modules are included in each stage. For example,
#'       the first stage includes module 1, the second stage includes
#'       modules 2 to 4, and so forth.}
#'       \item{pathway}{A matrix detailing all possible pathways through the MST panel.
#'       Each row represents a unique path a test taker might take, based on their
#'       performance and the cut scores defined.}
#'       \item{n.module}{A named vector indicating the number of modules available
#'       at each stage.}
#'       \item{n.stage}{A single numeric value representing the total number of
#'       stages in the MST panel.}
#'    }
#' }
#'
#' \item{item.by.mod}{A list where each entry represents a module
#' in the MST panel, detailing the item metadata within that module. Each module's
#' metadata includes item IDs, the number of categories, the IRT model
#' used (model), and the item parameters (e.g., par.1, par.2, par.3).}
#'
#' \item{item.by.path}{A list containing item metadata arranged
#' according to the paths through the MST structure. This detailed
#' breakdown allows for an analysis of item characteristics along specific
#' MST paths. Each list entry corresponds to a testing stage and path,
#' providing item metadata. This structure facilitates the examination of
#' how items function within the context of each unique path through the MST.}
#'
#' \item{eq.theta}{
#' Estimated ability levels (\eqn{\theta}) corresponding to the observed
#' scores, derived from the inverse TCC scoring method. This provides the
#' estimated \eqn{\theta} values for each potential pathway through the
#' MST stages. For each stage, \eqn{\theta} values are calculated for each
#' path, indicating the range of ability levels across the test takers.
#' For instance, in a three-stage MST, the \code{eq.theta} list may contain
#' \eqn{\theta} estimates for multiple paths within each stage, reflecting
#' the progression of ability estimates as participants move through the test.
#' The example below illustrates the structure of \code{eq.theta} output for a
#' 1-3-3 MST panel with varying paths:
#'      \describe{
#'        \item{stage.1}{\code{path.1} shows \eqn{\theta} estimates ranging from -7 to +7,
#'        demonstrating the initial spread of abilities.}
#'        \item{stage.2}{Multiple paths (\code{path.1}, \code{path.2}, ...) each with
#'        their own \eqn{\theta} estimates, indicating divergence in ability levels
#'        based on test performance.}
#'        \item{stage.3}{Further refinement of \eqn{\theta} estimates across paths,
#'        illustrating the final estimation of abilities after the last stage.}
#'   }
#' }
#'
#' \item{cdist.by.mod}{A list where each entry contains the conditional distributions
#' of the observed scores for each module given the ability levels.}
#'
#' \item{jdist.by.path}{Joint distributions of observed scores for different paths
#' at each stage in a MST panel. The example below outlines the organization of
#' \code{jdist.by.path} data in a hypothetical 1-3-3 MST panel:
#'      \describe{
#'        \item{stage.1}{Represents the distribution at the initial stage, indicating
#'        the broad spread of test-taker abilities at the outset.}
#'        \item{stage.2}{Represents the conditional joint distributions of the observed
#'        scores as test-takers move through different paths at the stage 2, based on
#'        their performance in earlier stages.}
#'        \item{stage.3}{Represents a further refinement of joint distribution of
#'        observed scores as test-takers move through different paths at the final
#'        stage 3, based on their performance in earlier stages.}
#'  }
#' }
#'
#' \item{eval.tb}{
#' A table summarizing the measurement precision of the MST panel. It contains the true
#' ability levels (\code{theta}) with the mean ability estimates (\code{mu}), variance
#' (\code{sigma2}), bias, and conditional standard error of measurement (CSEM) given
#' the true ability levels. This table highlights the MST panel's accuracy and
#' precision across different ability levels, providing insights into its effectiveness
#' in estimating test-taker abilities.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{shape_df}}, \code{\link{est_score}}
#'
#' @references
#' Magis, D., Yan, D., & Von Davier, A. A. (2017). \emph{Computerized adaptive
#' and multistage testing with R: Using packages catR and mstR}. Springer.
#'
#' Lim, H., Davey, T., & Wells, C. S. (2020). A recursion-based analytical
#' approach to evaluate the performance of MST. \emph{Journal of Educational
#' Measurement, 58}(2), 154-178.
#'
#' @examples
#' \donttest{
#' ## ------------------------------------------------------------------------------
#' # Evaluation of a 1-3-3 MST panel using simMST data.
#' # This simulation dataset was utilized in Lim et al.'s (2020) simulation study.
#' # Details:
#' #    (a) Panel configuration: 1-3-3 MST panel
#' #    (b) Test length: 24 items (each module contains 8 items across all stages)
#' #    (c) IRT model: 3-parameter logistic model (3PLM)
#' ## ------------------------------------------------------------------------------
#' # Load the necessary library
#' library(dplyr)
#' library(tidyr)
#' library(ggplot2)
#'
#' # Import item bank metadata
#' x <- simMST$item_bank
#'
#' # Import module information
#' module <- simMST$module
#'
#' # Import routing map
#' route_map <- simMST$route_map
#'
#' # Import cut scores for routing to subsequent modules
#' cut_score <- simMST$cut_score
#'
#' # Import ability levels (theta) for evaluating measurement precision
#' theta <- simMST$theta
#'
#' # Evaluate MST panel using the reval_mst() function
#' eval <-
#'   reval_mst(x,
#'     D = 1.702, route_map = route_map, module = module,
#'     cut_score = cut_score, theta = theta, range.tcc = c(-5, 5)
#'   )
#'
#' # Review evaluation results
#' # The evaluation result table below includes conditional biases and
#' # standard errors of measurement (CSEMs) across ability levels
#' print(eval$eval.tb)
#'
#' # Generate plots for biases and CSEMs
#' p_eval <-
#'   eval$eval.tb %>%
#'   dplyr::select(theta, bias, csem) %>%
#'   tidyr::pivot_longer(
#'     cols = c(bias, csem),
#'     names_to = "criterion", values_to = "value"
#'   ) %>%
#'   ggplot2::ggplot(mapping = ggplot2::aes(x = theta, y = value)) +
#'   ggplot2::geom_point(mapping = ggplot2::aes(shape = criterion), size = 3) +
#'   ggplot2::geom_line(
#'     mapping = ggplot2::aes(
#'       color = criterion,
#'       linetype = criterion
#'     ),
#'     linewidth = 1.5
#'   ) +
#'   ggplot2::labs(x = expression(theta), y = NULL) +
#'   ggplot2::theme_classic() +
#'   ggplot2::theme_bw() +
#'   ggplot2::theme(legend.key.width = unit(1.5, "cm"))
#' print(p_eval)
#' }
#'
#' @export
#' @import dplyr
reval_mst <- function(x, D = 1, route_map, module, cut_score, theta = seq(-5, 5, 1),
                      intpol = TRUE, range.tcc = c(-7, 7), tol = 1e-4) {
  ## -----------------------------------------------------
  # Extract all the panel information from the route map
  ## -----------------------------------------------------
  # Run the panel_info() function
  panel_data <- panel_info(route_map)

  # Pathways allowed
  pathway <- panel_data$pathway

  # Number of modules for each stage
  n.mod <- panel_data$n.module

  # Total number of modules across all stages
  tn.mod <- sum(n.mod)

  # Number of stages
  n.stg <- panel_data$n.stage

  ## -----------------------------------------------------
  # Create List of item metadata for modules and paths
  # across the stages
  ## -----------------------------------------------------
  # List containing item metadata of all modules
  meta_mod <-
    purrr::map(
      .x = 1:tn.mod,
      .f = ~ {
        x[module[, .x] == 1, ] %>%
          tibble::remove_rownames()
      }
    )
  names(meta_mod) <- paste0("m.", 1:tn.mod)

  # List containing item metadata for all possible (sub) pathways at each stage
  meta_path <-
    purrr::map(
      .x = 1:n.stg,
      .f = ~ {
        # Unique pathways upto the current stage
        path_uni <- unique(pathway[, 1:.x, drop = FALSE])

        # Item metadata for the unique pathways
        meta_uni <-
          purrr::map(
            .x = 1:nrow(path_uni),
            .f = ~ {
              meta_mod[path_uni[.x, ]] %>%
                dplyr::bind_rows()
            }
          )

        # Add the names and return the list
        names(meta_uni) <- paste0("path.", 1:length(meta_uni))
        return(meta_uni)
      }
    )
  names(meta_path) <- paste0("stage.", 1:n.stg)

  # Count the cumulative numbers of all possible (sub) pathways by each stage
  cum_path <- purrr::map_dbl(.x = meta_path,
                             .f = ~ {length(.x)})

  ## -----------------------------------------------------
  # Estimate the IRT theta scores using the Inverse-TCC method
  ## -----------------------------------------------------
  # Compute the IRT thetas corresponding to the all observed scores using
  # the Inverse-TCC scoring method for the all (sub) pathways
  eq_theta <-
    purrr::map(
      .x = meta_path,
      .f = ~ {
        # Implement the Inverse-TCC scoring for the current pathway
        tmp_meta <- .x
        tmp_eqtheta <-
          purrr::map(
            .x = tmp_meta,
            .f = ~ {
              inv_tcc_nr(
                x = .x, D = D,
                tol = tol, intpol = intpol,
                range.tcc = range.tcc
              )$est.theta
            }
          ) %>%
          do.call(what = "cbind")
      }
    )

  ## -----------------------------------------------------
  # Compute the conditional distributions of modules at each theta
  ## -----------------------------------------------------
  # List containing the conditional distributions of the modules given the ability levels across all modules.
  # Each element of the list includes the conditional distributions across all ability levels for each module.
  cdist_by_mod <-
    purrr::map(
      .x = meta_mod,
      .f = ~ {
        lwrc(x = .x, theta = theta, D = D)
      }
    )

  # Let's restructure the list above.
  # We need a list of conditional distributions of the modules across all ability levels.
  # Thus, each element of the list includes the conditional distributions for all
  # modules given each ability level.
  cdist_by_th <-
    purrr::map(
      .x = 1:length(theta),
      .f = ~ {
        num <- .x
        tmp_tb <-
          purrr::map(
            .x = cdist_by_mod,
            .f = ~ {
              .x[, num]
            }
          ) %>%
          bind.fill(type = "cbind", fill = 0L)
        colnames(tmp_tb) <- paste0("m.", 1:tn.mod)
        tmp_tb
      }
    )
  names(cdist_by_th) <- theta

  ## -----------------------------------------------------
  # Compute the joint conditional distributions at each theta
  ## -----------------------------------------------------
  # List of possible observed scores across all modules
  score_mod <-
    purrr::map(
      .x = meta_mod[pathway[1, ]],
      .f = ~ {
        0:sum(.x$cats - 1)
      }
    )
  names(score_mod) <- paste0("stage.", 1:n.stg)

  # Cumulative maximum sum scores across the stages
  score_maxcum <- cumsum(sapply(X = score_mod, FUN = "max"))

  # Empty list to include the cut scores used to route a test taker to next module across all stages
  cut4route <- vector("list", n.stg - 1)
  names(cut4route) <- paste0("stage.", 2:n.stg)

  # Empty list to include path labels across all stages
  path_label <- vector("list", n.stg - 1)
  names(path_label) <- paste0("stage.", 2:n.stg)

  # Empty list to include probabilities of the joint distributions across the stages
  # Each element of the list will contain the joint distributions across all ability levels
  joint_dist <- vector("list", n.stg)
  names(joint_dist) <- paste0("stage.", 1:n.stg)

  # For the routing module in the first stage,
  # assign the conditional distributions of the routing modules across the ability levels
  joint_dist[[1]] <-
    purrr::map(
      .x = cdist_by_th,
      .f = ~ {
        .x[, 1, drop = FALSE] %>%
          unname()
      }
    )
  names(joint_dist[[1]]) <- round(theta, 10)

  # Calculate the probabilities of the conditional joint distributions
  # given the ability levels across all stages
  for (s in 2:n.stg) {
    # The numbers of all possible (sub) pathways at the current stage
    n_path <- cum_path[s - 1]

    # Empty list to contain the cut scores used to route a test taker to next module
    cut4route[[s - 1]] <- vector("list", n_path)

    # A data.frame with the pathway by the current stage (first column)
    # and the modules (the current stage (second column)
    path_pre <-
      data.frame(pathway[, 1:s]) %>%
      tidyr::unite(col = "path", 1:(s - 1), sep = "_")

    # Unique pathways by the current stage
    uni_path_pre <- unique(path_pre$path)

    # Empty to contain the conditional joint distributions at the current stage
    jdist_temp <- vector("list", n_path)

    # Loop over the all possible (sub) pathways
    for (i in 1:n_path) {
      # Find cut scores for routing a next module
      # Firstly, find all unique next modules that will be assigned at the current pathway
      nmod_next <-
        path_pre %>%
        dplyr::filter(.data$path %in% uni_path_pre[i]) %>%
        dplyr::pull(2) %>%
        unique()

      # Secondly, find the indices which modules will be actually used among the all possible next modules
      idx_next <- match(nmod_next, unique(pathway[, s]))

      # Lastly, find the cut scores to be used to assign the next modules
      cut4route[[s - 1]][[i]] <- cut_score[[s - 1]][dplyr::lag(idx_next)[-1]]

      # Find the path labels to determine which module will be assigned to each
      # ability level based on the cut scores
      path_label[[s - 1]][[i]] <-
        give_path(
          score = unlist(eq_theta[[s - 1]][, i]),
          cut_sc = cut4route[[s - 1]][[i]]
        )$path

      # Cumulative maximum sum scores by the current stages
      max_score <- score_maxcum[s]

      # Number of all unique next modules (a.k.a. all unique possible routing numbers)
      n_route <- length(nmod_next)

      # All possible observed scores by the current stage
      score_cur <- 0:score_maxcum[s - 1]

      # Possible observed score for the next module
      score_next <- score_mod[[s]]

      # The conditional joint distributions given the ability levels by the current stage
      jdist_cur <-
        purrr::map(
          .x = joint_dist[[s - 1]],
          .f = ~ {
            .x[, i, drop = FALSE]
          }
        )

      # The conditional distribution of the next modules given the ability levels
      cdist_next <-
        purrr::map(
          .x = cdist_by_th,
          .f = ~ {
            .x[, nmod_next, drop = FALSE]
          }
        )

      # Column bind of the above two distributions given the ability levels
      cbind_dist <-
        purrr::map(
          .x = 1:length(theta),
          .f = ~ {
            bind.fill(
              List = list(jdist_cur[[.x]], cdist_next[[.x]]),
              type = "cbind",
              fill = 0L
            )
          }
        )

      # Compute the joint conditional distributions
      jdist_temp[[i]] <-
        jdist(
          cbind_dist = cbind_dist, path_lb = path_label[[s - 1]][[i]],
          n_route = n_route, max_score = max_score, score_cur = score_cur,
          score_next = score_next, theta = theta
        )
    }

    # Add to the list
    joint_dist[[s]] <-
      purrr::map(
        .x = 1:length(theta),
        .f = ~ {
          do.call("cbind", lapply(jdist_temp, "[[", .x))
        }
      )
    names(joint_dist[[s]]) <- round(theta, 10)
  }

  ## -----------------------------------------------------
  # Evaluate the measurement precision of the MST panel
  ## -----------------------------------------------------
  # Calculate the means and variances of the conditional joint distributions
  # across all ability levels
  cmoment <-
    lapply(
      X = joint_dist[[s]],
      FUN = cal_moment,
      node = eq_theta[[s]]
    ) %>%
    do.call(what = "rbind")

  # Create a table containing the moments and
  # conditional biases and standard errors of measurement (SEMs)
  eval_tb <-
    data.frame(theta = theta, cmoment) %>%
    dplyr::mutate(
      bias = .data$mu - .data$theta,
      csem = sqrt(.data$sigma2)
    ) %>%
    tibble::remove_rownames()


  ## -----------------------------------------------------
  # Return the results
  ## -----------------------------------------------------
  rst <- list(
    panel.info = panel_data,
    item.by.mod = meta_mod,
    item.by.path = meta_path,
    eq.theta = eq_theta,
    cdist.by.mod = cdist_by_th,
    jdist.by.path = joint_dist,
    eval.tb = eval_tb
  )
  rst
}


# This function computes the conditional joint distributions of the current pathways
# given the ability levels
jdist <- function(cbind_dist, path_lb, n_route, max_score, score_cur,
                  score_next, theta) {
  # Empty list across all ability levels
  cj_dist_all <- vector("list", length(theta))

  # Loop over the ability levels
  for (j in 1:length(theta)) {
    # Matrix to contain the conditional joint distribution at a given ability level
    cj_dist <- array(0, c(max_score + 1, n_route))
    rownames(cj_dist) <- 0:max_score

    # Loop over the observed scores of the (sub) pathway by the the current stage
    for (i in 1:length(score_cur)) {
      # Compute the conditional joint distribution
      cj_dist[i:(i + max(score_next)), path_lb[i]] <-
        cj_dist[i:(i + max(score_next)), path_lb[i]] +
        (cbind_dist[[j]][1:length(score_next), (path_lb[i] + 1)] * cbind_dist[[j]][i, 1])
    }

    # Add to the list
    cj_dist_all[[j]] <- cj_dist
  }

  # Return the results
  cj_dist_all
}

# This function determines the next path that each test taker with
# a given ability level needs to take
give_path <- function(score, cut_sc) {
  # Number of theta scores
  nscore <- length(score)

  # Number of categories (or labels) to which each theta will be given
  ncats <- length(cut_sc) + 1

  # A vector of labels to be given to each test taker
  labels <- 1:ncats

  # Create pathway vectors
  path <-
    cut(x = score, breaks = c(-Inf, cut_sc, Inf), labels = labels) %>%
    as.numeric()

  # Return the results
  rst <- list(path = path, score = score, cut.sc = cut_sc)
  rst
}

# This function extracts the panel information (e.g., available pathways,
# the number of modules, stages, etc) given the route map information
panel_info <- function(route_map) {
  # Transform the format of the route_map to data.frame
  route_map <- as.data.frame(route_map)

  ## Count the total number of modules
  end_col <- ncol(route_map)

  # Find the routing module (a.k.a the module in the first stage)
  mod_stg1 <- which(colSums(route_map) == 0)

  # Find the modules in the subsequent stages
  config <- list()
  config[[1]] <- mod_stg1
  col_bef <- mod_stg1
  i <- 1
  repeat{
    i <- i + 1
    col_aft <- which(colSums(route_map[col_bef, ] == 1) > 0)
    config[[i]] <- col_aft
    if (max(col_aft) == end_col) break
    col_bef <- col_aft
  }
  names(config) <- paste0("stage.", 1:length(config))
  config <- lapply(X = config, FUN = "unname")

  # Count the total number of stages and the number of modules at each stage
  n.mod <- sapply(X = config, FUN = "length")
  n.stg <- length(n.mod)

  # Create a data.frame including all possible pathways with no restriction
  path_all <- expand.grid(config)

  # Delete the pathways that are not allowed to used
  remain <- path_all
  for (s in 1:(n.stg - 1)) {
    for (m in 1:n.mod[s]) {
      tmp.1 <- config[[s]][m]
      tmp.2 <- which(route_map[tmp.1, ] == 1)
      delete <- (remain[, s] == tmp.1 & !(remain[, s + 1] %in% tmp.2))
      remain <- remain[!delete, ]
    }
  }

  # Arrange the pathway matrix
  remain <-
    dplyr::arrange(remain, dplyr::pick(dplyr::everything())) %>%
    as.matrix()
  rownames(remain) <- paste0("path.", 1:nrow(remain))

  # Return the results
  rst <- list(config = config, pathway = remain, n.module = n.mod, n.stage = n.stg)
  rst
}
