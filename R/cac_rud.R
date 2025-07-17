#' Classification Accuracy and Consistency Based on Rudner's (2001, 2005)
#' Approach
#'
#' This function computes classification accuracy and consistency indices using
#' the method proposed by Rudner in 2001 and 2005. This function supports both
#' scenarios: when the empirical ability distribution of the population is
#' available, and when individual ability estimates are used.
#'
#' @inheritParams cac_lee
#' @param se A numeric vector of the same length as `theta` representing the
#'   standard errors associated with each ability estimate. See the **Details**
#'   section for more information
#'
#' @details This function first validates the input arguments. If both `theta`
#' and `weights` are `NULL`, the function will stop and return an error message.
#' Either `theta` or `weights` must be specified. In addition, `se` must be
#' provided and must match the length of `theta` or the number of quadrature
#' points in `weights`.
#'
#' It then computes the probability that an examinee with a given ability is
#' classified into each performance level using the normal distribution function
#' centered at each `theta` (or quadrature point) with standard deviation `se`.
#' These probabilities are used to calculate conditional classification accuracy
#' (the probability of being correctly classified) and conditional classification
#' consistency (the probability of being consistently classified upon repeated
#' testing) for each ability value.
#'
#' Finally, the function computes marginal classification accuracy and
#' consistency across all examinees by aggregating the conditional indices with
#' the associated weights.
#'
#' @return A list containing the following elements:
#'  - confusion: A confusion matrix showing the cross table between true and expected levels.
#'  - marginal: A data frame showing the marginal classification accuracy and consistency indices.
#'  - conditional: A data frame showing the conditional classification accuracy and consistency indices.
#'  - prob.level: A data frame showing the probability of being assigned to each level category.
#'  - cutscore: A numeric vector showing the cut scores used in the analysis.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::gen.weight()], [irtQ::est_score()], [irtQ::cac_lee()]
#'
#' @references Rudner, L. M. (2001). Computing the expected proportions of
#'   misclassified examinees.
#' *Practical Assessment, Research, and Evaluation, 7*(1), 14.
#'
#'   Rudner, L. M. (2005). Expected classification accuracy. *Practical
#'   Assessment, Research, and Evaluation, 10*(1), 13.
#'
#' @examples
#' \donttest{
#' ## -------------------------------------------
#' # 1. Using the empirical ability distribution
#' ## -------------------------------------------
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameter estimates and convert them into item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Define cut scores on the theta scale
#' cutscore <- c(-2, -0.5, 0.8)
#'
#' # Create quadrature points and corresponding weights
#' node <- seq(-4, 4, 0.25)
#' weights <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = node)
#'
#' # Compute conditional standard errors across quadrature points
#' tif <- info(x = x, theta = node, D = 1, tif = TRUE)$tif
#' se <- 1 / sqrt(tif)
#'
#' # Compute classification accuracy and consistency
#' cac_1 <- cac_rud(cutscore = cutscore, se = se, weights = weights)
#' print(cac_1)
#'
#' ## -----------------------------------------
#' # 2. Using individual ability estimates
#' ## -----------------------------------------
#'
#' # Generate true abilities from N(0, 1)
#' set.seed(12)
#' theta <- rnorm(n = 1000, mean = 0, sd = 1)
#'
#' # Simulate item response data
#' data <- simdat(x = x, theta = theta, D = 1)
#'
#' # Estimate ability and standard errors using ML estimation
#' est_theta <- est_score(
#'   x = x, data = data, D = 1, method = "ML",
#'   range = c(-4, 4), se = TRUE
#' )
#' theta_hat <- est_theta$est.theta
#' se <- est_theta$se.theta
#'
#' # Compute classification accuracy and consistency
#' cac_2 <- cac_rud(cutscore = cutscore, theta = theta_hat, se = se)
#' print(cac_2)
#' }
#'
#' @import dplyr
#' @export
cac_rud <- function(cutscore,
                    theta = NULL,
                    se,
                    weights = NULL) {

  # check if the provided inputs are correct
  if (missing(se)) {
    stop("The standard errors must be provided in `se` argument.", call. = FALSE)
  }

  # check if the provided inputs are correct
  if (is.null(theta) & is.null(weights)) {
    stop("Eighter of `theta` or `weights` argument must not be NULL; both cannot be NULL",
         call. = FALSE
    )
  }

  # count the number of levels
  n.lev <- length(cutscore) + 1

  # add the bounds to the cut scores
  breaks <- c(-Inf, cutscore, Inf)

  # (1) when the quadrature points and the corresponding ses are provided
  if (!is.null(weights)) {
    # extract nodes and weights
    nodes <- weights[, 1]
    wts <- weights[, 2]

    # count the number of thetas
    n.theta <- length(wts)

    # check if the provided inputs are correct
    if (n.theta != length(se)) {
      stop("The numbers of weights and the standard errors must be equal.",
           call. = FALSE
      )
    }

    # assign the levels to each theta
    level <-
      cut(
        x = nodes, breaks = breaks, labels = FALSE,
        include.lowest = TRUE, right = FALSE, dig.lab = 7
      )

    # create an empty data frame to contain the conditional
    # classification accuracy and consistency for each theta value
    cond_tb <- data.frame(
      theta = nodes,
      weights = wts,
      level = level,
      accuracy = NA_real_,
      consistency = NA_real_
    )

    # create an empty matrix to contain the probability
    # that each examinee with a specific ability is assigned
    # to each level category
    ps_tb <- matrix(NA, nrow = n.theta, ncol = n.lev)
    colnames(ps_tb) <- paste0("p.level.", 1:n.lev)

    # compute the probabilities that is assigned to an examinee
    # with a specific ability
    for (i in 1:n.theta) {
      # compute the cumulative probability across all cut scores
      cum_ps <- stats::pnorm(q = breaks, mean = nodes[i], sd = se[i])

      # the probability that each examinee will be assigned to each level category
      ps <- diff(cum_ps)
      ps_tb[i, ] <- ps

      # conditional classification accuracy
      cond_ca <- ps[level[i]]
      cond_tb[i, 4] <- cond_ca

      # conditional classification consistency
      cond_cc <- sum(ps^2)
      cond_tb[i, 5] <- cond_cc
    }

    # compute the marginal accuracy and consistency
    margin_tb <-
      cond_tb %>%
      dplyr::group_by(.data$level, .drop = FALSE) %>%
      dplyr::summarise(
        accuracy = sum(.data$accuracy * .data$weights),
        consistency = sum(.data$consistency * .data$weights),
        .groups = "drop"
      ) %>%
      janitor::adorn_totals(where = "row", name = "marginal")

    # add more variables to ps_tb
    ps_tb2 <-
      data.frame(
        theta = nodes, weights = wts,
        level = level, ps_tb
      )

    # create a cross table between true and expected levels
    cross_tb <-
      ps_tb2 %>%
      dplyr::group_by(.data$level, .drop = FALSE) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with("p.level."),
          ~ {
            sum(.x * .data$weights)
          }
        ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$level) %>%
      tibble::column_to_rownames("level") %>%
      data.matrix()
    dimnames(cross_tb) <- list(True = 1:n.lev, Expected = 1:n.lev)
  } else {
    # (2) when individual ability estimates and ses are provided
    # count the number of thetas
    n.theta <- length(theta)

    # check if the provided inputs are correct
    if (n.theta != length(se)) {
      stop("The numbers of thetas and the standard errors must be equal.",
           call. = FALSE
      )
    }

    # assign uniform weights
    wts <- 1 / n.theta

    # assign the levels to each theta
    level <-
      cut(
        x = theta, breaks = breaks, labels = FALSE,
        include.lowest = TRUE, right = FALSE, dig.lab = 7
      )

    # create an empty data frame to contain the conditional
    # classification accuracy and consistency for each theta value
    cond_tb <- data.frame(
      theta = theta,
      weights = wts,
      level = level,
      accuracy = NA_real_,
      consistency = NA_real_
    )

    # create an empty matrix to contain the probability
    # that each examinee with a specific ability is assigned
    # to each level category
    ps_tb <- matrix(NA, nrow = n.theta, ncol = n.lev)
    colnames(ps_tb) <- paste0("p.level.", 1:n.lev)

    # compute the probabilities that is assigned to an examinee
    # with a specific ability
    for (i in 1:n.theta) {
      # compute the cumulative probability across all cut scores
      cum_ps <- stats::pnorm(q = breaks, mean = theta[i], sd = se[i])

      # the probability that each examinee will be assigned to each level category
      ps <- diff(cum_ps)
      ps_tb[i, ] <- ps

      # conditional classification accuracy
      cond_ca <- ps[level[i]]
      cond_tb[i, 4] <- cond_ca

      # conditional classification consistency
      cond_cc <- sum(ps^2)
      cond_tb[i, 5] <- cond_cc
    }

    # compute the marginal accuracy and consistency
    margin_tb <-
      cond_tb %>%
      dplyr::group_by(.data$level, .drop = FALSE) %>%
      dplyr::summarise(
        accuracy = sum(.data$accuracy * .data$weights),
        consistency = sum(.data$consistency * .data$weights),
        .groups = "drop"
      ) %>%
      janitor::adorn_totals(where = "row", name = "marginal")

    # add more variables to ps_tb
    ps_tb2 <-
      data.frame(
        theta = theta, weights = wts,
        level = level, ps_tb
      )

    # create a cross table between true and expected levels
    cross_tb <-
      ps_tb2 %>%
      dplyr::group_by(.data$level, .drop = FALSE) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::starts_with("p.level."),
          ~ {
            sum(.x * .data$weights)
          }
        ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$level) %>%
      tibble::column_to_rownames("level") %>%
      data.matrix()
    dimnames(cross_tb) <- list(True = 1:n.lev, Expected = 1:n.lev)
  }

  # return the results
  rst <- list(
    confusion = round(cross_tb, 7),
    marginal = margin_tb,
    conditional = cond_tb,
    prob.level = ps_tb2,
    cutscore = cutscore
  )
  return(rst)
}
