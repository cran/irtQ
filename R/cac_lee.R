#' Classification Accuracy and Consistency Using Lee's (2010) Approach
#'
#' This function computes classification accuracy and consistency indices for
#' complex assessments based on the method proposed by Lee (2010). This function
#' supports both dichotomous and polytomous item response theory (IRT) models.
#'
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function.
#' @param cutscore A numeric vector specifying the cut scores for
#'   classification. Cut scores are the points that separate different
#'   performance categories (e.g., pass vs. fail, or different grades).
#' @param theta A numeric vector of ability estimates. Ability estimates (theta
#'   values) are the individual proficiency estimates obtained from the IRT
#'   model. The theta parameter is optional and can be `NULL`.
#' @param weights An optional two-column data frame or matrix where the first
#'   column is the quadrature points (nodes) and the second column is the
#'   corresponding weights. This is typically used in quadrature-based IRT
#'   analysis.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.7 is commonly
#'   used for this purpose. Default is 1.
#' @param cut.obs Logical. If `TRUE`, it indicates the cutscores on the
#'   observed-summed score metric. If `FALSE`, it indicates they are on the IRT
#'   theta score metric. Default is `TRUE`.
#'
#' @details
#' This function first validates the input arguments. If both `theta` and `weights`
#' are `NULL`, the function will stop and return an error message. Either `theta`
#' (a vector of ability estimates) or `weights` (a quadrature-based weight matrix)
#' must be specified.
#'
#' If `cut.obs = FALSE`, the provided cut scores are assumed to be on the theta (ability)
#' scale, and they are internally converted to the observed summed score scale using the
#' test characteristic curve (TCC). This transformation allows classification to be carried
#' out on the summed score metric, even if theta-based cut points are provided.
#'
#' When `weights` are provided (D method), the function uses the Lord-Wingersky recursive
#' algorithm to compute the conditional distribution of observed total scores at each node.
#' These conditional distributions are used to compute:
#' - the probability of being classified into each performance level,
#' - conditional classification accuracy (probability of correct classification), and
#' - conditional classification consistency (probability of being assigned to the same level
#'   upon repeated testing).
#'
#' When `theta` values are provided instead (P method), the same logic applies, but using
#' an empirical distribution of examinees instead of quadrature-based integration.
#' In this case, uniform weights are assigned to all examinees.
#'
#' Finally, marginal classification accuracy and consistency are computed as weighted
#' averages of the conditional statistics across the ability distribution.
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
#' @seealso [irtQ::gen.weight()], [irtQ::est_score()], [irtQ::cac_rud()]
#'
#' @references Lee, W. C. (2010). Classification consistency and accuracy for
#' complex assessments using item response theory. *Journal of Educational
#' Measurement, 47*(1), 1-17.
#'
#' @examples
#' \donttest{
#' ## --------------------------------------------------------------------------
#' ## 1. When the empirical ability distribution of the population is available
#' ##    (D method)
#' ## --------------------------------------------------------------------------
#'
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameter data and convert it to item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Set the cut scores on the observed summed score scale
#' cutscore <- c(10, 20, 30, 50)
#'
#' # Create a data frame containing the quadrature points and corresponding weights
#' node <- seq(-4, 4, 0.25)
#' weights <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = node)
#'
#' # Calculate classification accuracy and consistency
#' cac_1 <- cac_lee(x = x, cutscore = cutscore, weights = weights, D = 1)
#' print(cac_1)
#'
#' ## -------------------------------------------------------------
#' ## 2. When individual ability estimates are available (P method)
#' ## -------------------------------------------------------------
#'
#' # Randomly draw true ability values from N(0, 1)
#' set.seed(12)
#' theta <- rnorm(n = 1000, mean = 0, sd = 1)
#'
#' # Simulate item response data
#' data <- simdat(x = x, theta = theta, D = 1)
#'
#' # Estimate ability parameters using maximum likelihood (ML)
#' est_th <- est_score(
#'   x = x, data = data, D = 1, method = "ML",
#'   range = c(-4, 4), se = FALSE
#' )$est.theta
#'
#' # Calculate classification accuracy and consistency
#' cac_2 <- cac_lee(x = x, cutscore = cutscore, theta = est_th, D = 1)
#' print(cac_2)
#'
#' ## ---------------------------------------------------------
#' ## 3. When individual ability estimates are available,
#' ##    but cut scores are specified on the IRT theta scale
#' ## ---------------------------------------------------------
#' # Set the cut scores on the theta scale
#' cutscore <- c(-2, -0.4, 0.2, 1.0)
#'
#' # Calculate classification accuracy and consistency
#' cac_3 <- cac_lee(
#'   x = x, cutscore = cutscore, theta = est_th, D = 1,
#'   cut.obs = FALSE
#' )
#' print(cac_3)
#' }
#'
#' @import dplyr
#' @export
cac_lee <- function(x,
                    cutscore,
                    theta = NULL,
                    weights = NULL,
                    D = 1,
                    cut.obs = TRUE) {

  # check if the provided inputs are correct
  if (is.null(theta) & is.null(weights)) {
    stop("Eighter of `theta` or `weights` argument must not be NULL; both cannot be NULL",
         call. = FALSE
    )
  }

  # if the cutscores are on the theta metric, compute the expected cutscores
  # on the observed score metric
  if (!cut.obs) {
    cutscore <- irtQ::traceline(x = x, theta = cutscore, D = D)$tcc
  }

  # count the number of levels
  n.lev <- length(cutscore) + 1

  # (1) when the quadrature points and the corresponding ses are provided
  if (!is.null(weights)) {
    # extract nodes and weights
    nodes <- weights[, 1]
    wts <- weights[, 2]

    # count the number of thetas
    n.theta <- length(wts)

    # estimate likelihood functions using lord-wingersky algorithm
    # for each theta, which is the conditional observed score distribution
    # given each theta; each column is the conditional score distribution for
    # a given theta
    lkhd <- lwrc(x = x, theta = weights[, 1], prob = NULL, D = D)

    # count the total number of observed scores for the test
    n.score <- nrow(lkhd)

    # check the maximum possible observed score
    max.score <- n.score - 1

    # add the bounds to the cut scores
    breaks <- c(0, cutscore, n.score)

    # compute the expected true score using the IRT models for each theta
    tscore <- traceline(x = x, theta = nodes, D = D)$tcc

    # assign the levels to each expected true score
    level <-
      cut(
        x = tscore, breaks = breaks, labels = 1:n.lev,
        include.lowest = TRUE, right = FALSE, dig.lab = 7
      )

    # create an empty data frame to contain the conditional
    # classification accuracy and consistency for each theta value
    cond_tb <- data.frame(
      theta = nodes,
      weights = wts,
      true.score = tscore,
      level = level,
      accuracy = NA_real_,
      consistency = NA_real_
    )

    # create an empty matrix to contain the probability
    # that each examinee with a specific ability is assigned
    # to each level category
    ps_tb <- matrix(NA, nrow = n.theta, ncol = n.lev)
    colnames(ps_tb) <- paste0("p.level.", 1:n.lev)

    # assign the level to each possible observed score
    loc.lev <-
      cut(
        x = 0:max.score, breaks = breaks, labels = 1:n.lev,
        include.lowest = TRUE, right = FALSE, dig.lab = 7
      )

    # compute the probabilities that is assigned to an examinee
    # with a specific ability
    for (i in 1:n.theta) {
      # the probability that each examinee will be assigned to each level category
      ps <-
        data.frame(level = loc.lev, prob = lkhd[, i]) %>%
        dplyr::summarise(prob = sum(.data$prob), .by = "level") %>%
        dplyr::pull("prob")
      ps_tb[i, ] <- ps

      # conditional classification accuracy
      cond_ca <- ps[level[i]]
      cond_tb[i, 5] <- cond_ca

      # conditional classification consistency
      cond_cc <- sum(ps^2)
      cond_tb[i, 6] <- cond_cc
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
        true.score = tscore, level = level,
        ps_tb
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

    # assign uniform weights
    wts <- 1 / n.theta

    # estimate likelihood functions using lord-wingersky algorithm
    # for each theta, which is the conditional observed score distribution
    # given each theta; each column is the conditional score distributon for
    # a given theta
    lkhd <- lwrc(x = x, theta = theta, prob = NULL, D = D)

    # count the total number of observed scores for the test
    n.score <- nrow(lkhd)

    # check the maximum possible observed score
    max.score <- n.score - 1

    # add the bounds to the cut scores
    breaks <- c(0, cutscore, n.score)

    # compute the expected true score using the IRT models for each theta
    tscore <- traceline(x = x, theta = theta, D = D)$tcc

    # assign the levels to each expected true score
    level <-
      cut(
        x = tscore, breaks = breaks, labels = 1:n.lev,
        include.lowest = TRUE, right = FALSE, dig.lab = 7
      )

    # create an empty data frame to contain the conditional
    # classification accuracy and consistency for each theta value
    cond_tb <- data.frame(
      theta = theta,
      weights = wts,
      true.score = tscore,
      level = level,
      accuracy = NA_real_,
      consistency = NA_real_
    )

    # create an empty matrix to contain the probability
    # that each examinee with a specific ability is assigned
    # to each level category
    ps_tb <- matrix(NA, nrow = n.theta, ncol = n.lev)
    colnames(ps_tb) <- paste0("p.level.", 1:n.lev)

    # assign the level to each possible observed score
    loc.lev <-
      cut(
        x = 0:max.score, breaks = breaks, labels = 1:n.lev,
        include.lowest = TRUE, right = FALSE, dig.lab = 7
      )

    # compute the probabilities that is assigned to an examinee
    # with a specific ability
    for (i in 1:n.theta) {
      # the probability that each examinee will be assigned to each level category
      ps <-
        data.frame(level = loc.lev, prob = lkhd[, i]) %>%
        dplyr::summarise(prob = sum(.data$prob), .by = "level") %>%
        dplyr::pull("prob")
      ps_tb[i, ] <- ps

      # conditional classification accuracy
      cond_ca <- ps[level[i]]
      cond_tb[i, 5] <- cond_ca

      # conditional classification consistency
      cond_cc <- sum(ps^2)
      cond_tb[i, 6] <- cond_cc
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
        true.score = tscore, level = level,
        ps_tb
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
