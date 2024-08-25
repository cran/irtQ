#' Classification accuracy and consistency using Lee's (2010) approach.
#'
#' @description
#' This function computes the classification accuracy and consistency indices
#' for complex assessments using the method proposed by Lee (2010). The function
#' can handle both dichotomous and polytomous item response theory (IRT) models.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of
#' score categories, models ...). See \code{\link{est_irt}}, \code{\link{irtfit}}, \code{\link{info}}
#' or \code{\link{simdat}} for more detail about the item metadata.
#' @param cutscore A numeric vector specifying the cut scores for classification.
#' Cut scores are the points that separate different performance categories
#' (e.g., pass vs. fail, or different grades).
#' @param theta A numeric vector of ability estimates. Ability estimates (theta values)
#' are the individual proficiency estimates obtained from the IRT model. The theta
#' parameter is optional and can be NULL.
#' @param weights An optional two-column data frame or matrix where the first column
#' is the quadrature points (nodes) and the second column is the corresponding weights.
#' This is typically used in quadrature-based IRT analysis.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible
#' to the normal ogive function (if set to 1.7). Default is 1.
#' @param cut.obs A logical value. If TRUE, it indicates the cutscores on the observed-summed score
#' metric. If FALSE, it indicates they are on the IRT theta score metric. Default is TRUE.
#'
#' @details
#' The function works by first checking the provided inputs. If both theta and weights are NULL,
#' the function will stop and return an error message. Depending on the provided inputs, the function
#' will then compute the classification accuracy and consistency indices using either the quadrature
#' points and corresponding weights (D method) or individual ability estimates (P method). The function
#' returns a list containing the confusion matrix, marginal and conditional classification accuracy and
#' consistency indices, probabilities of being assigned to each level category, and cut scores.
#'
#' @return
#' A list containing the following elements:
#' \itemize{
#' \item confusion: A confusion matrix showing the cross table between true and expected levels.
#' \item marginal: A data frame showing the marginal classification accuracy and consistency indices.
#' \item conditional: A data frame showing the conditional classification accuracy and consistency indices.
#' \item prob.level: A data frame showing the probability of being assigned to each level category.
#' \item cutscore: A numeric vector showing the cut scores used in the analysis.
#' }
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{gen.weight}}, \code{\link{est_score}}, \code{\link{cac_rud}}
#'
#' @references
#' Lee, W. C. (2010). Classification consistency and accuracy for complex assessments
#' using item response theory. \emph{Journal of Educational Measurement, 47}(1), 1-17.
#'
#' @examples
#' \donttest{
#' ## ------------------------------------------------------------------------------
#' # 1. When the empirical ability distribution of the population is available
#' #    (D method)
#' ## ------------------------------------------------------------------------------
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # read item parameter data and transform it to item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # set the cut scores on the observed-summed score metric
#' cutscore <- c(10, 20, 30, 50)
#'
#' # create the data frame including the quadrature points
#' # and the corresponding weights
#' node <- seq(-4, 4, 0.25)
#' weights <- gen.weight(dist = "norm", mu = 0, sigma = 1, theta = node)
#'
#' # calculate the classification accuracy and consistency
#' cac_1 <- cac_lee(x = x, cutscore = cutscore, weights = weights, D = 1)
#' print(cac_1)
#'
#' ## ------------------------------------------------------------------------------
#' # 2. When individual ability estimates are available (P method)
#' ## ------------------------------------------------------------------------------
#' # randomly select the true abilities from N(0, 1)
#' set.seed(12)
#' theta <- rnorm(n = 1000, mean = 0, sd = 1)
#'
#' # simulate the item response data
#' data <- simdat(x = x, theta = theta, D = 1)
#'
#' # estimate the ability parameters using the ML estimation
#' est_th <- est_score(
#'   x = x, data = data, D = 1, method = "ML",
#'   range = c(-4, 4), se = FALSE
#' )$est.theta
#'
#' # calculate the classification accuracy and consistency
#' cac_2 <- cac_lee(x = x, cutscore = cutscore, theta = est_th, D = 1)
#' print(cac_2)
#'
#' ## ------------------------------------------------------------------------------
#' # 3. When individual ability estimates are available,
#' #    but, the cutscores are on the IRT theta score metric
#' ## ------------------------------------------------------------------------------
#' # set the theta cut scures
#' cutscore <- c(-2, -0.4, 0.2, 1.0)
#'
#' # calculate the classification accuracy and consistency
#' cac_3 <- cac_lee(
#'   x = x, cutscore = cutscore, theta = est_th, D = 1,
#'   cut.obs = FALSE
#' )
#' print(cac_3)
#' }
#'
#' @import dplyr
#' @export
cac_lee <- function(x, cutscore, theta = NULL, weights = NULL, D = 1,
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
