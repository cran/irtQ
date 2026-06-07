#' Lord-Wingersky Recursion Formula
#'
#' This function computes the conditional distributions of number-correct (or
#' observed) scores given either the probabilities of category responses for
#' each item or a set of theta values, using the Lord and Wingersky recursion
#' formula (1984).
#'
#' @param x A data frame containing item metadata (e.g., item parameters, number
#'   of categories, IRT model types, etc.). See [irtQ::est_irt()] or
#'   [irtQ::simdat()] for more details about the item metadata. This data frame
#'   can be easily created using the [irtQ::shape_df()] function. If `prob` is
#'   `NULL`, the item metadata in `x` is used in the recursion formula. See
#'   **Details** below.
#' @param theta A vector of theta values at which the conditional distributions
#'   of observed scores are computed. This argument is required only when item
#'   metadata is provided via the `x` argument.
#' @param prob A matrix containing the category response probabilities for each
#'   item. Each row corresponds to an item, and each column represents a score
#'   category. If items have different numbers of categories, empty cells should
#'   be filled with zeros or `NA` values. If `x` is `NULL`, this matrix is used
#'   in the recursion formula.
#' @param cats A numeric vector specifying the number of score categories for
#'   each item. For example, a dichotomous item has two categories. This
#'   argument is required only when a probability matrix is provided via the
#'   `prob` argument.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.7 is commonly
#'   used for this purpose. Default is 1.
#'
#' @details The Lord and Wingersky recursive algorithm provides an efficient
#'   method for calculating the compound probabilities of all possible
#'   number-correct (i.e., observed) scores on a test, based on IRT models. This
#'   algorithm is particularly useful for obtaining the IRT model-based
#'   distribution of observed scores.
#'
#'   The conditional distributions of observed scores can be computed using
#'   either the item metadata specified in `x` or the category probability
#'   matrix specified in `prob`.
#'
#' @return When the `prob` argument is provided, the function returns a vector
#'   of probabilities for all possible observed scores on a test.
#'
#'   When the `x` argument is specified, it returns a matrix of conditional
#'   probabilities for each observed score across all specified theta values.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Kolen, M. J. & Brennan, R. L. (2004) *Test Equating, Scaling, and
#'   Linking* (2nd ed.). New York: Springer.
#'
#'   Lord, F. & Wingersky, M. (1984). Comparison of IRT true score and
#'   equipercentile observed score equatings. *Applied Psychological Measurement,
#'   8*(4), 453-461.
#'
#' @examples
#' ## Example 1: Using a matrix of category probabilities
#' ## This example is from Kolen and Brennan (2004, p. 183)
#' # Create a matrix of probabilities for correct and incorrect responses to three items
#' probs <- matrix(c(.74, .73, .82, .26, .27, .18), nrow = 3, ncol = 2, byrow = FALSE)
#'
#' # Create a vector specifying the number of score categories for each item
#' cats <- c(2, 2, 2)
#'
#' # Compute the conditional distribution of observed scores
#' lwrc(prob = probs, cats = cats)
#'
#' ## Example 2: Using a matrix of category probabilities for a mixed-format test
#' # Category probabilities for a dichotomous item
#' p1 <- c(0.2, 0.8, 0, 0, 0)
#'
#' # Category probabilities for another dichotomous item
#' p2 <- c(0.4, 0.6, NA, NA, NA)
#'
#' # Category probabilities for a polytomous item with five categories
#' p3 <- c(0.1, 0.2, 0.2, 0.4, 0.1)
#'
#' # Category probabilities for a polytomous item with three categories
#' p4 <- c(0.5, 0.3, 0.2, NA, NA)
#'
#' # Combine the probability vectors into a matrix
#' p <- rbind(p1, p2, p3, p4)
#'
#' # Create a vector specifying the number of score categories for each item
#' cats <- c(2, 2, 5, 3)
#'
#' # Compute the conditional distribution of observed scores
#' lwrc(prob = p, cats = cats)
#'
#' ## Example 3: Using a data frame of item metadata for a mixed-format test
#' # Import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Read item parameters and convert them to item metadata
#' x <- bring.flexmirt(file = flex_prm, "par")$Group1$full_df
#'
#' # Compute the conditional distribution of observed scores for a range of theta values
#' lwrc(x = x, theta = seq(-4, 4, 0.2), D = 1)
#'
#' @export
lwrc <- function(x = NULL, theta, prob = NULL, cats, D = 1) {
  if (is.null(x)) {
    max.cat <- ncol(prob)
    n.cats <- length(cats)
    n.prob <- nrow(prob)

    logic <- is.numeric(prob)
    if (!logic) stop("All numbers in 'prob' should be numeric.", call. = FALSE)

    if (n.prob != n.cats) {
      stop(paste0(
        "There are ", n.prob, " items in the probability matrix (or data.frame), ",
        "whereas there are ", n.cats, " items in the category vector."
      ), call. = FALSE)
    }

    if (min(cats) < 2) {
      stop("Minimum number of categories for each item is 2", call. = FALSE)
    }

    # Probabilities for each category at the 1st item
    p <- as.double(prob[1, 1:cats[1]])

    # Create a temporary vector to contain probabilities
    tmp <- c()

    # Possible observed score range for a first item
    obs.range <- 0:(cats[1] - 1)

    # proceed the further analysis when # of the items > 1
    if (n.cats > 1) {
      # Calculate probabilities to earn an observed score by accumulating over remaining items
      for (j in 2:n.prob) { # should start from a 2nd item

        # Probability to earn zero score. This is a special case
        tmp[1] <- p[1] * prob[j, 1]

        # The range of category scores for an added item
        cat.range <- 0:(cats[j] - 1)

        # Possible minimum and maximum observed scores when the item is added
        # but, except zero and perfect score
        obs_rg <- range(obs.range)
        cat_rg <- range(cat.range)
        min.obs <- obs_rg[1]
        max.obs <- obs_rg[2]
        min.cat <- cat_rg[1]
        max.cat <- cat_rg[2]
        min.s <- min.obs + min.cat + 1
        max.s <- max.obs + max.cat - 1

        # possible observed score range except zero and perfect score
        poss.score <- min.s:max.s
        length.score <- length(poss.score)
        # score.mat <- matrix(poss.score, nrow=length(min.cat:max.cat), ncol=length.score, byrow=TRUE)
        score.mat <- col(array(NA, c(cats[j], length.score)))

        # Difference between the possible observed score and each category score if the added item
        poss.diff <- score.mat - min.cat:max.cat

        # The difference above should be greater than or equal to the minimum observed score
        # where the item is added, and less than or equal to the maximum observed score where
        # the item is added
        cols <- poss.diff >= min.obs & poss.diff <= max.obs

        # Final probability to earn the observed score
        prob.score <- c()
        for (k in 1:length.score) {
          tmp.cols <- which(cols[, k])
          tmp.diff <- poss.diff[tmp.cols, k]
          prob.score[k] <- sum(p[(tmp.diff + 1)] * prob[j, tmp.cols])
        }
        tmp[1 + poss.score] <- prob.score

        # Probability to earn perfect score. This is a special case.
        tmp[max.s + 2] <- p[length(p)] * prob[j, cats[j]]

        # Update probabilities
        p <- tmp

        # Update the range of possible observed scores
        # obs.range <- (min.s-1):(max.s+1)
        obs.range <- c((min.s - 1), poss.score, (max.s + 1))

        # Reset the temporary vector
        tmp <- c()
      }
    }

    # return the results
    names(p) <- paste0("score.", obs.range)
    p
  } else {
    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # count N of thetas
    n.theta <- length(theta)

    # prepare the probability matrices
    prepdat <- prep4lw2(x = x, theta = theta, D = D)

    # apply the recursion formula
    p <- lwRecurive(
      prob.cats = prepdat$prob.cat,
      cats = prepdat$cats,
      n.theta = n.theta
    )

    # return the results
    p
  }
}


### --------------------------------------------------------------------------------------------------------------------
# "lwRecursive" function
# Arg:
# prob.cats: (list) each element of a list includes a probability matrix (or data.frame) for an item.
# 					In each probability matrix, rows indicate theta values (or weights) and columns indicate
# 					score categories. Thus, each elements in the matrix represents the probability that a person
# 					who has a certain ability earns a certain score.
# cats: (vector) containes score categories for all items
# n.theta: (integer) number of thetas
#' @importFrom Rfast rowsums
lwRecurive <- function(prob.cats, cats, n.theta) {
  # if(length(unique(sapply(prob.cats, nrow))) != 1L) {
  #   stop("At least, one probability matrix (or data.frame) has the difference number of rows across all marices in a list", call.=FALSE)
  # }

  # if(length(prob.cats) != length(cats)) {
  #   stop(paste0("There are ", nrow(prob.cats), " items in the probability matrix (or data.frame) at each element of a list, ",
  #               "whereas there are ", length(cats), " items in the category vector."), call.=FALSE)
  # }

  if (min(cats) < 2) {
    stop("Minimum number of categories for each item is 2", call. = FALSE)
  }

  # Probabilities for each category at 1st item
  p <- prob.cats[[1]]

  # the number of theta values
  # n.theta <- nrow(prob.cats[[1]])

  # possible observed score range for a first item
  obs.range <- 0:(cats[1] - 1)

  # proceed the further analysis when # of the items > 1
  if (length(cats) > 1) {
    # create a temporary matrix to contain all probabilities
    tScore.range <- 0:sum(cats - 1)
    # tmp <- matrix(0, nrow=n.theta, ncol=length(tScore.range))
    tmp <- array(0, c(n.theta, length(tScore.range)))

    # Calculate probabilities to earn an observed score by accumulating over remaining items
    for (j in 2:length(cats)) {
      # Probability to earn zero score. This is a special case.
      tmp[, 1] <- p[, 1] * prob.cats[[j]][, 1]

      # The range of category scores for an added item
      cat.range <- 0:(cats[j] - 1)

      # Possible minimum and maximum observed scores when the item is added
      # but, except zero and perfect score
      obs_rg <- range(obs.range)
      cat_rg <- range(cat.range)
      min.obs <- obs_rg[1]
      max.obs <- obs_rg[2]
      min.cat <- cat_rg[1]
      max.cat <- cat_rg[2]
      min.s <- min.obs + min.cat + 1
      max.s <- max.obs + max.cat - 1

      # possible observed score range except zero and perfect score
      poss.score <- min.s:max.s
      length.score <- length(poss.score)
      # score.mat <- matrix(poss.score, nrow=length(min.cat:max.cat), ncol=length.score, byrow=TRUE)
      # score.mat <- matrix(poss.score, nrow=cats[j], ncol=length.score, byrow=TRUE)
      score.mat <- col(array(NA, c(cats[j], length.score)))

      # Difference between the possible observed score and each category score if the added item
      # poss.diff <- score.mat - min.cat:max.cat
      poss.diff <- score.mat - cat.range

      # The difference above should be greater than or equal to the minimum observed score
      # where the item is added, and less than or equal to the maximum observed score where
      # the item is added
      cols <- poss.diff >= min.obs & poss.diff <= max.obs

      # Final probability to earn the observed score
      prob.score <- array(0, c(n.theta, length.score))
      for (k in 1:length.score) {
        tmp.cols <- cols[, k]
        tmp.diff <- poss.diff[tmp.cols, k]
        prob.score[, k] <- Rfast::rowsums(p[, (tmp.diff + 1), drop = FALSE] * prob.cats[[j]][, tmp.cols, drop = FALSE])
      }
      tmp[, 1 + poss.score] <- prob.score

      # cols2 <- c(cols)
      # diff2 <- poss.diff[cols2]
      # cols3 <- rep(c(cat.range + 1), length.score)[cols2]
      # prob.score2 <- p[, (diff2 + 1), drop=FALSE] * prob.cats[[j]][, cols3, drop=FALSE]
      # length.col <- Rfast::colsums(cols)
      # end.col <- cumsum(length.col)
      # start.col <- end.col - length.col + 1
      # prob.score2 <-
      #   purrr::map2(.x=start.col, .y=end.col,
      #               .f=function(x, y) Rfast::rowsums(prob.score2[, x:y])) %>%
      # prob.score2 <- do.call(what="cbind", prob.score2)
      # prob.score2 <-
      #   mapply(FUN = function(x, y) {Rfast::rowsums(prob.score2[, x:y])}, x=start.col, y=end.col)

      # Probability to earn perfect score. This is a special case.
      tmp[, max.s + 2] <- p[, ncol(p)] * prob.cats[[j]][, cats[j]]

      # Update the range of possible observed scores
      # obs.range <- (min.s-1):(max.s+1)
      obs.range <- c((min.s - 1), poss.score, (max.s + 1))

      # Update probabilities
      # p <- tmp[, 1:length(obs.range), drop=FALSE]
      p <- tmp[, obs.range + 1, drop = FALSE]
    }
  }

  # return the results
  # colnames(p) <- paste0("score.", 0:sum(cats-1))
  colnames(p) <- paste0("score.", obs.range)
  rownames(p) <- paste0("theta.", 1:n.theta)
  t(p)
}

# "prep4lw2" function
# This function is used only for "lwrc" function
prep4lw2 <- function(x, theta, D) {
  # break down the item metadata into several components
  elm_item <- breakdown(x)

  # compute category probabilities for all items
  prob.cats <- trace(elm_item, theta, D = D, tcc = FALSE)$prob.cats

  # extract score categories for all items
  cats <- elm_item$cats

  # Return results
  list(prob.cats = prob.cats, cats = cats)
}


# "lw_extend" function
# Extend a score distribution by convolving it with one item's probability matrix.
# This is the core polynomial-multiplication step shared by both the forward pass,
# the backward pass, and the final combine step in lwrc_noitem().
#
# p         : n.theta x (Sa+1)     current cumulative score distribution
# prob_item : n.theta x cats_item  category probabilities for the item being added
#             (or, in the combine step, a score distribution treated as probabilities)
# cats_item : integer              number of categories (= ncol of prob_item)
# n.theta   : integer              number of quadrature points
# Returns   : n.theta x (Sa + cats_item)   updated score distribution
#' @importFrom Rfast rowsums
lw_extend <- function(p, prob_item, cats_item, n.theta) {
  Sa     <- ncol(p) - 1L                  # current maximum score
  Sc     <- Sa + cats_item - 1L           # maximum score after adding this item
  result <- matrix(0, nrow = n.theta, ncol = Sc + 1L)

  # score 0: only achievable when both the prefix and the new item score zero
  result[, 1L] <- p[, 1L] * prob_item[, 1L]

  # score Sc (perfect): only achievable when both parts are at their maximum
  result[, Sc + 1L] <- p[, Sa + 1L] * prob_item[, cats_item]

  # middle scores s = 1 ... Sc-1: sum over all valid (prefix_score, item_score) pairs
  if (Sc > 1L) {
    for (s in seq_len(Sc - 1L)) {
      k_min  <- max(0L, s - Sa)           # smallest valid item category score
      k_max  <- min(cats_item - 1L, s)    # largest valid item category score
      k_seq  <- k_min:k_max
      p_idx  <- s - k_seq + 1L           # column indices into p  (prefix score = s - k)
      pb_idx <- k_seq + 1L               # column indices into prob_item (item score = k)
      # Rfast::rowsums vectorises the weighted sum across all n.theta simultaneously
      result[, s + 1L] <- Rfast::rowsums(
        p[, p_idx, drop = FALSE] * prob_item[, pb_idx, drop = FALSE]
      )
    }
  }

  result
}


# "lwrc_noitem" function
# Compute lkhd_noitem for all J items using a single forward-backward pass instead
# of J separate lwRecurive() calls.  Reduces the dominant cost in sx2_fit() from
# O(J^3 K^2 Q) to approximately O(J^2 K^2 Q / 6) while preserving exact results.
#
# Algorithm:
#   Forward pass  : fwd[[i]] = score distribution for items 1 ... (i-1)
#   Backward pass : bwd[[i]] = score distribution for items i ... J
#   Combine       : lkhd_noitem[[i]] = convolve(fwd[[i]], bwd[[i+1]])
#
# prob.cats : length-J list, each element n.theta x cats[j]
# cats      : integer vector of category counts (length J)
# n.theta   : integer, number of quadrature points (nquad)
# Returns   : length-J list; element i is n.theta x (t.score - cats[i] + 2)
#             (same format as t(lwRecurive(prob.cats[-i], cats[-i], n.theta)))
#' @importFrom Rfast rowsums
lwrc_noitem <- function(prob.cats, cats, n.theta) {
  J          <- length(cats)
  identity_m <- matrix(1, nrow = n.theta, ncol = 1L)  # empty-set distribution: P(score=0)=1

  # FORWARD PASS -----------------------------------------------------------
  # fwd_list[[i]] holds the score distribution for items 1 ... (i-1);
  # fwd_list[[1]] is the identity (empty prefix, score 0 with probability 1)
  fwd_list <- vector("list", J + 1L)
  fwd_list[[1L]] <- identity_m
  for (i in seq_len(J)) {
    fwd_list[[i + 1L]] <- lw_extend(fwd_list[[i]], prob.cats[[i]], cats[i], n.theta)
  }

  # BACKWARD PASS ----------------------------------------------------------
  # bwd_list[[i]] holds the score distribution for items i ... J;
  # bwd_list[[J+1]] is the identity (empty suffix)
  # Because score distributions are additive (convolution is commutative),
  # processing items in reverse order yields the same marginal distribution.
  bwd_list <- vector("list", J + 1L)
  bwd_list[[J + 1L]] <- identity_m
  for (i in seq(J, 1L, by = -1L)) {
    bwd_list[[i]] <- lw_extend(bwd_list[[i + 1L]], prob.cats[[i]], cats[i], n.theta)
  }

  # COMBINE ----------------------------------------------------------------
  # lkhd_noitem[[i]] = convolution of the prefix (items before i) and
  #                    the suffix (items after i).
  # The right distribution is passed to lw_extend() as if it were a
  # probability matrix; the polynomial-multiplication formula is identical.
  lkhd_noitem <- vector("list", J)
  for (i in seq_len(J)) {
    left  <- fwd_list[[i]]        # items 1 ... (i-1)
    right <- bwd_list[[i + 1L]]   # items (i+1) ... J

    if (ncol(left) == 1L) {
      # i = 1: no items to the left; result is just the suffix distribution
      lkhd_noitem[[i]] <- right
    } else if (ncol(right) == 1L) {
      # i = J: no items to the right; result is just the prefix distribution
      lkhd_noitem[[i]] <- left
    } else {
      # general case: polynomial convolution of two score distributions
      lkhd_noitem[[i]] <- lw_extend(left, right, ncol(right), n.theta)
    }
  }

  lkhd_noitem
}
