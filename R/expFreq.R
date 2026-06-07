# This function returns a contingency table of the expected frequencies
# to be used to compute S-X2 fit statistic
expFreq <- function(t.score, cats, prob.cats, lkhd_noitem, lkhd, wts, score.freq) {
  # number of possible scores for the test without this item
  t_noitem <- t.score - cats + 2L

  # compute joint[k, j] = sum_theta w_theta * P(cat=j-1|theta) * P(rest_score=k-1|theta)
  # one BLAS DGEMM replaces K separate Rfast::colsums() calls
  joint <- crossprod(lkhd_noitem, prob.cats * wts[, 2])  # (t_noitem x cats)

  # fill the (t.score+1) x cats staircase matrix in one vectorised assignment:
  # category j occupies rows j:(j + t_noitem - 1) of column j
  row_idx <- rep(seq_len(t_noitem), cats) + rep(seq_len(cats) - 1L, each = t_noitem)
  col_idx <- rep(seq_len(cats), each = t_noitem)
  tmp1 <- array(0, c(t.score + 1, cats))
  tmp1[cbind(row_idx, col_idx)] <- joint

  # divide each row by the marginal score distribution (one BLAS DGEMV);
  # R recycles denom column-wise so each row i is divided by denom[i]
  denom <- as.vector(crossprod(lkhd, wts[, 2]))  # length t.score+1
  tmp2  <- tmp1 / denom

  colnames(tmp2) <- paste0("score.", 0:(cats - 1))
  rownames(tmp2) <- paste0("score.", 0:t.score)

  tmp2 <- score.freq * tmp2
  tmp2 <- tmp2[c(-1, -nrow(tmp2)), ]

  row.first <- purrr::map_dbl(1:cats, .f = function(i) sum(tmp2[1:(cats - 1), i]))
  row.end <- purrr::map_dbl(1:cats, .f = function(i) sum(tmp2[nrow(tmp2):(nrow(tmp2) - cats + 2), i]))

  first.name <- rownames(tmp2)[cats - 1]
  last.name <- rownames(tmp2)[nrow(tmp2) - cats + 2]

  tmp2 <- tmp2[-c(1:(cats - 1), nrow(tmp2):(nrow(tmp2) - cats + 2)), ]
  tmp3 <- rbind(row.first, tmp2, row.end)
  rownames(tmp3) <- c(first.name, rownames(tmp2), last.name)

  data.frame(tmp3)
}

# This function returns a contingency table of the observed frequencies
# to be used to compute S-X2 fit statistic
obsFreq <- function(rawscore, response, t.score, cats) {
  # Encode (rawscore, response) as a row-major linear index into the
  # (t.score+1) x cats matrix, then rebuild with a single tabulate() call;
  # avoids factor allocation and table() overhead called J times
  lin_idx <- rawscore * cats + response + 1L
  tmp2 <- matrix(tabulate(lin_idx, nbins = (t.score + 1L) * cats),
                 nrow = t.score + 1L, ncol = cats, byrow = TRUE)

  colnames(tmp2) <- paste0("score.", 0:(cats - 1))
  rownames(tmp2) <- paste0("score.", 0:t.score)

  tmp2 <- tmp2[c(-1, -nrow(tmp2)), ]

  row.first <- purrr::map_dbl(1:cats, .f = function(i) sum(tmp2[1:(cats - 1), i]))
  row.end <- purrr::map_dbl(1:cats, .f = function(i) sum(tmp2[nrow(tmp2):(nrow(tmp2) - cats + 2), i]))

  first.name <- rownames(tmp2)[cats - 1]
  last.name <- rownames(tmp2)[nrow(tmp2) - cats + 2]

  tmp2 <- tmp2[-c(1:(cats - 1), nrow(tmp2):(nrow(tmp2) - cats + 2)), ]
  tmp3 <- rbind(row.first, tmp2, row.end)
  rownames(tmp3) <- c(first.name, rownames(tmp2), last.name)

  data.frame(tmp3)
}
