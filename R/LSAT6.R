#' LSAT6 Data
#'
#' A well-known dichotomous response dataset from the Law School Admission Test (LSAT),
#' Section 6, as used in Thissen (1982).
#'
#' @usage LSAT6
#'
#' @format A data frame with 1,000 rows and 5 columns, where each row represents a unique
#' examinee's response pattern to five dichotomously scored items (0 = incorrect, 1 = correct).
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
#' *Psychometrika, 47*, 175-186.
#'
#' @examples
#' # structure of the data
#' head(LSAT6)
#'
#' # fit the 2PL model to the LSAT6 data
#' est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2, verbose = FALSE)
#'
"LSAT6"
