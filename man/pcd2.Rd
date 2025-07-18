% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcd2.R
\name{pcd2}
\alias{pcd2}
\title{Pseudo-count D2 method}
\usage{
pcd2(
  x,
  data,
  D = 1,
  item.skip = NULL,
  missing = NA,
  Quadrature = c(49, 6),
  weights = NULL,
  group.mean = 0,
  group.var = 1,
  crit.val = NULL,
  min.resp = NULL,
  purify = FALSE,
  max.iter = 10,
  verbose = TRUE
)
}
\arguments{
\item{x}{A data frame containing item metadata (e.g., item parameters, number
of categories, IRT model types, etc.). See \code{\link[=est_irt]{est_irt()}} or
\code{\link[=simdat]{simdat()}} for more details about the item metadata. This data frame
can be easily created using the \code{\link[=shape_df]{shape_df()}} function.}

\item{data}{A matrix of examinees' item responses corresponding to the items
specified in the \code{x} argument. Rows represent examinees and columns
represent items.}

\item{D}{A scaling constant used in IRT models to make the logistic function
closely approximate the normal ogive function. A value of 1.7 is commonly
used for this purpose. Default is 1.}

\item{item.skip}{A numeric vector of item indices to exclude from IPD analysis.
If \code{NULL}, all items are included. Useful for omitting specific items based on
prior insights.}

\item{missing}{A value indicating missing responses in the data set. Default
is \code{NA}.}

\item{Quadrature}{A numeric vector of length two:
\itemize{
\item first element: number of quadrature points
\item second element: symmetric bound (absolute value) for those points
For example, \code{c(49, 6)} specifies 49 evenly spaced points from –6 to 6.
These points are used in the E-step of the EM algorithm. Default is \code{c(49,   6)}.
}}

\item{weights}{A two-column matrix or data frame containing the quadrature
points (in the first column) and their corresponding weights (in the second
column) for the latent variable prior distribution. If not \code{NULL}, the
scale of the latent ability distribution is fixed to match the scale of the
provided quadrature points and weights. The weights and points can be
conveniently generated using the function \code{\link[=gen.weight]{gen.weight()}}.

If \code{NULL}, a normal prior density is used instead, based on the
information provided in the \code{Quadrature}, \code{group.mean}, and \code{group.var}
arguments. Default is \code{NULL}.}

\item{group.mean}{A numeric value specifying the mean of the latent variable
prior distribution when \code{weights = NULL}. Default is 0. This value is fixed
to resolve the indeterminacy of the item parameter scale during
calibration.}

\item{group.var}{A positive numeric value specifying the variance of the
latent variable prior distribution when \code{weights = NULL}. Default is 1.
This value is fixed to resolve the indeterminacy of the item parameter
scale during calibration.}

\item{crit.val}{A critical value applied in hypothesis testing using
the Pseudo-count \eqn{D^{2}} statistic. Default is \code{NULL}.}

\item{min.resp}{A positive integer specifying the minimum required number of
responses for each evaluated item. Defaults to \code{NULL}.}

\item{purify}{Logical. Indicates whether to apply a purification procedure.
Default is \code{FALSE}.}

\item{max.iter}{A positive integer specifying the maximum number of
iterations allowed for the purification process. Default is \code{10}.}

\item{verbose}{Logical. If \code{TRUE}, progress messages from the purification
procedure will be displayed; if \code{FALSE}, the messages will be suppressed.
Default is \code{TRUE}.}
}
\value{
This function returns a list containing four main components:

\item{no_purify}{A list containing the results of Pseudo-count \eqn{D^{2}} analysis
without applying the purification procedure. It includes:
\describe{
\item{ipd_stat}{A data frame summarizing the Pseudo-count \eqn{D^{2}} statistics
for all items. The columns include:
\code{id} (item ID),
\code{pcd2} (the computed \eqn{D^{2}} value), and
\code{N} (the number of valid examinee responses per item).}
\item{ipd_item}{A numeric vector of item indices that were flagged as exhibiting
item parameter drift (IPD), based on the specified critical value \code{crit.val}.
If no items are flagged or \code{crit.val = NULL}, this is \code{NULL}.}
}
}

\item{purify}{A logical value indicating whether the iterative purification
procedure was applied (\code{TRUE}) or not (\code{FALSE}).}

\item{with_purify}{A list containing the results of Pseudo-count \eqn{D^{2}} analysis
after applying the purification procedure. This list is populated only when both
\code{purify = TRUE} and \code{crit.val} is not \code{NULL}. It includes:
\describe{
\item{ipd_stat}{A data frame reporting the final Pseudo-count \eqn{D^{2}} statistics
after purification. Columns include:
\code{id} (item ID),
\code{pcd2} (the computed \eqn{D^{2}} value),
\code{N} (the number of valid responses), and
\code{n.iter} (the iteration number in which each item was evaluated).}
\item{ipd_item}{A numeric vector of item indices flagged as IPD items during
purification. Items are ordered by the iteration in which they were flagged.}
\item{n.iter}{An integer indicating the number of purification iterations
completed.}
\item{complete}{A logical value indicating whether the purification procedure
converged before reaching the maximum number of iterations (\code{max.iter}).
If \code{FALSE}, the iteration limit was reached before convergence.}
}
}

\item{crit.val}{A numeric value indicating the critical threshold used to flag
items for parameter drift. If not specified by the user, this will be \code{NULL}.}
}
\description{
This function calculates the Pseudo-count \eqn{D^{2}} statistic
to evaluate item parameter drift, as described by Cappaert et al. (2018) and
Stone (2000). The Pseudo-count \eqn{D^{2}} statistic is designed to detect
item parameter drift efficiently without requiring item recalibration, making
it especially valuable in computerized adaptive testing (CAT) environments.
This method compares observed and expected response frequencies across
quadrature points, which represent latent ability levels. The expected
frequencies are computed using the posterior distribution of each examinee's
ability (Stone, 2000), providing a robust and sensitive measure of item
parameter drift, ensuring the stability and accuracy of the test over time.
}
\details{
The Pseudo-count \eqn{D^{2}} statistic quantifies item parameter drift (IPD) by
computing the weighted squared differences between the observed and expected
response frequencies for each score category across ability levels. The expected
frequencies are determined using the posterior distribution of each examinee's
ability (Stone, 2000).

The Pseudo-count \eqn{D^{2}} statistic is calculated as:
\deqn{
Pseudo-count D^{2} = \sum_{k=1}^{Q} \left( \frac{r_{0k} + r_{1k}}{N}\right)
\left( \frac{r_{1k}}{r_{0k} + r_{1k}} - E_{1k} \right)^2
}

where \eqn{r_{0k}} and \eqn{r_{1k}} are the pseudo-counts for the incorrect
and correct responses at each ability level \eqn{k}, \eqn{E_{1k}} is the
expected proportion of correct responses at each ability level \eqn{k},
calculated using item parameters from the item bank, and \eqn{N} is the total
count of examinees who received each item.

\strong{Critical Value (\code{crit.val})}:
The \code{crit.val} argument specifies the threshold used to flag an item as
exhibiting potential parameter drift. If an item's Pseudo-count \eqn{D^{2}}
value exceeds this threshold, it is identified as a drifted item. If
\code{crit.val = NULL}, the function reports the raw statistic without flagging.

\strong{Minimum Response Count (\code{min.resp})}:
The \code{min.resp} argument sets a lower bound on the number of responses required
for an item to be included in the analysis. Items with fewer responses than
\code{min.resp} are automatically excluded by replacing all their responses
with \code{NA}. This avoids unreliable estimates based on small sample sizes.

\strong{Purification Procedure}:
Although Cappaert et al. (2018) did not incorporate purification into their method,
\code{\link[=pcd2]{pcd2()}} implements an optional iterative purification process similar
to Lim et al. (2022). When \code{purify = TRUE} and a \code{crit.val} is provided:
\itemize{
\item The procedure begins by identifying items flagged for drift using the initial
Pseudo-count \eqn{D^{2}} statistics.
\item In each subsequent iteration, the item with the highest flagged Pseudo-count
\eqn{D^{2}} value is removed from the item set, and the statistics are recalculated
using only the remaining items.
\item The process continues until no additional items are flagged or the
number of iterations reaches \code{max.iter}.
\item All flagged items and statistics are saved, and convergence status is reported.
}
This process ensures that drift detection is not distorted by already-flagged items,
improving the robustness of the results.
}
\examples{
## Example 1: No critical value specified
## Compute the Pseudo-count D² statistics for dichotomous items
## Import the "-prm.txt" output file generated by flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Extract metadata for the first 30 3PLM items
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[1:30, 1:6]

# Generate abilities for 500 examinees from N(0, 1)
set.seed(25)
score <- rnorm(500, mean = 0, sd = 1)

# Simulate response data using the item metadata and ability values
data <- simdat(x = x, theta = score, D = 1)

# Compute the Pseudo-count D² statistics (no purification applied)
ps_d2 <- pcd2(x = x, data = data)
print(ps_d2)

## Example 2: Applying a critical value with purification
# Compute the Pseudo-count D² statistics with purification enabled
ps_d2_puri <- pcd2(x = x, data = data, crit.val = 0.002, purify = TRUE)
print(ps_d2_puri)

}
\references{
Cappaert, K. J., Wen, Y., & Chang, Y. F. (2018). Evaluating CAT-adjusted
approaches for suspected item parameter drift detection. \emph{Measurement:
Interdisciplinary Research and Perspectives, 16}(4), 226-238.

Stone, C. A. (2000). Monte Carlo based null distribution for an alternative
goodness-of-fit test statistic in IRT models. \emph{Journal of educational
measurement, 37}(1), 58-75.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
