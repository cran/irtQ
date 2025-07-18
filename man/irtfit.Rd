% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irtfit.R
\name{irtfit}
\alias{irtfit}
\alias{irtfit.default}
\alias{irtfit.est_item}
\alias{irtfit.est_irt}
\title{Traditional IRT Item Fit Statistics}
\usage{
irtfit(x, ...)

\method{irtfit}{default}(
  x,
  score,
  data,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  D = 1,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  pcm.loc = NULL,
  ...
)

\method{irtfit}{est_item}(
  x,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  pcm.loc = NULL,
  ...
)

\method{irtfit}{est_irt}(
  x,
  score,
  group.method = c("equal.width", "equal.freq"),
  n.width = 10,
  loc.theta = "average",
  range.score = NULL,
  alpha = 0.05,
  missing = NA,
  overSR = 2,
  min.collapse = 1,
  pcm.loc = NULL,
  ...
)
}
\arguments{
\item{x}{A data frame containing item metadata (e.g., item parameters,
number of categories, IRT model types, etc.); or an object of class
\code{est_irt} obtained from \code{\link[=est_irt]{est_irt()}}, or \code{est_item} from
\code{\link[=est_item]{est_item()}}.

See \code{\link[=est_irt]{est_irt()}} or \code{\link[=simdat]{simdat()}} for more details about the item
metadata. This data frame can be easily created using the
\code{\link[=shape_df]{shape_df()}} function.}

\item{...}{Further arguments passed to or from other methods.}

\item{score}{A numeric vector containing examinees' ability estimates
(theta values).}

\item{data}{A matrix of examinees' item responses corresponding to the items
specified in the \code{x} argument. Rows represent examinees and columns
represent items.}

\item{group.method}{A character string specifying the method used to group
examinees along the ability scale when computing the \eqn{\chi^{2}} and
\eqn{G^{2}} fit statistics. Available options are:
\itemize{
\item \code{"equal.width"}: Divides the ability scale into intervals of equal width.
\item \code{"equal.freq"}: Divides the examinees into groups with (approximately)
equal numbers of examinees.
}

Note that \code{"equal.freq"} does not guarantee exactly equal group sizes
due to ties in ability estimates. Default is \code{"equal.width"}.
The number of groups and the range of the ability scale are controlled by
the \code{n.width} and \code{range.score} arguments, respectively.}

\item{n.width}{An integer specifying the number of intervals (groups) into which
the ability scale is divided for computing the fit statistics. Default is 10.}

\item{loc.theta}{A character string indicating the point on the ability scale
at which the expected category probabilities are calculated for each group.
Available options are:
\itemize{
\item \code{"average"}: Uses the average ability estimate of examinees within each group.
\item \code{"middle"}: Uses the midpoint of each group's ability interval.
}

Default is \code{"average"}.}

\item{range.score}{A numeric vector of length two specifying the lower and upper
bounds of the ability scale. Ability estimates below the lower bound or above the
upper bound are truncated to the respective bound. If \code{NULL}, the observed minimum
and maximum of the \code{score} vector are used. Note that this range restriction is
independent of the grouping method specified in \code{group.method}.
Default is \code{NULL}.}

\item{D}{A scaling constant used in IRT models to make the logistic function
closely approximate the normal ogive function. A value of 1.7 is commonly
used for this purpose. Default is 1.}

\item{alpha}{A numeric value specifying the significance level (\eqn{\alpha}) for
the hypothesis tests of the \eqn{\chi^{2}} and \eqn{G^{2}} item fit statistics.
Default is \code{0.05}.}

\item{missing}{A value indicating missing responses in the data set. Default
is \code{NA}. See \strong{Details} below.}

\item{overSR}{A numeric threshold used to identify ability groups (intervals)
whose standardized residuals exceed the specified value. This is used to
compute the proportion of misfitting groups per item. Default is 2.}

\item{min.collapse}{An integer specifying the minimum expected frequency required
for a cell in the contingency table. Neighboring groups will be merged if any
expected cell frequency falls below this threshold when computing the
\eqn{\chi^{2}} and \eqn{G^{2}} statistics. Default is 1.}

\item{pcm.loc}{A vector of integers indicating the locations (indices) of
partial credit model (PCM) items for which slope parameters are fixed.}
}
\value{
This function returns an object of class \code{irtfit}, which includes
the following components:

\item{fit_stat}{A data frame containing the results of three IRT item fit statistics—
\eqn{\chi^{2}}, \eqn{G^{2}}, infit, and outfit—for all evaluated items.
Each row corresponds to one item, and the columns include:
the item ID; \eqn{\chi^{2}} statistic; \eqn{G^{2}} statistic;
degrees of freedom for \eqn{\chi^{2}} and \eqn{G^{2}};
critical values and p-values for both statistics;
outfit and infit values;
the number of examinees used to compute these statistics;
and the proportion of ability groups (prior to cell collapsing) that have
standardized residuals greater than the threshold specified in the \code{overSR} argument.}

\item{contingency.fitstat}{A list of contingency tables used to compute the
\eqn{\chi^{2}} and \eqn{G^{2}} fit statistics for all items. Note that the
cell-collapsing strategy is applied to these tables to ensure sufficient
expected frequencies.}

\item{contingency.plot}{A list of contingency tables used to generate raw
and standardized residual plots (Hambleton et al., 1991) via the
\code{\link[=plot.irtfit]{plot.irtfit()}}. Note that these tables are based on the original,
uncollapsed groupings.}

\item{individual.info}{A list of data frames containing individual residuals
and corresponding variance values. This information is used to compute infit
and outfit statistics.}

\item{item_df}{A data frame containing the item metadata provided in the
argument \code{x}.}

\item{ancillary}{A list of ancillary information used during the item fit
analysis.}
}
\description{
This function computes traditional IRT item fit statistics, including the
\eqn{\chi^{2}} fit statistic (e.g., Bock, 1960; Yen, 1981),
the log-likelihood ratio \eqn{\chi^{2}} fit statistic (\eqn{G^{2}}; McKinley
& Mills, 1985), and the infit and outfit statistics (Ames et al., 2015). It
also returns contingency tables used to compute the \eqn{\chi^{2}} and
\eqn{G^{2}} statistics.
}
\details{
To compute the \eqn{\chi^2} and \eqn{G^2} item fit statistics, the \code{group.method}
argument determines how the ability scale is divided into groups:
\itemize{
\item \code{"equal.width"}: Examinees are grouped based on intervals of equal width a
long the ability scale.
\item \code{"equal.freq"}: Examinees are grouped such that each group contains
(approximately) the same number of individuals.
}

Note that \code{"equal.freq"} does not guarantee \emph{exactly} equal frequencies across
all groups, since grouping is based on quantiles.

When dividing the ability scale into intervals to compute the \eqn{\chi^2}
and \eqn{G^2} fit statistics, the intervals should be:
\itemize{
\item \strong{Wide enough} to ensure that each group contains a sufficient number of
examinees (to avoid unstable estimates),
\item \strong{Narrow enough} to ensure that examinees within each group are relatively
homogeneous in ability (Hambleton et al., 1991).
}

If you want to divide the ability scale into a number of groups other than
the default of 10, specify the desired number using the \code{n.width} argument.
For reference:
\itemize{
\item Yen (1981) used 10 fixed-width groups,
\item Bock (1960) allowed for flexibility in the number of groups.
}

Regarding degrees of freedom (\emph{df}):
\itemize{
\item The \eqn{\chi^2} statistic is approximately chi-square distributed with
degrees of freedom equal to the number of ability groups minus the number
of item parameters (Ames et al., 2015).
\item The \eqn{G^2} statistic is approximately chi-square distributed with
degrees of freedom equal to the number of ability groups (Ames et al., 2015;
Muraki & Bock, 2003).

Note that if \code{"DRM"} is specified for an item in the item metadata set,
the item is treated as a \code{"3PLM"} when computing the degrees of freedom
for the \eqn{\chi^{2}} fit statistic.
}

Note that infit and outfit statistics should be interpreted with caution when
applied to non-Rasch models. The returned object—particularly the contingency
tables—can be passed to \code{\link[=plot.irtfit]{plot.irtfit()}} to generate raw and standardized
residual plots (Hambleton et al., 1991).
}
\section{Methods (by class)}{
\itemize{
\item \code{irtfit(default)}: Default method for computing traditional IRT item fit
statistics using a data frame \code{x} that contains item metadata.

\item \code{irtfit(est_item)}: An object created by the function \code{\link[=est_item]{est_item()}}.

\item \code{irtfit(est_irt)}: An object created by the function \code{\link[=est_irt]{est_irt()}}.

}}
\examples{
\donttest{
## Example 1
## Use the simulated CAT data
# Identify items with more than 10,000 responses
over10000 <- which(colSums(simCAT_MX$res.dat, na.rm = TRUE) > 10000)

# Select items with more than 10,000 responses
x <- simCAT_MX$item.prm[over10000, ]

# Extract response data for the selected items
data <- simCAT_MX$res.dat[, over10000]

# Extract examinees' ability estimates
score <- simCAT_MX$score

# Compute item fit statistics
fit1 <- irtfit(
  x = x, score = score, data = data, group.method = "equal.width",
  n.width = 10, loc.theta = "average", range.score = NULL, D = 1, alpha = 0.05,
  missing = NA, overSR = 2
)

# View the fit statistics
fit1$fit_stat

# View the contingency tables used to compute fit statistics
fit1$contingency.fitstat


## Example 2
## Import the "-prm.txt" output file from flexMIRT
flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")

# Select the first two dichotomous items and the last polytomous item
x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df[c(1:2, 55), ]

# Generate ability values from a standard normal distribution
set.seed(10)
score <- rnorm(1000, mean = 0, sd = 1)

# Simulate response data
data <- simdat(x = x, theta = score, D = 1)

# Compute item fit statistics
fit2 <- irtfit(
  x = x, score = score, data = data, group.method = "equal.freq",
  n.width = 11, loc.theta = "average", range.score = c(-4, 4), D = 1, alpha = 0.05
)

# View the fit statistics
fit2$fit_stat

# View the contingency tables used to compute fit statistics
fit2$contingency.fitstat

# Plot raw and standardized residuals for the first item (dichotomous)
plot(x = fit2, item.loc = 1, type = "both", ci.method = "wald",
     show.table = TRUE, ylim.sr.adjust = TRUE)

# Plot raw and standardized residuals for the third item (polytomous)
plot(x = fit2, item.loc = 3, type = "both", ci.method = "wald",
     show.table = FALSE, ylim.sr.adjust = TRUE)
}

}
\references{
Ames, A. J., & Penfield, R. D. (2015). An NCME Instructional
Module on Item-Fit Statistics for Item Response Theory Models.
\emph{Educational Measurement: Issues and Practice, 34}(3), 39-48.

Bock, R.D. (1960), \emph{Methods and applications of optimal scaling}. Chapel
Hill, NC: L.L. Thurstone Psychometric Laboratory.

Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of
item response theory}. Newbury Park, CA: Sage.

McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit
statistics.
\emph{Applied Psychological Measurement, 9}, 49-57.

Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test
scoring for rating scale data (Computer Software). Chicago, IL: Scientific
Software International. URL http://www.ssicentral.com

Wells, C. S., & Bolt, D. M. (2008). Investigation of a nonparametric
procedure for assessing goodness-of-fit in item response theory. \emph{Applied
Measurement in Education, 21}(1), 22-40.

Yen, W. M. (1981). Using simulation results to choose a latent trait model.
\emph{Applied Psychological Measurement, 5}, 245-262.
}
\seealso{
\code{\link[=plot.irtfit]{plot.irtfit()}}, \code{\link[=shape_df]{shape_df()}}, \code{\link[=est_irt]{est_irt()}},
\code{\link[=est_item]{est_item()}}
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
