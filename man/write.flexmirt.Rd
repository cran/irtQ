% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_flexmirt.R
\name{write.flexmirt}
\alias{write.flexmirt}
\title{Write a "-prm.txt" file for flexMIRT}
\usage{
write.flexmirt(
  x,
  file = NULL,
  norm.pop = c(0, 1),
  rePar = TRUE,
  mgroup = FALSE,
  group.name = NULL
)
}
\arguments{
\item{x}{A data frame containing the item metadata (e.g., item parameters, number of categories, models ...) for a single group or
a list of the item metadata for multiple groups. See \code{\link{est_irt}}, \code{\link{irtfit}}, \code{\link{info}},
or \code{\link{simdat}} for more details about the item metadata. The item metadata can be easily created using the function
\code{\link{shape_df}}.}

\item{file}{The destination file name.}

\item{norm.pop}{A numeric vector of two components specifying a mean and standard deviation of the normal
population ability distribution for a single group or a list of the numeric vectors of length two for multiple groups.
When a list is provided, each internal numeric vector should contain a mean and standard deviation of the ability
distribution of each group (e.g., \code{norm.pop = list(c(0, 1), c(0, 0.8), c(0.5, 1.2)}) for three groups). When \code{mgroup = TRUE}
and a single vector of length two is provided (e.g., \code{norm.pop = c(0, 1)}), the same vector will be recycled across all groups.
Default is c(0,1).}

\item{rePar}{A logical value indicating whether the item parameters in the item metadata
are the reparameterized item parameters. If TRUE, the item intercepts and logits of item guessing parameters
should be included in the item metadata. If FALSE, the item difficulty and item guessing parameters
should be included in the item metadata.}

\item{mgroup}{A logical value indicating whether a "-prm.txt" file is created for a single group or multiple groups.
Default is FALSE.}

\item{group.name}{A character vector of group names. If NULL, the group names are automatically generated (e.g., Group1).}
}
\value{
A "-prm.txt" file.
}
\description{
This function writes an output file of "-prm.txt" for flexMIRT (Cai, 2017). The current version of this function
can be used only for the unidimensional IRT models. This function was written by modifying the function \code{read.flexmirt}
(Pritikin & Falk, 2020).
}
\examples{
\donttest{
## 1. Create "-prm.txt" file for a single group
##    using the simulated CAT data
# 1-(1) extract the item metadata
x <- simCAT_MX$item.prm

# 1-(2) set a name of "-prm.txt" file
temp_prm <- file.path(tempdir(), "single_group_temp-prm.txt")

# 1-(3) write out the "-prm.txt" file
write.flexmirt(x, file = temp_prm, norm.pop = c(0, 1), rePar = FALSE)

## 2. Create "-prm.txt" file for multiple groups
##    using the simulated three multiple group data
# 2-(1) extract the item metadata
x <- simMG$item.prm

# set a name of "-prm.txt" file
temp_prm <- file.path(tempdir(), "mg_group_temp-prm1.txt")

# write out the "-prm.txt" file
write.flexmirt(x,
  file = temp_prm, norm.pop = list(c(0, 1), c(0.5, 0.8), c(-0.3, 1.3)),
  rePar = FALSE, mgroup = TRUE, group.name = c("GR1", "GR2", "GR3")
)

# or write out the "-prm.txt" file so that
# all groups have the same ability distributions
# and the group names are generate autoumatically
temp_prm <- file.path(tempdir(), "mg_group_temp-prm2.txt")
write.flexmirt(x,
  file = temp_prm, norm.pop = c(0, 1),
  rePar = FALSE, mgroup = TRUE, group.name = NULL
)
}

}
\references{
Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
Chapel Hill, NC: Vector Psychometric Group.

Pritikin, J. N., & Falk, C. F. (2020). OpenMx: A modular research environment for item response theory
method development. \emph{Applied Psychological Measurement, 44}(7-8), 561-562.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
