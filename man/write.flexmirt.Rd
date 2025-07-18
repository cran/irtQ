% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_flexmirt.R
\name{write.flexmirt}
\alias{write.flexmirt}
\title{Write a "-prm.txt" File for flexMIRT}
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
\item{x}{A data frame of item metadata (e.g., item parameters, number of
categories, model types) for a single group, or a list of such data frames
for multiple groups. See \code{\link[=est_irt]{est_irt()}} or \code{\link[=simdat]{simdat()}} for item
metadata format. You can also create metadata using \code{\link[=shape_df]{shape_df()}}.}

\item{file}{A character string specifying the destination file path (with a
".txt" extension).}

\item{norm.pop}{A numeric vector of length two specifying the mean and
standard deviation of the normal population ability distribution for a
single group, or a list of such vectors for multiple groups. When a list is
provided, each internal vector must contain the mean and standard deviation
for a group's ability distribution (e.g., \code{norm.pop = list(c(0, 1), c(0, 0.8), c(0.5, 1.2))} for three groups). If \code{mgroup = TRUE} and a single
vector is provided (e.g., \code{norm.pop = c(0, 1)}), it will be recycled across
all groups. The default is \code{c(0, 1)}.}

\item{rePar}{A logical value indicating whether the item parameters are
reparameterized. If \code{TRUE}, item intercepts and logits of guessing
parameters are assumed. If \code{FALSE}, item difficulty and guessing parameters
are assumed.}

\item{mgroup}{A logical value indicating whether the file includes multiple
groups. Default is \code{FALSE}.}

\item{group.name}{A character vector of group names. If \code{NULL}, group names
are automatically generated (e.g., "Group1", "Group2", ...).}
}
\value{
This function creates a flexMIRT-style "-prm.txt" file at the
specified path.
}
\description{
This function writes a flexMIRT-compatible "-prm.txt" file (Cai,
2017). It currently supports only unidimensional IRT models. This function
was developed by modifying \code{read.flexmirt()} from Pritikin & Falk (2020).
}
\examples{
\donttest{
## 1. Create a "-prm.txt" file for a single group
##    using the simulated CAT data
# 1-(1) Extract the item metadata
x <- simCAT_MX$item.prm

# 1-(2) Set the name of the "-prm.txt" file
temp_prm <- file.path(tempdir(), "single_group_temp-prm.txt")

# 1-(3) Write the "-prm.txt" file
write.flexmirt(x, file = temp_prm, norm.pop = c(0, 1), rePar = FALSE)

## 2. Create a "-prm.txt" file for multiple groups
##    using simulated multi-group data
# 2-(1) Extract the item metadata
x <- simMG$item.prm

# Set the name of the "-prm.txt" file
temp_prm <- file.path(tempdir(), "mg_group_temp-prm1.txt")

# Write the "-prm.txt" file
write.flexmirt(x,
  file = temp_prm, norm.pop = list(c(0, 1), c(0.5, 0.8), c(-0.3, 1.3)),
  rePar = FALSE, mgroup = TRUE, group.name = c("GR1", "GR2", "GR3")
)

# Or write the "-prm.txt" file so that
# all groups share the same ability distribution
# and group names are generated automatically
temp_prm <- file.path(tempdir(), "mg_group_temp-prm2.txt")
write.flexmirt(x,
  file = temp_prm, norm.pop = c(0, 1),
  rePar = FALSE, mgroup = TRUE, group.name = NULL
)
}

}
\references{
Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional
item analysis and test scoring (Computer Software). Chapel Hill, NC: Vector
Psychometric Group.

Pritikin, J. N., & Falk, C. F. (2020). OpenMx: A modular research
environment for item response theory method development. \emph{Applied
Psychological Measurement, 44}(7-8), 561-562.
}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
