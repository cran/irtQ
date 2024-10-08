% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{bind.fill}
\alias{bind.fill}
\title{Bind Fill}
\usage{
bind.fill(List, type = c("rbind", "cbind"), fill = NA)
}
\arguments{
\item{List}{A list containing different length of numeric vectors}

\item{type}{A character string specifying whether rbind is used or cbind is used.}

\item{fill}{The value used to fill in missing data when aligning datasets. For \code{type = "cbind"},
this fills missing rows in shorter columns. For \code{type = "rbind"}, this fills missing columns in shorter rows.
Accepts any R object (e.g., numeric, character, logical). Defaults to NA.}
}
\value{
A matrix.
}
\description{
This function creates a cbind matrix or rbind matrix using a list containing different length
of numeric vectors.
}
\examples{
# sample list
score_list <- list(item1=c(0:3), item2=c(0:2), item3=c(0:5), item3=c(0:4))

# examples
# 1) create a rbind with the sample score list
bind.fill(score_list, type="rbind")

# 2) create a cbind with the sample score list
bind.fill(score_list, type="cbind")

# 3) create a cbind with the sample score list,
#    and fill missing data with 0s.
bind.fill(score_list, type="cbind", fill = 0L)

}
\author{
Hwanggyu Lim \email{hglim83@gmail.com}
}
