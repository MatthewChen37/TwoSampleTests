% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests.R
\name{EBT}
\alias{EBT}
\title{EBT Function}
\usage{
EBT(dst, ix, sizes)
}
\arguments{
\item{dst}{is the distance matrix of the data set}

\item{ix}{is a permutation of row indices of z}

\item{sizes}{is a vector of the sample sizes}
}
\value{
The EBT test statistic
}
\description{
Calculates the test statistic for the two-sampled Energy Distance Test
}
\examples{
x = matrix(rexp(100 * 2), 100, 2)
y = matrix(rexp(200 * 2), 200, 2)
z <- rbind(x,y)
dst <- as.matrix(dist(z))
EBT_statistic <- EBT(dst, 1:nrow(dst), c(100, 200))
}
