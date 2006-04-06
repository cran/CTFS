\name{abundance.quad}
\alias{abundance.quad}
\concept{dynamics}
\concept{abundance}
\concept{population size}

\title{Abundance of Trees by Species and User Defined Quadrates}

\description{
 Computes the number of trees for each species in user specified
 quadrates.  The user can also specify any portion of the plot (or all
 of it) be considered.  The default quadrate size is 1 hectare.
}
\usage{
abundance.quad(census1, mindbh = 10, 
	plotdim = c(1000, 500), gridsize = 100)
}

\arguments{
  \item{census1}{name of census datafile for a single census}
  \item{mindbh}{Minimum DBH for computing population size in each 
  census}
	\item{plotdim}{dimensions of the plot: east-west 1000m and
	north-south 500m}
  \item{gridsize}{side of the square quadrate, 100 x 100m (1 ha) by
	default}
}
\details{
	This function produces the same results as running \code{\link{abundance}} using
	the split vars of species and quadrate.
	
	See \code{\link{CTFS.abundance}} for details on the computation methods
	of abundance and associated functions.  Values for only 1 census can be
	computed at a time.  A tree is included for the computation based on its
	value for \code{status} and \code{dbh} for a single census.
	
}

\value{
  \code{abundance.quad} returns a \strong{list of arrays} with the values
  of \strong{species} as the first dimension and the values of
  \strong{gridsize} as the second dimension of the array.  The array
  contains the following named components:
  
  \item{$N}{the abundance of each category, population size}
  \item{$dbhmean}{the mean dbh in mm for trees used to compute 
  abundance}
  \item{$meandate}{mean date of census for trees used to compute abundance}
  
  The first \strong{array} of the returned \strong{list} is a matrix
  of the abundance of each \strong{species} by \strong{gridsize}.
  Abundance values are returned for all levels of each vector and if
  no value can be computed then 0 or NA is returned as appropriate.
}


\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.abundance}}
}
\examples{\dontrun{
1. Default use of abundance.quad
abund.ha.out <- abundance.quad(tst.bci90.full)
abund.ha.out[[1]]

2. Compute abundance for each quadrate.  If some quadrates are missing
from the input data file, then they will also be missing from the
output.
abund.quad.out <- abundance.quad(tst.bci90.full,gridsize=20)
abund.quad.out[[1]][,1:10]
}
}
\keyword{univar}
