\name{abundance.spp}
\alias{abundance.spp}
\concept{population size}
\concept{species density}

\title{Abundance of Trees by Species and DBH Classes}
\description{
  Computes the number of trees for each species and each DBH class.
  DBH classes can be user defined.  The dataset may have 1 or more
  censuses but only 1 census is used by this function.  See details
  for how the variables should be named.
}
\usage{
abundance.spp(census1, alivecode=c("A"), dbhclass=c(10,100,300))
}

\arguments{
\item{census1}{name of census datafile for a single census}
\item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  Valid categories are: "A" and 
  "AB" and "AS".}
\item{dbhclass}{lower class boundaries for DBH classes, classes
   are nonoverlapping}
}
\details{
	This function returns the same values as \code{\list{abundance}} using
	species and dbh classes as the split variables. 

	See \code{\link{CTFS.abundance}} for details on the computation methods 
	of abundance and associated functions. Values for only 1 census
	can be computed at a time.
	
	A tree is included for the computation based on its value for
	\code{status} and \code{dbh} for a single census. 

}
\value{
  \code{abundance.spp}returns a \strong{list of arrays} with the values
  of \strong{species} as the first dimension and the values of
  \strong{dbh category} as the second dimension of the array.  The array
  contains the following named components:
  
  \item{$N}{the abundance of each category, population size}
  \item{$dbhmean}{the mean dbh in mm for trees used to compute 
  abundance}
  \item{$meandate}{mean date of census for trees used to compute abundance}
  
  The first \strong{array} of the returned \strong{list} is a matrix
  of the abundance of each \strong{species} by \strong{dbh class}.
  Abundance values are returned for all levels of each vector and if
  no value can be computed then 0 or NA is returned as appropriate.
}


\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.abundance}}
}
\examples{\dontrun{
# 1. Default use of abundance.spp
abund.sp.out <- abundance.spp(tst.bci90.full)
abund.sp.out[[1]]

# 2. Using different dbh classes
abund.sp.out1 <-
abundance.spp(tst.bci90.full,dbhclass=c(100,200,300,400,500))

}
}
\keyword{univar}
