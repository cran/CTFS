\name{findborderquads}
\alias{findborderquads}
\concept{neighbors}
\concept{border}
\concept{gx}
\concept{gy}
\concept{density dependent}

\title{Returns Indices of Neighboring Quadrats Given a Quadrat Index}
\description{
  Identified neighboring quadrates from a given quadrate number,
  refered to as the \code{index}.  Returns the quadrat numbers within a given distance.}
\usage{
findborderquads(index, dist = 20, gridsize = 20, plotdim = c(1000, 500))
}
\arguments{
  \item{index}{quadrate number, between 0 and 1249 in the standard
  plot}
  \item{dist}{distance in m within which the neighboring quadrates are
  located.  Distance is measured from any side of the index quadrate.}
 	\item{gridsize}{side of the square quadrate, 20 x 20m by default}
	\item{plotdim}{dimensions of the plot: east-west 1000m and
	north-south 500m}
}
\details{
   A quadrate is considered a neighbor if
  any part of it lies within \code{dist} of the center of the quadrate
  indicated by \code{index}.  At the minimum the quadrates that share 
  edges, including corners, will always be returned.  This is 8
  quadrates.  The \code{index} quadrate is NOT returned.  All values
  >=0 and <20 give the same 8 quadrates for a centrally located
  quadrate.
  
  As \code{dist} increases, more quadrates can be included up to the
  edges of the plot.  For \code{gridsize}=20, 8 quadrates are
  neighbors up to \code{dist}<20.  This doesn't change until
  \code{dist} >= 20, then there are 24 neighbors, etc.
  
  This function depends upon \code{\link{index.to.rowcol}} to identify
  the neighboring quadrates. 
}
\value{
	Returns a vector of numbers, the quadrate indices for all
	surrounding quadrates.}

\author{Rick Condit}
\seealso{\code{\link{CTFS.quadfunctions}}}
\examples{\dontrun{
sort(findborderquads(231,dist=10))
sort(findborderquads(231,dist=20))
sort(findborderquads(231,dist=0))

}
}
\keyword{manip}
