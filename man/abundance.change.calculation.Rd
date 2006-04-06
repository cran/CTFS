\name{abundance.change.calculation}
\alias{abundance.change.calculation}
\concept{dynamics}
\concept{abundance}
\concept{population size}

\title{Calculation of Annual Change in Population Size and Probability 
that Population Size has Changed}
\description{
  Computes annual change in population size (rate) and probability that 
	that rate is = 0 (prob) for values provided.  It is 
	called from \code{\link{abundance.change}} but can also be used 
	independently with scalar, vector or 2 dimensional arrays.
}
\usage{
abundance.change.calculation(N0 = N0, N1 = N1, S = S, meantimeint = meantimeint)
}
\arguments{
  \item{N0}{abundance in first census}
  \item{N1}{abundance in second census}     
	\item{S}{number of surviors from first census}
	\item{meantimeint}{mean number of years between census for trees used in 
  abundance change computation}
}
\details{
	See \code{\link{CTFS.abundance}} for details on the computation of 
	abundance and associated functions.

	If arguments are not arrays, they are coerced into arrays and any names
	associated with their dimensions are retained.  This appears to be quite
	robust regardless of the object form of the arguments.  But if problems
	asired, then make arguments into arrays with dimnames.  Alternatively,
	call this from this from \code{\link{abundance.change}}.
	}
\value{
  Returns a list of:
  \item{rate }{annual rate of change in abundance}
  \item{prob }{probability that N0 == N1}
}
\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.abundance}}, \code{\link{abundance.change}}
}
\examples{\dontrun{
1.  use with scalar
abundance.change.calculation(500,510,490,5)

2.  use with vectors that have no names associated
N0 <- seq(100,1000,by=100)
N1 <- seq(110,1100,by=100)
S <- seq(90,990,by=100)
meantimeint <- rep(5,10)

}}
\keyword{htest}
