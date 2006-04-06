\name{mortality.calculation}
\alias{mortality.calculation}
\concept{dynamics}
\concept{mortality}
\concept{confidence limits}

\title{ Compute Annual Mortality Rates and Confidence Limits given N, S and Time }
\description{
A simple function to calculate mortality and confidence limits from N0
(number alive at the outset), S (number of survivors), and time (time
interval).  Values can be vectors or scalars.
}
\usage{
mortality.calculation(N0, S, meantimeint)
}

\arguments{
  \item{N0}{number of trees at the first census}
  \item{S}{number of trees that survived to the second census, N - D}
  \item{meantimeint}{interval between censuses in years for trees in the 
	first census}
}

\details{
	This function is called many times by \code{\link{mortality}}.  It can be
	used directly but it is primarily intended as an internal
	function. \code{N0, S} and \code{meantimeint} can be single values or 
	vectors.  For full details on computation of mortality rate and 
	confidence limits see \code{\link{CTFS.mortality}}.
}	
\value{
  \code{mortality.calculation} returns a \strong{dataframe} with the
  following components:

  \item{N0}{number of trees at the first census}
  \item{S}{number of trees that survived to the second census, N - D}
  \item{D}{number of dead trees}
  \item{rate}{annual mortality rate}
  \item{lowerCI}{lower 95\% confidence interval}
  \item{upperCI}{upper 95\% confidence interval}
  \item{meantimeint}{mean number of years between intervals}
}
\author{Rick Condit, Pamela Hall}
\seealso{\code{\link{CTFS.mortality}}}
}

\examples{\dontrun{
# 1. Default use of mortality.calculation()

mortality.calculation(500, 490, 5)

# N0   S  D        rate     lowerCI     upperCI meantimeint
# 1 500 490 10 0.004040541 0.002214243 0.007415848           5

# 2. Input values as vectors
N0 = c(100, 100, 100)
S = c(90, 80, 70)
time = c(5, 5, 5)
mortality.calculation(N0,S,time)

# N0  S  D       rate    lowerCI    upperCI meantimeint
# 1 100 90 10 0.02107210 0.01144898 0.03836600           5
# 2 100 80 20 0.04462871 0.02867126 0.06825114           5
# 3 100 70 30 0.07133499 0.04943058 0.10088494           5
}}

\keyword{univar}
