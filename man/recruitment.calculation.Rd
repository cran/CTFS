\name{recruitment.calculation}
\alias{recruitment.calculation}
\concept{dynamics}
\concept{mortality}
\concept{confidence limits}

\title{Compute Annual Recruitment Rates and Confidence Limits given N, S and Time}
\description{
A simple function to calculate recruitment and confidence limits from N1
(number alive at the end of interval), S (number of survivors), and time (time
interval).  Values can be vectors or scalars.
}
\usage{
recruitment.calculation(N1, S, meantimeint)
}

\arguments{
  \item{N1}{number of trees at the first census}
  \item{S}{number of trees that survived to the second census, N - D}
  \item{meantimeint}{interval between censuses in years for trees in the 
	first census}
}
\details{
	This function is called many times by \code{\link{recruitment}}.  It can be
	used directly but it is primarily intended as an internal
	function. \code{N1, S} and \code{meantimeint} can be single values or 
	vectors.  For full details on computation of recruitment rate and 
	confidence limits see \code{\link{CTFS.recruitment}}.
}	
\value{
  \code{recruitment.calculation} returns a \strong{dataframe} with the
  following components:

  \item{N1}{number of trees at the first census}
  \item{S}{number of trees that survived to the second census, N - D}
  \item{R}{number of recruited trees}
  \item{rate}{annual recruitment rate}
  \item{lowerCI}{lower 95\% confidence interval}
  \item{upperCI}{upper 95\% confidence interval}
	\item{meantimeint}{mean number of years between intervals}
}
\author{Rick Condit, Pamela Hall}
\seealso{\code{\link{CTFS.recruitment}}}
}
\examples{\dontrun{
1. Default use of rexcruitment.calculation()

recruitment.calculation(510, 490, 5)
   N1   S  R        rate     lowerCI    upperCI meantimeint
1 510 490 20 0.008001067 0.005189957 0.01233286           5

2. Input values as vectors
N = c(110, 110, 110)
S = c(90, 80, 70)
time = c(5, 5, 5)
recruitment.calculation(N,S,time)
   N1  S  R       rate    lowerCI    upperCI meantimeint
1 110 90 20 0.04013414 0.02581525 0.06143132           5
2 110 80 30 0.06369075 0.04420141 0.09014424           5
3 110 70 40 0.09039702 0.06561979 0.12211671           5
}}

\keyword{univar}
