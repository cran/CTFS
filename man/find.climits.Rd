\name{find.climits}
\alias{find.climits}
\concept{dynamics}
\concept{mortality}
\concept{recruitment}
\concept{confidence limits}

\title{Computes 95\% CI for a Survivorship}
\description{ Computes the confidence intervals for the number of 
surviving trees.  This can be used to compute confidence intervals for 
rates of mortality and recruitment.  The intervals are based on the beta 
distribution.
}

\usage{
find.climits(N, S, alpha = 0.05, kind = "upper")
}
\arguments{
  \item{N}{number of individuals from a census}
  \item{S}{number of survivors}
  \item{alpha}{type I error }
  \item{kind}{type of confidence limit}
}
\details{
	This function is called from \code{\link{mortality.calculation}} and
	\code{\link{recruitment.calculation}}.  On its own it computes the number
	of surviving trees for one confidence interval of a given \option{kind}
	for given values of \option{N}.  For mortality, \option{N} should be
	\eqn{N0} the number of individuals alive in the first census.  For
	recruitment, \option{N} should be \eqn{N1} the number of individuals
	alive in the second census.  For both mortality and recruitment,
	\option{S} is the number of individuals that were alive in \eqn{N0} AND
	in \eqn{N1}.  Since the confidence intervals provided are limits are
	survivors, the lower confidence interval for survivorship provides the
	upper confidence interval for mortality rate (high survivorship = lower
	mortality) and upper confidence interval for recruitment rate (high
	survivorship = lower number of recruits).  And the upper confidence
	interval for survivorship provides the lower confidence for mortality
	(low survivorship = higher mortality) and the lower confidence interval
	for recruitment rate (low survivorship = higher number of recruits).
	
	The confidence limits for survivorship are computed using a beta 
	distribution.
	}

\value{
  Returns the confidence limits in number of surviving trees:
  
  \item{S }{numeric, the best estimate of the number of surviving trees for 
  the given \option{kind} of limit.}
}

\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.mortality}}, \code{\link{CTFS.recruitment}}}

\examples{\dontrun{
# N0 = 100, Dead = 5, Recruits = 10, time = 5 yrs, therefore: 
# S = N0 - Dead = 100 - 5 = 95
# N1 = N0 - Dead + Recruits = 100 - 5 + 10 = 105

# find CI limits for mortality
SM.upper <- find.climits(100,95,kind="upper")
SM.lower <- find.climits(100,95,kind="lower")
mort.rate.upper <- ( log(N0) - log(SM.lower) ) / 5
mort.rate.lower <- ( log(N0) - log(SM.upper) ) / 5

# find CI limits for recruitment
SR.upper <- find.climits(105,95,kind="upper")
SR.lower <- find.climits(105,95,kind="lower")
rec.rate.upper <- ( log(N1) - log(SR.lower) ) / 5
rec.rate.lower <- ( log(N1) - log(SR.upper) ) / 5

}}

\keyword{univar}


