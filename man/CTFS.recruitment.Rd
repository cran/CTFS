\name{CTFS.recruitment}
\alias{CTFS.recruitment}
\concept{dynamics}
\concept{recruitment}

\title{Overview of Analysis of Recruitment Rates in Tropical Tree Populations}
\description{
 Provides an overview for the analysis of recruitment rates for tree
 populations by a variety of categories including the functions that 
 are available, supporting functions and ways for using the options of 
 the main functions.
}
\details{
	FUNCTIONS TO COMPUTE RECRUITMENT RATES
	
	\tabular{ll}{
	\code{\link{recruitment}} \tab Annual Recruitment Rates by Categories (User 
	defined groups)\cr
	\code{\link{recruitment.eachspp}} \tab Annual Recruitment Rates by Species\cr
	}
	
	FUNCTIONS FOR FORMATTING RESULTS
	
	\tabular{ll}{
	\code{\link{assemble.demography}} \tab Reformat the Output from
	Demographic Functions from List to Dataframe\cr
	}
    
	FUNCTIONS CALLED BY USER FUNCTIONS\cr

	\tabular{ll}{
	\code{\link{find.climits}} \tab Calculates confidences limits for 
	recruitment rates\cr
	\code{\link{fill.dimension}} \tab Fills all the dimensions of a 2 dimensional array\cr
	\code{\link{fill.1dimension}} \tab Fills all the dimensions of a 1 dimensional array\cr
}

	COMPUTATION OF RECRUITMENT
	
	The annual recruitment rate is calculated as
  \deqn{r = \frac{(logN1 - logS)}{mean(time1 - time0)}}{%
	r = (logN1-logS) / mean(time1-time0)}
	
	where \code{N1} is the number of live individuals at the second census,\cr 
	where \code{S} is the number of surviving individuals: trees that were 
	alive in the first census and alive in the second census,\cr 
	where \code{time1} and \code{time0} are expressed in years and is the 
	length of the census interval as recorded for trees in the second census (which 
	includes all of the surviving trees.
  
	The confidence limits are computed using \code{\link{find.climits}} which
	returns the number of survivors, \code{S}, out of \code{N} individuals for
	each confidenc interval in turn.  The beta distribution is used to
	determine the number of S for the upper 95\% and lower 5\% (default
	probability level).  Confidence limits for mortality rate are computed as
	from these CI for \code{S} as:
	\deqn{rate.CI.upper = \frac{( logN1 - log(S.lowerCI) )}{mean(time1 - time0)}\cr
	rate.CI.lower = \frac{( logN1 - log(S.upperCI) )}{mean(time1 - time0)}}{%
	rate.CI.upper = ( logN1 - log(S.lowerCI) ) / mean(time1 - time0)\cr
	rate.CI.lower = ( logN1 - log(S.upperCI) ) / mean(time1 - time0)}
	
	Note that \code{S.lowerCI} is a lower value of \code{S} which results in a higher
	number of recruits and hence a higher recruitment rate.  And \code{S.upperCI} is a
	higher value of \code{S} which results in a lower number of deaths and hence a
	lower recruitment rate.
  }

\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.groupcats}}}
\keyword{univar}
