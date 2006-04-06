\name{CTFS.mortality}
\alias{CTFS.mortality}
\concept{dynamics}
\concept{mortality}

\title{Overview of Analysis of Mortality Rates in Tropical Tree Populations}

\description{
 Provides an overview for the analysis of mortality rates for tree
 populations by a variety of categories.  Overview includes the functions that 
 are available, supporting functions and use of the options of 
 the main functions.
}
\details{
	FUNCTIONS TO COMPUTE MORTALITY RATES
	
	\tabular{ll}{
	\code{\link{mortality}} \tab Annual Mortality Rates by Categories (User 
	defined groups)\cr
	\code{\link{mortality.dbh}} \tab Annual Mortality Rates by DBH Classes\cr
	\code{\link{mortality.eachspp}} \tab Annual Mortality Rate for Species by DBH Class\cr
	\code{\link{mortality.calculation}} \tab  A Single Annual Mortality Rate given N, S and Time \cr
	}
	
	FUNCTIONS FOR FORMATTING RESULTS
	
	\tabular{ll}{
	\code{\link{assemble.demography}} \tab Reformat the Output from
	Demographic Functions from List to Dataframe\cr
	\code{\link{assemble.table.each}} \tab Reformat the Output from
	Demographic Functions from List to Dataframe\cr
	}
    
	FUNCTIONS CALLED BY MORTALITY FUNCTIONS

	\tabular{ll}{
	\code{\link{find.climits}} \tab Calculates confidences limits for 
	mortality rates\cr
	\code{\link{fill.dimension}} \tab Fills all the dimensions of a 2 dimensional array\cr
	\code{\link{fill.1dimension}} \tab Fills all the dimensions of a 1 dimensional array\cr
	}

	COMPUTATION OF MORTALITY
	
	The annual mortality rate is calculated as 
	\deqn{m = \frac{(logN0-logS)}{ mean(time1-time0)}}{%
	m = (logN0 - logS) / mean(time1 - time0)}
	
	where \code{N0} is the number of live individuals at the first census, \cr
	where \code{S} is the number of surviving individuals of N0 at the second
	census which is equal to \eqn{N0 - D}, \cr 
	where \code{time1} and \code{time0} are expressed in years.
  
	The confidence limits are computed using \code{\link{find.climits}} which
	returns the number of survivors, \code{S}, out of \code{N} individuals for
	each confidence interval in turn.  The beta distribution is used to
	determine the number of S for the upper 95\% and lower 5\% (default
	probability level).  Confidence limits for mortality rate are computed as
	from these CI for \code{S} as:
	\deqn{	rate.CI.upper = \frac{( logN0 - log(S.lowerCI) )}{mean(time1 - time0)}\cr
	rate.CI.lower = \frac{( logN0 - log(S.upperCI) )}{mean(time1 - time0)}}{%

	rate.CI.upper = ( logN0 - log(S.lowerCI) ) / mean(time1 - time0)\cr
	rate.CI.lower = ( logN0 - log(S.upperCI) ) / mean(time1 - time0)}
	
	Note that \code{S.lowerCI} is a lower value of \code{S} which results in a higher
	number of deaths and hence a higher mortality rate.  And \code{S.upperCI} is a
	higher value of \code{S} which results in a lower number of deaths and hence a
	lower mortality rate.
  
	}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.groupcats}}}
\keyword{univar}
