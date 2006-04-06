\name{CTFS.abundance}
\alias{CTFS.abundance}
\concept{dynamics}
\concept{abundance}
\concept{population size}

\title{Overview of Analysis of Abundance and Rate of Change in Population 
Size in Tropical Tree Populations}

\description{
 Provides an overview for the analysis of abundance and change in 
 population size by a variety of categories.  This overview includes functions that 
 are available, supporting functions and ways for using the options of 
 the main functions.
}

\details{

	FUNCTIONS TO COMPUTE ABUNDANCE AND CHANGE IN POPULATION SIZE 
RATE
	
	\tabular{ll}{
	\code{\link{abundance}} \tab Abundance of Tree by Categories (User Defined Groups)\cr
	\code{\link{abundance.spp}} \tab Abundance of Trees by Species and DBH 
	Classes, DBH determined separately for each year\cr
	\code{\link{abundance.quad}} \tab Abundance of Trees by Species and User Defined Quadrates\cr
	\code{\link{abundance.change}} \tab Annual Change in Population Size by Categories (User Defined
Groups\cr	
	}

	FUNCTIONS FOR FORMATTING RESULTS
	
	\tabular{ll}{
	\code{\link{assemble.demography}} \tab Reformat the Output from
	Demographic Functions from List to Dataframe\cr
	}

	COMPUTATION OF ABUNDANCE AND CHANGE IN ABUNDANCE\cr
	
	Abundance is the number of individuals alive at a given census.  
	Trees that are broken at the point of dbh measurement are not 
	included.  By convention, the first census is referred to as \code{N0} 
	and the second census as \code{N1}, as in an initial and any subsequent 
	census.  Censuses do NOT have to be in immediate order though they must 
	be in the order of \dQuote{earlier, more recent}.  The change in 
	abundance can be computed for any two censuses no matter how many years 
	apart they are.
	
	The annual population rate is calculated as
  \deqn{p = \frac{logN1 - logN0)}{(time1 - time0)}}{%
	 p = (logN1-logN0) / (time1-time0)}
	 
	where \code{N0} is the number of live individuals at the first census,\cr 
	where \code{N1} is the number of live individuals at the	second census 
	which includes recruits at that census,\cr 
	where \code{time1} and \code{time0} are expressed in years.
  
	The probabiilty that \eqn{N0 = N1} (test of hypothesis) is computed as 
	the combined probability of the observed mortality and recruitment.  
	Probabilities for mortality are computed using a binomial distribution 
	and a poisson distribution is used for recruitment.  The number of 
	survivors, \code{S}, is computed based \code{N0}.  

	DETERMINATION OF WHETHER TO INCLUDE A TREE IN THE COMPUTATION
	
	By default, trees are included in the computation of abundance for a
	single census if their \code{status} = \dQuote{A} and if their
	\code{mindbh} >= 10 mm for the census being used.  The user can alter these
	values.
	
	For computation of changes in abundance, the inclusion of a tree is
	determined independently for each census.  The user can only define
	\code{status} and \code{mindbh} once.  These values are used for
	assessing a tree for inclusion for each census.
	
	Survivors are trees that were alive in census 1 and are also alive in 
	census 2.
}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.groupcats}}}
\keyword{univar}
