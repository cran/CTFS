\name{CTFS.growth}
\alias{CTFS.growth}
\concept{dynamics}
\concept{growth}
\concept{relative growth}
\title{Overview of Analysis of Growth Rates in Tropical Tree Populations}
\description{
 Provides an overview for the analysis of growth rates for tree
 populations by a variety of categories including the functions that 
 are available, supporting functions and ways for using the options of 
 the main functions.
}
\details{
	FUNCTIONS TO COMPUTE GROWTH RATES
	
	\tabular{ll}{
	\code{\link{growth}} \tab Annual Growth Rates by Categories (User 
	defined groups)\cr
	\code{\link{growth.dbh}} \tab Annual Growth Rates by DBH Classes\cr
	\code{\link{growth.eachspp}} \tab Annual Growth Rates by Species\cr
	\code{\link{growth.indiv}} \tab Annual Growth Rates for each Tree\cr
	}
	
	FUNCTIONS FOR FORMATTING RESULTS
	
	\tabular{ll}{
	\code{\link{assemble.demography}} \tab Reformat the Output from
	Demographic Functions from List to Dataframe\cr
	\code{\link{assemble.table.each}} \tab Reformat the Output from
	Demographic Functions from List to Dataframe\cr
	}
     
	LIST OF FUNCTIONS CALLED BY USER FUNCTIONS

	\tabular{ll}{
	\code{\link{rndown5}} \tab Rounds down DBH where measurement was less 
	precise than 1 mm\cr
	\code{\link{trim.growth}} \tab Removes extreme values of growth using user 
	defined limits\cr
	\code{\link{fill.dimension}} \tab Fills all the dimensions of a 2 dimensional array\cr
	\code{\link{fill.1dimension}} \tab Fills all the dimensions of a 1 dimensional array\cr
	}
	
	COMPUTATION OF GROWTH
	
	The annual growth rate for a single tree is computed as:
	\deqn{growth = \frac{( dbh1 - dbh0 )}{(time1 - time0)}}{%
	growth = ( dbh1 - dbh0 ) / (time1 - time0)}
	
	where \code{dbh0} and \code{dbh1} are the diameters of the tree at the 
	initial and subsequent census, respectively;
	\cr 
	where \code{time0} and \code{time1} are the census dates for the initial 
	and subsequent census, respectively expressed in years;
	\cr 
	\cr The mean growth rate, \code{meangrow}, is the mean of the growth rate for all trees 
	in the population as defined by the provided file.  The mean growth 
	rate is expressed in mm per year.
	\cr
	\cr The relative growth rate for a single tree is computed as:
	\deqn{RGR  = \frac{( log(dbh1) - log(dbh0) )}{(time1 - time0)}}{%
	RGR  = ( logdbh1 - logdbh0 ) / (time1 - time0)}
	
	with variables as described above.
	\cr 
	\cr The mean relative growth rate, \code{meanRGR}, is the mean of the relative growth 
	rate for all trees in the population as defined by the provided file.  
	The mean relative growth rate is expressed as 1 / years.
	\cr
	\cr The \code{meandbh} and \code{meantime} is the mean of the dbh and census 
	intervals, respectively, for the trees used for the mean growth 
	rate computation.
	\cr
	\cr Before computing the growth rate, the ``status'' of a tree is 
	determined.  The tree must have been recorded as ``A'', alive at both 
	censuses. Therefore, recruits and trees that died are not used for 
	computing growth rates.  Also, the time of the census must be 
	recorded for both censuses.
	\cr
	\cr The 95\% confidence intervals are computed using a 
	t-distribution.
}
	

\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.groupcats}}}
\keyword{univar}
