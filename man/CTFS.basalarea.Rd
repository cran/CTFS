\name{CTFS.basalarea}
\alias{CTFS.basalarea}
\concept{dynamics}
\concept{basal area}

\title{Overview of Analysis of the Change in Basal Area of Tropical Tree Populations}

\description{
 Provides an overview for the analysis of basalarea and change in
 basal area for tree grouped into a variety of categories.  This
 overview includes functions that are available, supporting functions
 and ways for using the options of the main functions.
 
 Basal area is computed on the tagged tree from a census file and from
 the sum of any multiple stems of the tagged tree that exist.  The
 census and mulit-stems datasets are used for these analyses.  If
 the multi-stem datasets do not exist, basal area is computed only
 from the dbh of the tagged tree.
}

\details{

	FUNCTIONS TO COMPUTE BASAL AREA AND CHANGE IN BASAL AREA
	
	\tabular{ll}{
	\code{\link{basal.area}} \tab Basal Area of Tree by Categories (User Defined Groups)\cr
	\code{\link{ba.indiv}} \tab Basal Area of Individual Trees\cr
	\code{\link{ba.change}} \tab Annual Change in Basal Area by
	Categoreis (User Defined Groups)\cr
	\code{\link{convert.ba.mult}} \tab Computes Basal Area of all Stems
	of a Tagged Tree in the Multi-Stems dataset\cr
	\code{\link{ba}} \tab Computes Basal Area \eqn{(m^2)} from a
	dbh vector\cr
	\code{\link{basum}} \tab Computes the sum of the Basal Area from a 
	dbh vector\cr
	}

	FUNCTIONS FOR FORMATTING RESULTS
	
	\tabular{ll}{
	\code{\link{assemble.demography}} \tab Reformat the Output from
	Demographic Functions from List to Dataframe\cr
	}

	COMPUTATION OF BASAL AREA AND CHANGE IN BASAL AREA\cr
	
	Basal area for an individual is the sum of the basal area of the main
	stem and the basal area of each living multiple stem.  Trees that are
	alive in any given census and are larger than the mindbh in a census are
	used for the computation.  This means that a tree can be included in one
	census and not the other if its dbh is below the minimum in one census
	and not the other.  For instance, trees that are broken at the point of
	dbh measurement are not included in the census when they are broken.
	
	The basal area, in \eqn{m^2} is calculated as
    \deqn{sumba = \sum(pi * ( \frac{dbh[i]}{2000} )^2)}{%
		sumba = sum(pi * (dbh[i] / 2000)^2 )}
    
	where \code{dbh[i]} is the dbh in mm of each stem of an individual 
	tree
    \cr 

	The change in basal area is calculated as
	\deqn{bachange = \frac{(\sum(logba2)-\sum(logba1))}{(mean(time2-time1))}}{%
	bachange = ( sum(logba2) - sum(logba1) ) / ( mean(time2-time1) )}
	
	where \code{sum(logba2)} is the sum of the basal area of all the
	individuals alive at the second census
    \cr 
	where \code{sum(logba1)} is the sum of the basal area of all the
	individuals alive at the first census
    \cr 
	where \code{time1} and \code{time0} are expressed in years.

	DETERIMINATION OF WHETHER TO INCLUDE A TREE IN THE COMPUTATION
	
	By default, trees are included in the computation of basal area for a
	single census if their \code{status} = \dQuote{A} and if their
	\code{mindbh} >= 10 for the census being used.  The user can alter these
	values.
	
	For computation of changes in basal area, the inclusion of a tree is
	determined independently for each census.  The user can only define
	\code{status} and \code{mindbh} once.  These values are used for
	assessing a tree for inclusion for each census.


}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.groupcats}}}
\keyword{univar}
