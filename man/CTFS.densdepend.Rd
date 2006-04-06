\name{CTFS.densdepend}
\alias{CTFS.densdepend}
\concept{density dependence}
\concept{neighborhood}
\concept{dynamics}

\title{Counting Neighbors and their Attributes - Investigating Density 
Dependence}

	\description{ In order to identify the presence of density dependence for
	a focal species, attributes of the neighborhood surrounding the individuals
	of that species need to be determined.  In the CTFS package, this is done
	using \code{\link{neighbordens.1sp}}, described in greater detail
	below.  The result of this function is a \strong{dataframe} containing measures of
	neighborhood attributes for each defined neighborhood for each tree of a
	focal species.  This file can be used to analyse the effect of
	neighborhood conditions on the size, growth, mortality and recruitment of
	individuals for each species by combining it with the data for each tree 
	of the focal species. }
	
\details{

PREPARATIONS

	In order to run these functions, a dataset structured by quadrate needs
	to be created.  It must be a \strong{list} containing \strong{dataframes}
	for each quadrate which, in turn, have the data for a single census for
	all the trees in each quadrate.  This \strong{list} is created from a
	CTFS census file.  Use \code{\link{sep.quadinfo}} first to create a
	vector containing the quadrate number for each tree.  Then
	\code{\link{sep.quadinfo}}, using this vector, will create an
	appropriately structured \strong{list}.  The example data sets are: \code{\link{tst.bci90.quad}} and
	\code{\link{tst.bci95.quad}}.
	

FUNCTIONS FOR COMPUTING NEIGHBORHOOD ATTRIBUTES

	\code{\link{neighbordens.1sp}}
	
	Cover function that returns a \strong{dataframe} with neighborhood attributes for
	individuals of a single species (focal species).  This function can be
	called from other functions to produce \strong{dataframes} for many species.
	
	\code{\link{neighbordens.1quad} }
	
	Called by \code{neighborhdens.1sp} for each quadrate in the plot.
	Calls \code{findborderquads} to locate surrounding quadrates that
	constitute the area of the neighborhoods.  This limits the searching for
	neighbors to only surrounding quadrates within the radius of the
	neighborhood.  Then it calls \code{countdens.1quad} to compute
	neighborhood attributes for individuals of the focal species in a
	single quad.
	
	\code{\link{countdens.1quad}}

	Called by \code{neighbordens.1quad} for each quadrate using the subset of
	the data that contains the trees in the neighboring quadrates and the
	trees of the focal species.  Using \code{borderdist}, the neighborhood
	attributes are set to NA for trees that are sufficiently close to a
	border of the plot such that for a given neighborhood distance, that
	neighborhood is not entirely within the plot.
		
	\code{\link{countdens.1tree}}

	This function actually computes the number and basal area of 
	conspecifics and heterospecifics in a given neighborhood for a given 
	individual tree of the focal species.
	
	LIST OF FUNCTIONS CALLED BY USER FUNCTIONS 
	
	\code{\link{findborderquads}} identified bordering quadrates\cr
	
	\code{\link{xydistvect}} computes xy distance between points for vectors
	of points\cr
	
	\code{\link{fill.distclasses}} fills a matrix of values, structured for
	density dependent functions\cr
	
	\code{\link{borderdist}} computes distance to the plot border for a given
	coordinate pair\cr


NEIGHBORHOOD ATTRIBUTES COMPUTED AND HOW TO COMPUTE OTHER ATTRIBUTES
	
	Number of conspecifics, \code{Ncon} and sum of the basal area of all 
	conspecifics, \code{BAcon}
	Number of heterospecifics, \code{Nhet} and sum of the basal area of all 
	heterospecifics, \code{BAhet}
	
	If the user wishes to compute other attributes of neighbors or select 
	neighbors on other basis than being alive (or alive and broken), then 
	the user should write new versions of \code{countdens.1tree} and/or 
	\code{countdens.1tree.1quad}.  The tree data for the focal species is 
	contained in the variable \code{focaldata} and the data for the trees in the 
	neighboring quadrates that could be valid neighbors is in the variable 
	\code{neighdata}.  These are both \strong{dataframes}.
	
	Changing the types of attributes computed from all living neighbors
	involves rewriting \code{\link{countdens.1tree}}.  At this time
	\code{countdens.1tree} is passed only the \code{dbh} of a neighbor tree.  If
	the desired attribute can be computed from that, then only
	\code{countdens.1tree} needs to be rewritten.  If other information about the
	neighbor tree is needed, then the user needs to rewrite \code{countdens.1tree}
	and \code{countdens.1tree.1quad} so that the calling statements provide the
	needed variables for the user defined neighbor attribute.
	
	If the user wishes to select neighbors based on other attributes than
	species name and status of alive (or alive but broken), it is suggested
	that the user rewrite \code{countdens.1tree.1quad} as this is where the
	subsetting of the focal species tree data and the neighbor tree data take
	place.  The logical vector \code{inc} defines which neighbor trees to
	use.  Redefining this vector would, most likely, be the simpliest
	adjustment. There is a comment in the code at the location of the 
	definition of \code{inc}.
	
	
}
\author{Rick Condit and Pamela Hall}

\seealso{\code{\link{CTFS.mortality}}, \code{\link{CTFS.growth}}, 
\code{\link{CTFS.recruitment}}, \code{|link{CTFS-internal}}
}

\examples{\dontrun{
data(tst.bci90.quad)
data(tst.bci90.spp)
socrex90.neigh <- neighbordens.1sp("socrex",tst.bci90.quad, tst.bci90.spp)
str(socrex90.neigh)
plot(socrex90.neigh$Ncon.0.5,tst.bci90.spp$socrex$dbh)
}}

\keyword{univar}

