\name{neighbordens.1sp}
\alias{neighbordens.1sp}
\alias{neighbordens.1quad}
\alias{countdens.1quad}
\alias{countdens.1tree}
\concept{density dependence}
\concept{neighborhood}

	\title{Create File of Neighbor Attributes for a Single Species}
	\description{ 
	Calculates the number and summed basal area of conspecific
	and heterospecific trees within defined distances of all individuals of a
	given species.  The resulting data frame can be used for analysis of the
	effect of neighbors on the population characteristics of a species.  }
	
	\usage{neighbordens.1sp(spname, census.quad, census.spp, distbreak =
	c(0, 5, 10, 20, 30), gridsize = 20, plotdim = c(1000, 500)) }

\arguments{
  \item{spname}{species code name (aka target species), in quotes}
  \item{census.quad}{data file as a list with quadrates as names}
  \item{census.spp}{data file as a list with species as names}
  \item{distbreak}{vector of radii for neighborhood annuli}
  \item{gridsize}{gridsize that \code{census.quad} was created with}
  \item{plotdim}{dimensions of plot}
	} 
	
\details{ 

	This function finds all trees within a given distance for each individual
	tree of the target species.  It computes the number of individuals and
	their total basal area for two groups: conspecifics and heterospecifics.
	This is done for each neighborhood size.  These values are returned for
	each individuals of the target species in a data frame with the 
	identification number of the target species trees.
	
	The function requires a two versions of the tree data file.  One is a
	\strong{list of data frames} of the trees in each quadrate.  The second
	is a \strong{list of data frames} of trees in each species.  See
	\code{\link{CTFS.datafiles}} for structure and definitions of these data
	files and \code{link{CTFS.groupcats}} for how to make these data files.
	Briefly it requires creating a vector of quadrate number for each tree
	using \code{\link{sep.quadinfo}} and then using 
	\code{\link{sep.vct.data}} to create the list.
	
	In order to minimize the time it takes to identify trees that are within
	the defined neighborhood distances, the potential neighbors of a target
	species are first identified by using the quadrate number of the target
	individual and its bordering quadrates.  This is why \code{census.quad}
	is needed to run the functions.  It is the source of the neighboring 
	tree data.  \code{census.spp} provides the
	data for individuals for each target species.
	
	Distance classes are defined as larger than and including the low value 
	to smaller than and not including the high value of the distance classes 
	provided by \code{distbreak}.
	
	If the neighborhood annulus does not lie entirely within the plot then no
	information is accumulated for this target tree at that neighborhood 
	distance (values are set to NA).
	
	Remember that the neighborhoods computed are ANNULI so that the 
	cummulative number of individuals (or basal area) within a given 
	neighborhood distance is computed by summing the annuli values up to 
	that distance.
	
	Remember that the returned values are NOT densities, but counts.  
	Therefore, if densities are needed, then the area of the annuli needs to 
	be taken into account.  However, if the same neighborhood distances are 
	used for all species being investigated, then the counts are comparable 
	since the neighborhood annuli are all of the same sizes for all 
	species.
	
	All trees that are alive or alive but broken are included as valid 
	neighbors for the computation of neighborhood attributes.
	
	For more information on how to change the selection criteria for 
	neighbors or the attributed computed, see 
	\code{\link{CTFS.densdepend}}.
}
\value{
  Returns a \strong{dataframe} of neighborhood attributes for each tree of 
	a single species.  This data frame matches the order and number of trees 
	of a species as found in \code{census.spp}.  The following values are 
	produced for each distance class in \code{distbreak}.  Numbers and basal 
	area are for the annuli as defined by \code{distbreak}.  See 
	\code{\link{socrex90.neigh}} for an example result.
	
  \item{Ncon}{Number of conspecific individuals}
  \item{BAcon}{Total basal area \eqn{(m^2)} of conspecific individuals}
  \item{Nhet}{Number of heterospecific individuals}
  \item{BAhet}{Total basal area \eqn{(m^2)} of heterospecific individuals}
	
	The labels for each variable are created from the dist breaks.  }

\author{Rick Condit and Pamela Hall}
\seealso{
\code{\link{CTFS.densdepend}}, 
\code{\link{CTFS.datafiles}}, \code{\link{CTFS.groupcats}}, 
\code{\link{socrex90.neigh}}
}

\examples{\dontrun{

socrex90.neigh <- neighbordens.1sp("socrex",tst.bci90.quad,tst.bci90.spp)
socrex90.neigh[1:10,]
}}
\keyword{univar}
