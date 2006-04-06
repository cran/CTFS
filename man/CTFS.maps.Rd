\name{CTFS.maps}
\alias{CTFS.maps}
\concept{map}
\concept{topography}
\concept{elevation}
\concept{habitat}
\concept{point pattern}

\title{Making Maps of Tree Positions, Topography and Area Characteristics}

\description{
	Functions that draw (and save) maps of tree location, topography,
	and area characteristics.  Trees locations can be represented by
	symbol sized relative to DBH. Area characteristics include habitat
	types.  Any characteristics that can be expressed as a coded
	variable for a given area or sized quadrate can be mapped.  All
	three of these types of maps (individual tree, topography and area
	characteristics) can be superimposed.  Maps can be saved in a
	variety of formats and generated for any subset of the dataset.
	These functions can be used to create the maps for each species in
	the site stand tables.
	}

\details{
	
	LIST OF FUNCTIONS FOR DRAWING MAPS
	
	\tabular{ll}{
	\code{\link{map}} \tab Draws a map of the location of each tree of 
	given species, more than one can be plotted.  Elevation contours can be included \cr
	\code{\link{maptopo}} \tab Draws a map of elevation contours \cr
	\code{\link{map1species}} \tab Draws a map of a single species.  This is 
	called by \code{\link{map}}. \cr
	}
	
  WHAT MAPS CAN BE DRAWN
	
	Elevation data has been either measured or kriged to 5 
	x 5 m at each site.  A contour map can be drawn at this scale or 
	larger.  The elevation data must be converted to a matrix and only the 
	matrix is used in \code{\link{map}}.  See \code{\link{bcitopo.info}} for 
	further information on how to make the matrix from the measured 
	elevation data.
	
	Any number of species can be plotted.  The default list of colours works 
	best on 3 or 4 species.  If plotting more than that then a colour scheme 
	such as \code{rainbow} or \code{terrain.colors} should be used.  Shades 
	of gray can also be used.
	
	The \code{census.spp}	argument can be constructed of any interesting 
	group of trees and/or species as long as the file is a list and 
	\code{spname} are the names of the tree census data in the list.  Each 
	name is taken as a color (or shade of gray) for plotting.  
	
}
\author{Pamela Hall}
\seealso{\code{\link{CTFS.diversity}}}
\keyword{hplot}

