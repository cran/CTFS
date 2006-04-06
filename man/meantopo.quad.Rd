\name{meantopo.quad}
\alias{meantopo.quad}
\concept{topography}
\concept{elevation}
\concept{habitat}

\title{Mean Topographic Values for Quadrates}
\description{
  Computes the mean elevation, mean convexity and mean slope for quadrates 
	in CTFS plots.  Edge corrections are made.  These values are computed 
	from elevation measured on 5 x 5 m grid.  These values preserved in the 
	data file \code{\link{bciquad.info}}.
}
\usage{
meantopo.quad(elevfile, gridsize = 20, plotdim = c(1000, 500))
}
\arguments{
  \item{elevfile}{\strong{list} composed of a \strong{data frame} and a 
	\strong{matrix} each which contain elevation data}
  \item{gridsize}{the scale at which the measurement in elevation were
  made.  Must be regular.}
  \item{plotdim}{dimensions of plot}
}
\details{
  This function requires the elevation data to be in a \strong{list}.  
	This data structure can be found in the data set 
	\code{\link{bcitopo.info}} or created using \code{\link{elev.to.list}}.  
	The values returned by this function have been saved in 
	\code{\link{bciquad.info}}
	
	On exit, this function returns the number of columns and rows of the 
	elevation matrix.  In the default case this should be 25 and 30.
}
\value{
  Returns a \strong{data frame} with:
	
	\item{meanelev}{mean elevation (m) for quadrate}
  \item{convex}{convexity for quadrate, position is concave (hill), negative is 
	convex (valley)}
	\item{slope}{mean slope (degrees) for quadrate}	
}

\author{Rick Condit}
\seealso{\code{\link{bcitopo.info}}, 
\code{\link{bciquad.info}},\code{\link{CTFS.readwrite}}, \code{\link{map}}
}
\examples{\dontrun{
data(bcitopo.info)
meanquad.info <- meantopo.quad(bcitopo.info$mat)
}}
\keyword{univar}
