\name{bcielev.info}
\alias{bcielev.info}
\concept{topography}
\concept{elevation}
\docType{data}
\title{Dataset of elevation for 5 x 5 m points in the BCI CTFS Plot}
\description{
	The measured elevation for points in the permanent plot in order
	to map the topography of the site.  The BCI plot was measured at 5
	x 5 m intervals.  This will be replaced with a \strong{list} version of 
	this data.  See \code{\link{bcitopo.info}}
}
\usage{data(bcielev.info)}
\format{
	A data frame of 20301 observation with 3 variables.
	\describe{
		\item{x }{x coordinate, east-west axis}
		\item{y }{y coordinate, north-south axis}
		\item{elev }{elevation at the x, y coordinate}
}
}
\details{
	These data are used to make topographic maps with
	\code{CTFS.maps}.  See \code{\link{bcitopo.info}} for the newer 
	version of this data structure.
}
\source{
The Center for Tropical Forest Science\cr
\url{http://www.ctfs.si.edu}
}


\examples{\dontrun{
data(bcielev.info)
str(bcielev.info)
}
}
\keyword{datasets}

