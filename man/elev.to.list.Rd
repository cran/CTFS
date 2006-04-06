\name{elev.to.list}
\alias{elev.to.list}
\concept{elevation}

\title{Converts a Dataframe of Elevation Values to a Matrix}
\description{
  Converts the dataframe of x, y coordinates and elevation to a list with 
	two objects: the original data frame and the same data in matrix 
	format.  The matrix format is useful for making maps and calculating 
	mean topographic values for quadrates.
}
\usage{
elev.to.list(elevfile, gridsize = 5)
}

\arguments{
  \item{elevfile}{\strong{data frame} of x, y coordinates and elevation, 
	\code{\link{bcielev.info}}}
  \item{gridsize}{scale of elevation values (m)}
}
\details{
This function creates two versions of the same information, as the 
original \strong{data frame} and as a \strong{list}.  The point of this is 
that other functions rely upon one or the other object structure.
  }
\value{
 A \strong{list} composed of
  \item{col}{\strong{data frame} of x, y coordinates and elevation}
		\item{x}{x coordinate, East-West}
		\item{y}{y coordinate, North-South}
		\item{elevation}{elevation, in meters}
  \item{mat}{\strong{matrix} in the orientation of the plot}
}

\author{Rick Condit}
\seealso{\code{\link{CTFS.readwrite}}, \code{\link{map}}, 
\code{\link{meantopo.quad}}}
}
\examples{\dontrun{
bcitopo.info <- elev.to.list(bcielev.info)
str(bcitopo.info)
}}
\keyword{manip}

