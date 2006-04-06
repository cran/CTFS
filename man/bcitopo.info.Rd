\name{bcitopo.info}
\alias{bcitopo.info}
\concept{topography}
\concept{elevation}
\docType{data}
\title{Topographic Data for CTFS Plot : Elevation, Convexity and Slope}
\description{
 A \strong{list} composed of two objects.  One is a \strong{data frame} of 
 x, y coordinates and elevation at that point.  The second is the same 
 data in a \strong{matrix} format.  The grid size is 5 x 5 m.
}
\usage{data(bcitopo.info)}
\format{
 A \strong{list} with 2 components:
  \describe{
   \item{col}{\strong{dataframe} of x, y coordinates and elevation}
	    \item{x}{a numeric vector, location on East-West axis, in meters}
	    \item{y}{a numeric vector, location on North-South axis, in 
    meters}
	    \item{elev}{a numeric vector, elevation in m}
	\item{mat}{\strong{matrix} of the \code{elev}}		
	}}
\details{
  These data are used to make topographic maps and to compute mean 
	topographic values for quadrates.  This data set contains the same 
	information as \code{\link{bcielev.info}}.  
}
\source{
 The Center for Tropical Forest Science\cr
\url{http://www.ctfs.si.edu}
}

\examples{\dontrun{
data(bcitopo.info)
names(bcitopo.info)
str(bcitopo.info$col)
maptopo(bcitopo.info$mat,plotdim=c(1000,500))
}}
\keyword{datasets}
