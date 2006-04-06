\name{gxgy.to.index}
\alias{gxgy.to.index}
\concept{tree position}
\concept{quadrate}
\concept{index}
\concept{gx}
\concept{gy}

\title{Converts GX GY Coordinates to a Quadrate Number}
\description{
  This function takes an x, y plot location and identifies 
  the quadrate number.  The quadrate is based on a gridsize that by 
  default divides the plot into 20 by 20 m squares.  The gridsize can be 
  defined by the user so other quadrate sizes can be used. The inverse 
  function is \code{index.to.gxgy}.
}
\usage{
gxgy.to.index(gx, gy, gridsize = 20, plotdim = c(1000, 500))
}

\arguments{
	\item{gx}{tree x coordinate}
	\item{gy}{tree y coordinate}
	\item{gridsize}{side of the square quadrate, 20 x 20m by default}
	\item{plotdim}{dimensions of the plot: east-west 1000m and
	north-south 500m}
	}
  
\details{
 The function should work for any 4 sided plot, square or 
 rectangular.  The convention is 20 by 20 m quadrates, 1250 of them in a 
 1000 by 500 m plot.  The convention is to number the quadrates starting in 
 the southwest corner and going northward up a "column" of 25, and 
 returning to the southern border for 26-50, etc.  The columns are 
 numbered 1 to 50, west to east (left to right).  The rows are 
 numbered 1 to 25, south to north (bottom to top).
}
\value{
  Returns a vector, the index number of the quadrate (default values 
  range from 1 to 1250).  Returns NA when an invalid gx or gy is 
  provided.
}
\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}}}
	
\examples{\dontrun{

# the corners:
tst.gx=c(0,0,999,999)
tst.gy=c(0,499,0,499)
gxgy.to.index(tst.gx,tst.gy)

# the far edge from 0,0 is out of bounds
gxgy.to.index(500,1000)

# any location in the conventional plot
gxgy.to.index(321,123)

# same location in smaller plot yields same index
gxgy.to.index(321,123,plotdim=c(500,500))

# index for hectare sized quadrates
gxgy.to.index(321,123,gridsize=100)
}
}
\keyword{manip}

