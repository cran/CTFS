\name{index.to.gxgy}
\alias{index.to.gxgy}
\concept{tree position}
\concept{quadrate}
\concept{index}
\concept{gx}
\concept{gy}
\concept{location}

\title{Converts a Quadrate Number to GX GY Coordinates}
\description{
  This function takes a quadrate (index) number and returns the x, y plot 
  location of the quadrate corner.  The corner is the southwest or 
  lower left of the quadrate.  The quadrate is based on a gridsize that by 
  default divides the plot into 20 by 20 m squares.  The gridsize can be 
  defined by the user so other quadrate sizes can be used. The far 
  edges of the plot are considered out-of-bounds. The inverse 
  function is \code{gxgy.to.index}.
}
\usage{
index.to.gxgy(index, gridsize = 20, plotdim = c(1000, 500))
}
\arguments{
	\item{index}{quadrate, valid default values 1 to 1250}
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
  Returns a data frame which contains the x and y coordinate of the 
  southwest (or lower left) corner of the quadrate provided.  Returns 
  NA when an invalid quadrate is given.  
}
\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}}}
	
\examples{\dontrun{

tst.corners=c(1,25,1226,1250)
index.to.gxgy(tst.corners)

# invalid index
index.to.gxgy(1251)
}
}
\keyword{manip}
