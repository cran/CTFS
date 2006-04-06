\name{gxgy.to.hectindex}
\alias{gxgy.to.hectindex}
\concept{tree position}
\concept{quadrate}
\concept{index}
\concept{gx}
\concept{gy}
\concept{hectare}

\title{Converts GX GY Coordinates to a Hectare Number}
\description{
  This function takes an x, y plot location and identifies the hectare
  number.  The far edges of the plot is considered out-of-bounds. The 
  same results can be obtained by using a gridsize of 100 by 100 with 
  the function \code{gxgy.to.index}.
}
\usage{
gxgy.to.hectindex(gx, gy, plotdim = c(1000, 500))
}
\arguments{
	\item{gx}{tree x coordinate}
	\item{gy}{tree y coordinate}
	\item{plotdim}{dimensions of the plot: east-west 1000m and
	north-south 500m}
	}
\details{
 The function should work for any 4 sided plot, square or 
 rectangular. The convention is to number the hectares starting in 
 the southwest corner and going northward up a "column" of 5, and 
 returning to the southern border for 6-10, etc.
}
\value{
  Returns a vector, the hectare number (default values 
  range from 1 to 50).  Returns a -1 when an invalid gx or gy value 
  is given.  
}
\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}}}
	
\examples{\dontrun{
tst.gy=c(0,499,0,499)
tst.gx=c(0,0,999,999)
gxgy.to.hectindex(tst.gx,tst.gy)

# the same values can be obtained using gxgy.to.index with a gridsize 
of 100
gxgy.to.hectindex(783,456)

gxgy.to.index(783,456,gridsize=100)
}
}
\keyword{manip}
