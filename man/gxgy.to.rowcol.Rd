\name{gxgy.to.rowcol}
\alias{gxgy.to.rowcol}
\concept{tree position}
\concept{quadrate}
\concept{index}
\concept{gx}
\concept{gy}

\title{Converts a GX GY Coordinates to Row and Column Numbers}
\description{
  This functions takes gx gy tree coordinate position and identifies
  the row and column numbers in which it is located.  The row and
  column numbers are based on a gridsize that by default divides the
  plot into 20 by 20 m squares.  The gridsize can be defined by the
  user so other quadrate sizes can be used. The inverse function is 
  \code{rowcol.to.gxgy}.
}
\usage{
gxgy.to.rowcol(gx, gy, gridsize = 20, plotdim = c(1000, 500))
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
 Returns a data frame of row and column numbers.  Returns NA when an 
 invalid gx or gy is given.
}

\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}}}
	
\examples{\dontrun{
tst.gy=c(0,499,0,499)
tst.gx=c(0,0,999,999)
gxgy.to.rowcol(tst.gx,tst.gy)

# the far edges of the plot are out of bounds
tst.gx=c(0,0,1000,1000)
tst.gy=c(0,500,0,500)
gxgy.to.rowcol(tst.gx,tst.gy)
}
}
\keyword{manip}
