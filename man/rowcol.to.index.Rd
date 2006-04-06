\name{rowcol.to.index}
\alias{rowcol.to.index}
\concept{tree position}
\concept{quadrate}
\concept{index}
\concept{gx}
\concept{gy}

\title{Converts a Row and Column to a Quadrate Number}
\description{
  This functions takes row and column numbers and identifies the
  quadrate number (index).  The row and column numbers are based on a
  gridsize that by default divides the plot into 20 by 20 m squares.
  The gridsize can be defined by the user so other quadrate sizes can
  be used.
}
\usage{
rowcol.to.index(rowno, colno, gridsize = 20, plotdim = c(1000, 500))
}

\arguments{
  \item{rowno}{ row number, range from 1 to 25 for default plot}
  \item{colno}{ column number, range from 1 to 50 for default plot}
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
  range from 1 to 1250).  Returns NA when invalid row or column 
  number is provided.
}
 
\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}}}
	
\examples{\dontrun{

rowcol.to.index(21,43)

# out of bounds
rowcol.to.index(26,39)

The 4 corner quadrates

tstrow=c(1,25,1,25)
tstcol=c(1,1,50,50)
rowcol.to.index(tstrow,tstcol)

}
}
\keyword{manip}
