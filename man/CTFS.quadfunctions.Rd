\name{CTFS.quadfunctions}
\alias{CTFS.quadfunctions}
\concept{quadrates}
\concept{gx}
\concept{gy}
\concept{index}
\concept{location}

\title{Overview of Functions that provide Information on the location 
of trees in a Plot}
\description{
A series of functions that convert from one form of location to 
another, such as from a tree coordinate (gx, gy) to a quadrate, from 
a quadrate to a row and col index, etc.  }   

\details{
	FUNCTIONS FOR CONVERTING AMONG INDEX, COORDINATES AND QUADRATES

	\tabular{ll}{
	\code{\link{rowcol.to.index}} \tab Converts a Row and Column to a Quadrate Number\cr
	\code{\link{index.to.rowcol}} \tab Converts a Quadrate Number (Indext) to Row and Column Numbers\cr
	\code{\link{gxgy.to.index}} \tab Converts GX GY Coordinates to a Quadrate Number\cr
	\code{\link{index.to.gxgy}} \tab Converts a Quadrate Number to GX GY Coordinates\cr
	\code{\link{gxgy.to.rowcol}} \tab Converts a GX GY Coordinates to Row and Column Numbers\cr
	\code{\link{gxgy.to.hectindex}} \tab Converts GX GY Coordinates to a Hectare Number\cr
	\code{\link{findborderquads}} \tab Returns Indices of Neighboring Quadrats Given a Quadrat
	Index\cr
}

 PLOT LAYOUT
 
 The standard plot is 1000 by 500 m, running West to East and South to
 North, like the conventional cooridnates of a graph.
 
 The \bold{0,0} corner and edge going from \bold{0,0 to 0,500} and going from
 \bold{0,0 to 1000,0} are included in the plot. The outer edges, going from
 \bold{500,0 to 1000,500} and \bold{1000,0 to 1000,500} are \bold{NOT} included in the
 plot.  This holds true for any portion of the plot.\cr
 
 Therefore, hectare 1 includes the sides 0,0 to 0,100 and 0,0 to
 100,0 but NOT 100,0 to 100,100 and 100,0 to 100,100.\cr
 Quadrate 1 includes sides 0,0 to 0,20 and 0,0 to 20, 0 but NOT 20,0
 to 20,20 and 0,20 to 20,20.\cr

 The convention is 20 by 20 m quadrates, 1250 of them in
 a 1000 by 500 m plot.  The convention is to number the quadrates
 starting in the southwest corner and going northward up a "column" of
 25, and returning to the southern border for 26-50, etc:
 
 These function should work for any 4 sided plot, square or
 rectangular.  But the dimension of the plot if different from 1000,
 500 must be provided for each function.
 
 Quadrate Number, Index
 \tabular{rrrr}{
         25\tab 50 \tab 75 \ldots \tab 1250\cr
         24\tab 49 \tab 74 \ldots \tab 1249\cr
         . \tab . \tab . \ldots \tab .\cr
         . \tab . \tab . \ldots \tab . \cr
         . \tab . \tab . \ldots \tab . \cr
         3\tab 28 \tab 53 \ldots \tab 1228\cr
         2\tab 27 \tab 52 \ldots \tab 1227\cr
         1\tab 26 \tab 51 \ldots \tab 1226\cr
       }


Rows and Columns
 
 cols\tabular{rrrrrr}{
	1\tab 2\tab 3 \ldots \tab 48 \tab 49 \tab 50\cr}
 rows
 \tabular{rr}{
 	25\tab \ldots \cr
 	24\tab \ldots \cr
 	.\tab \ldots \cr
 	.\tab \ldots \cr	
 	.\tab \ldots \cr
 	2\tab \ldots \cr
 	1\tab \ldots \cr
	}
	
 Hectare Number, Index
 
 \tabular{rrrrrrrrrr}{
         5\tab 10 \tab 15 \tab 20 \tab 25 \tab 30 \tab 35 \tab 40 \tab 45 \tab 50\cr
         4\tab 9 \tab 14 \tab 19 \tab 24 \tab 29 \tab 34 \tab 39 \tab 
         44 \tab 49\cr
         3\tab 8 \tab 13 \tab 18 \tab 23 \tab 28 \tab 33 \tab 38 \tab 
         43 \tab 48\cr
         2\tab 7 \tab 12 \tab 17 \tab 22 \tab 27 \tab 32 \tab 37 \tab 
         42 \tab 47\cr
         1\tab 6 \tab 11 \tab 16 \tab 21 \tab 26 \tab 31 \tab 36 \tab 
         41 \tab 46\cr
       }

 DEFINITIONS
 
 \itemize{
 \item{index}{quadrate number, default 1 to 1250}
 \item{gxgy}{coordinates (m) of a tree relative to plot dimensions, 
 gx determines the column number and gy determines the row number}
 \item{rowcol}{row and column number when viewing plot as a 
 rectangular matrix with quadrate 1 in lower left corner.  This is NOT 
 the usual visual layout of a mathematical matrix, see above}
 \item{hectindex}{hectare number, default 1 to 50}
 \item{plotdim}{dimensions of plot in m, x then y. The default is
 \code{plotdim=c(1000,500)}}
 \item{gridsize}{dimensions of the gridsize, quadrates or hectares in 
 most cases.  The default is \code{gridsize=20} for a quadrate and
 \code{gridsize=100} for a hectare}
 }
 
 \code{gxgy.to.hectindex} has a fixed \code{gridsize=100} but it gives
 the same result as using \code{gxgy.to.index} with \code{gridsize=100}.
  
 EXAMPLES
 
 Examples of how to use these functions, showing the relationship 
 between \strong{index}, \strong{rowcol}, \strong{gxgy}, and 
 \strong{hectindex}.
 
 \preformatted{
 index.to.rowcol(43)
    row col
 1  18   2
 index.to.gxgy(43)
   gx  gy
 1 20 340
 gxgy.to.index(20,340)
 [1] 43
 rowcol.to.index(18,2)
 [1] 43
 gxgy.to.rowcol(20,340)
    row col
 1  18   2
 gxgy.to.hectindex(20,340)
 [1] 4
 
 gxgy.to.index(1000,500)
 [1] NA

 }
}
 
\author{Rick Condit, Suzanne Lao and Pamela Hall}

\seealso{\code{\link{CTFS.groupcats}},
	\code{\link{CTFS.readwrite}}, \code{\link{CTFS.utilities}}
	}

	
\keyword{univar}
