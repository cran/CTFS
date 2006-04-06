\name{readelevdata}
\alias{readelevdata}
\concept{elevation}

\title{Read a Text or RData Dataframe and Make a Matrix, Specifically 
for Elevation Data}
\description{
	Read a tab delimited text file or an RData \strong{dataframe}
	containing the x, y coordinates and elevation at those coordinates
	 The standard file has measurements at every 5 by 5 m point in the
	 plot.  A matrix is returned for use in topographic mapping
	 functions.
}
\usage{
readelevdata(elevfile, gridsize = 5, text = FALSE)
}

\arguments{
  \item{elevfile}{file containing elevation data.  \code{bcielev.info}
 is the RData \strong{dataframe}.  Or specify a text, tab delimited
 file.}
  \item{gridsize}{the scale at which the measurement in elevation were
  made.  Must be regular.}
  \item{text}{indicates where the \code{elevfile} is a \strong{text}
  file or a \strong{dataframe}. }
}
\details{
    This function is used to convert a spreadsheet text file or a
	\strong{dataframe} to a matrix so that topographic maps can be
	made of the plots.
}
\value{
	Returns a \strong{matrix} with rows and columns in the same
	pattern of the actual measurement in the plot.  In the case of the
	standard plot, 1000 by 500 m, measured at 5 m intervals, there are
	101 rows and 201 columns.  All edges of the plot are included.\cr  
	
	1000 / 5 = 200 + 1 (outer edge) = 201 columns\cr
	500 / 5 = 100 + 1 (outer edge) = 101 rows\cr

}
\author{Rick Condit and Suzanne Lao}
\seealso{\code{\link{CTFS.readwrite}}}

\examples{\dontrun{
readelevdata(bcielev.info)->elev.mat
dim(elev.mat)
bcielev.info[1:10,]
elev.mat[1:10,1:10]
  }
}
\keyword{manip}