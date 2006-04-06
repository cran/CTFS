\name{CTFS.utilities}
\alias{CTFS.utilities}
\concept{utilities}
\concept{dates}
\concept{quadrates}
\concept{coordinates}

\title{Overview of Functions that are CTFS Utilities}
\description{
 Provides an overview of a variety of functions that perform tasks
 associated with other CTFS functions.  Only those that are
 appropriate to use alone are listed here.  Utilities that are only
 called by other functions are found in \code{CTFS-internal}.
}
\details{

CTFS UTILITIES

1.  Functions for Miscellaneous Computations

	\item{\code{\link{find.climits}}}{Computes 95\% CI for a Single Annual Mortality
	Rate.}
	
	\item{\code{\link{fill.dimension}}}{Fills all the dimensions of a 2 dimensional array.}
	
	\item{\code{\link{fill.1dimension}}}{Fills all the dimensions of a 1 dimensional array.}
	
	\item{\code{\link{rndown5}}}{Rounds down DBH where measurement was less 
	precise than 1 mm.}
	
	\item{\code{\link{trim.growth}}}{Removes extreme values of growth using user 
	defined limits.}
	
	\item{\code{\link{unwind.matrix}}}{Turns a matrix into a vector,
	\dQuote{unwinding} it by columns.}
	
	\item{\code{\link{xydist}}}{Distance between x,y coordinates.}

	\item{\code{\link{xydistvect}}}{Distance between vector of x,y
	coordinates, x as first element, y as second.  Used with
	\code{apply}.}
	
2.  Functions for Converting Dates
	\item{\code{\link{tojulian}}}{Converts Calendar Dates to Julian Dates and Back.}
	\item{\code{\link{fromjulian}}}{Converts Calendar Dates to Julian Dates and Back.}
	
3.  Dunctions for Creating Grouping Variables for Trees 
	see \code{link{CTFS.groupcats}}

	\item{\code{\link{sep.quadinfo}}}{Creates a Vector of Quadrate
	Characteristics for each Tree in a Dataset.}
	
	\item{\code{\link{sep.sppinfo}}}{Creates a Vector of Species
	Characteristics for each Tree in a Dataset.}
	
	\item{\code{\link{sep.dbh}}}{Creates a Vector of DBH Classes for
	each Tree in a Dataset.}

4.  Functions for Manipulating Quadrates and X, Y Coordinates
	see \code{\link{CTFS.quadfunctions}}\cr
	
	\item{\code{\link{rowcol.to.index}}}{Converts a Row and Column to a Quadrate Number}
	\item{\code{\link{index.to.rowcol}}}{Converts a Quadrate Number (Indext) to Row and Column Numbers}
	\item{\code{\link{gxgy.to.index}}}{Converts GX GY Coordinates to a Quadrate Number}
	\item{\code{\link{index.to.gxgy}}}{Converts a Quadrate Number to GX GY Coordinates}
	\item{\code{\link{gxgy.to.rowcol}}}{Converts a GX GY Coordinates to Row and Column Numbers}
	\item{\code{\link{gxgy.to.hectindex}}}{Converts GX GY Coordinates to a Hectare Number}
	
5.  Functions for Reading, Writing and Converting CTFS Data Files
	see \code{\link{CTFS.readwrite}}\cr
	
	\item{\code{\link{sep.data}}}{Re-structures a \strong{dataframe} into a
	\strong{list} of \strong{dataframes} by values of a Column in the \strong{dataframe} using a
	variable in the initial \strong{dataframe}.  Replacement for:
	\code{\link{splitdata}}.}

	\item{\code{\link{sep.vct.data}}}{Re-structures a \strong{dataframe}
	into a \strong{list} of \strong{dataframes} by values of a Column in the \strong{dataframe}
	using a vector of a variable not in the initial \strong{dataframe}.
	Replacement for: \code{\link{splitdata}}.}
	
	\item{\code{\link{mergecensus}}}{Merges two Censuses into one
	\strong{dataframe} and sets the value of \code{status}.}
	\item{\code{\link{create.status}}}{Sets the value of \code{status}
	for a single census dataset.}
	
	\item{\code{\link{texttordata}}}{Creates CTFS R \strong{dataframe} from any tab delimited text file.}
	
	\item{\code{\link{sep.species}}}{Separate One or More Dataframes of Tree Records from a List.
  If a single species name is provided, a data frame of the tree records 
	of that species is returned.  This is the same as \code{load.species}.  
	If a vector of species names is provided, then a \strong{list} of the 
	\strong{dataframes} of those species is returned.}

	\item{\code{\link{load.species}}}{Separate One or More Dataframes of Tree Records from a List.
	A more limited function than \code{sep.species}.  The same effect can be
	achieved with \code{siteyear.spp$sppname} where \code{siteyear.spp} is a
	\strong{list} created using \code{sep.data}, eg.  \code{tst.bci90.spp}
	and \code{sppname} is one of the \strong{list} names.}

	\item{\code{\link{readelevdata}}}{Read a Text or RData \strong{dataframe} and
	Make a \strong{matrix}, Specifically for Elevation Data.}
	

6.  Functions for Organizing Results from other Runctions
	see \code{\link{CTFS.readwrite}}

	\item{\code{\link{assemble.demography}}}{Reformat the Output from
	Demographic Functions from \strong{list} to \strong{dataframe}.}
	\item{\code{\link{assemble.table.each}}}{Reformat the Output from
	Demographic Functions from \strong{list} to \strong{dataframe}.}

	}
	
\note{
For developers:  The date conversion, quad and  functions are found in
\dQuote{utilities.R}.  The read and write and creating groups functions are found in
\dQuote{readandwrite.R}. 
}
\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}},
\code{\link{CTFS.readwrite}}, \code{\link{CTFS.groupcats}}
}
\keyword{manip}
