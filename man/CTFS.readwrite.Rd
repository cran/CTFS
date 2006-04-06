\name{CTFS.readwrite}
\alias{CTFS.readwrite}
\concept{read data files}
\concept{write data files}

\title{Overview of Functions that Read, Write and Restructure CTFS
data files}

\description{
 Provides an overview of the functions that allow asci data files to
 be read into R using CTFS variable defintions, to be restructured for
 use in various analyses, to be accessed with ease and to be saved for
 further use.
}

\details{
	
	FUNCTIONS THAT READ, WRITE AND RESTRUCTURE CTFS DATA FILES
	
	\item{\code{\link{sep.data}}}{Re-structures a \strong{dataframe} into a
	\strong{list} of \strong{dataframes} by values of a Column in the \strong{dataframe} using a
	variable in the initial \strong{dataframe}.  Replacement for:
	\code{\link{splitdata}}.}

	\item{\code{\link{sep.vct.data}}}{Re-structures a \strong{dataframe}
	into a \strong{list} of \strong{dataframes} by values of a Column in the \strong{dataframe}
	using a vector of a variable not in the initial \strong{dataframe}.
	Replacement for: \code{\link{splitdata}}.}
	
	\item{\code{\link{mergecensus}}}{Merges two Censuses into one
	\strong{dataframe} and sets the value of \code{status}.\cr
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

	FUNCTIONS FOR FORMATTING RESULTS
	
	\item{\code{\link{assemble.demography}}}{Reformat the Output from
	Demographic Functions from \strong{list} to \strong{dataframe}.}
	\item{\code{\link{assemble.table.each}}}{Reformat the Output from
	Demographic Functions from \strong{list} to \strong{dataframe}.}
	}
	
	FUNCTIONS FOR CREATING VECTORS
	
	See \code{\link{CTFS.groupcats}}.
}
\note{
For developers:  The CTFS readwrite functions are found in a file approximately named:
\code{readandwrite.R}
}


\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}},
\code{\link{CTFS.utilities}}, \code{\link{CTFS.groupcats}}}

\keyword{manip}
