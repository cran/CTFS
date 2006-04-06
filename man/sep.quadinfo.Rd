\name{sep.quadinfo}
\alias{sep.quadinfo}
\concept{habitat}
\concept{restructure data files}

\title{ Creates a Vector of Quadrate Characteristics for each Tree in 
a Dataset }
\description{
	\code{sep.quadinfo} creates a vector of the same length as the
	input data file that contains a categorical value that can be
	found in the \code{bciquad.info} data file.  This vector can be
	used with \code{\link{mortality}}, \code{\link{growth}},
	\code{\link{recruitment}}, \code{\link{abundance}} and other functions to provide for
	dynamic rates for habitat, quad, etc. This function can also be used as a 
	template for writing other functions that use "table-lookup": using a 
	key variable from the data file to look up the value of a column on a 
	table file based on the row with the key variable value.  This function 
	is also used to create the data file structured needed for density 
	dependence analysis, see \code{\link{CTFS.densdepend}} and 
	\code{\link{neighbordens.1sp}}
	}
\usage{
sep.quadinfo(census1, table.in, by.col = "hab", 
   gridsize = 20, plotdim = c(1000,500))
}
\arguments{
  \item{census1}{ tree data file, must be a \strong{data.frame}}
  \item{table.in}{ table file, must be a \strong{data.frame} and is 
  usually the species information file e.g. \code{bciquad.info}}
  \item{by.col}{ the column name of the variable that is to be used 
  for classifying the trees into groups.}
	\item{gridsize}{dimension of square quadrate, side length in meters}
	\item{plotdim}{plot dimensions, x,y lengths in meters}

}
\details{
  The \code{table.in} file must contain a variable with the
  same name as a variable on the \code{data.in} file.  The
  \code{table.in} file is composed of rows of quadrates with
  information about each quadrate including its habitat category under
  some given scheme of classification.  The gx and gy value of the
  tree is used to locate the quadrate number (using
  \code{\link{gxgy.to.index}}.  Then the quadrate number is used is
  used to locate a row in \code{table.in} and then the value returned
  is the column with the habitat value for that quadrate.
  
  Note that this function can be used to extract any value from 
  \code{bciquad.info} for grouping trees.  More variables can be added 
  to \code{bciquad.info} for grouping trees.
	
	Use this function with \code{by.col} = \dQuote{quad} to create a vector 
	containing the quadrate number of each tree.  This can then be used with 
	\code{\link{sep.vct.data}} to create the \strong{list} of 
	\strong{dataframes} of quadrates containing tree data.
}
\value{
 \code{sep.quadinfo} returns a \strong{vector} of values. The mode of 
 the value depends on what column is extracted from \code{table.in}. 
 }
\author{Pamela Hall}
\seealso{\code{\link{CTFS.utilities}},
\code{\link{CTFS.quadfunctions}}, \code{\link{CTFS.datafiles}}
}

\examples{\dontrun{
1. Default use of sep.quadinfo()
hab.vct <- sep.quadinfo(tst.bci90.full,bciquad.info)
length(hab.vct)

2. Using sep.quadinfo for calculation of mortality and recruitment by 
habitat
hab.vct <- sep.quadinfo(tst.bci90.full,bciquad.info)
mort.hab.out <- mortality(bci90.full,split1=hab.vct)
rec.hab.out <- recruit(bci90.full,split1=hab.vct)
}
}
\keyword{manip}
