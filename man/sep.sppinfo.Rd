\name{sep.sppinfo}
\alias{sep.sppinfo}
\concept{growth form}
\concept{restructure data files}
\concept{species}

\title{ Creates a Vector of Species Characteristics for each Tree in 
a Dataset }
\description{
	\code{sep.sppinfo} creates a vector of the same length as the
	input data file that contains a categorical value that can be
	found in the \code{bcispp.info} data file.  This vector can be
	used with \code{\link{mortality}}, \code{\link{growth}},
	\code{\link{recruitment}}, \code{\link{abundance}} and other functions to provide for
	dynamic rates for growth form, family, etc.  This function can
	also be used as a template for writing other functions that use
	"table-lookup": using a key variable from the data file to look up
	the value of a column on a table file based on the row with the
	key variable value.
	}
\usage{
sep.sppinfo(census1, table.in, by.col = "grform")
}

\arguments{
  \item{census1}{ usually a tree data file which is a \strong{dataframe}, 
	but any object with a dim=NULL is also allowed}
  \item{table.in}{ table file, must be a \strong{data.frame} and is 
  usually the species information file e.g. \code{bcispp.info}}
  \item{by.col}{ the column name of the variable that is to be used 
  for classifying the trees into groups.}
}
\details{
  The \code{table.in} file must contain a variable with the 
  same name as a varible on the \code{\census1} file.  The 
  \code{table.in} file is composed of rows of species with information 
  about each species, including growth form.  The species name of the 
  tree from the \code{census1} is used to locate a row in 
  \code{table.in} and then the value returned is the column with the 
  growth formm.  
  
	This function is usually used to create a vector of species names for each 
	tree in a census file.  However, it can also be used, usually within 
	another function, to generate a vector of species from a vector of 
	species names values.  This usage can make such code much more simple as the cut 
	function, labels and class limits are handled by \code{sep.spp.info}.

  Note that this function can be used to extract any value from 
  \code{bcispp.info} for grouping trees.  More variables can be added 
  to \code{bcispp.info} for grouping trees.
}
\value{
 \code{sep.sppinfo} returns a \strong{vector} of values. The mode of 
 the value depends on what column is extracted from \code{table.in}. 
}
\author{Pamela Hall}
\seealso{\code{\link{CTFS.utilities}},
\code{\link{CTFS.quadfunctions}}, \code{\link{CTFS.datafiles}}
}
\examples{\dontrun{
# 1. Default use of sep.sppinfo
grform.vct <- sep.sppinfo(tst.bci90.full,bcispp.info)

# 2. User defined column "family", then used as a split vector in 
growth()
fam.vct <- sep.sppinfo(tst.bci90.full,bcispp.info,by="fam")
growth.fam.out <- growth("tst.bci90.full",split1=fam.vct)

# 3. Using a vector of species names
grform1.vct <-sep.sppinfo(tst.bci90.full$sp,bcispp.info)

}
}
\keyword{manip}
