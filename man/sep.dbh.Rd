\name{sep.dbh}
\alias{sep.dbh}
\concept{DBH class}
\concept{restructure data files}

\title{Creates a Vector of DBH Classes for each Tree in a Dataset}
\description{
	\code{sep.dbh} creates a vector of the same length as the input 
	data file that contains the value for the DBH class of the tree.  
	This vector can be used with \code{\link{mortality}}, 
	\code{\link{growth}}, \code{\link{recruitment}}, and 
	\code{\link{abundance}} to compute rates and counts for each DBH 
	class.  This function can also be used as a 
	template for writing other functions that use the function \code{cut} 
	which "bins" values into user defined bins (classes).
}


\usage{
sep.dbh(census1, dbhclass = c(10, 100, 300))
}

\arguments{
  \item{census1}{ usually a tree data file which is a \strong{dataframe}, 
	but any object with a dim=NULL is also allowed}
  \item{dbhclass}{lower class boundaries for DBH classes, classes 
  are nonoverlapping, no maximum dbh needed}
}

\details{

  This function is usually used to create a vector of dbh classes for each 
	tree in a census file.  However, it can also be used, usually within 
	another function, to generate a vector of dbh classes from a vector of 
	dbh values.  This usage can make such code much more simple as the cut 
	function, labels and class limits are handled by \code{sep.dbh}.
	}
	
\value{
 \code{sep.dbh} returns a \strong{vector} of values which have the 
 mode \code{factor}
}

\author{Pamela Hall}
\seealso{\code{\link{CTFS.utilities}},
\code{\link{CTFS.quadfunctions}}, \code{\link{CTFS.datafiles}}
}


\examples{\dontrun{
1. Default use of sep.dbh
dbh.vct <- sep.dbh(tst.bci90.full)
dbh.vct[1:5]

2. User defined census and dbhclasses
dbh.vct <- sep.dbh(tst.bci90.full,dbhclass=c(10,20,30,40,50))
dbh.vct[1:5]

3. Using a vector of dbh values
dbh1.vct <- sep.dbh(tst.bci90.full$dbh)
}
}
\keyword{manip}
