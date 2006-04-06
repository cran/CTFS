\name{splitdata}
\alias{splitdata}
\concept{separate data}
\concept{restructure data file}

\title{
	Re-structures a Dataframe into a List of Dataframes by values of a Column 
	in the Dataframe
}
\description{
	THIS FUNCTION WILL BE DELETED IN FUTURE RELEASES AND REPLACED BY :
	\code{\link{sep.data}} and \code{\link{sep.vct.data}}.  See 
	\code{\link{CTFS.readwrite}} and \code{\link{CTFS.groupcats}} for a 
	description of replacement functions.
	
	Creates a list from a dataframe with each dataframe being the rows 
	from the orignial file that have the same value for the split 
	variable.  There are as many dataframes in the list as there are 
	valid values for the split variable.
	}
\usage{
splitdata(census1, split = census1$sp, allsplit = unique(split))
}
\arguments{
  \item{census1}{any \strong{\dQuote{dataframe}}}
  \item{split}{any column name in the \strong{dataframe}}
  \item{allsplit}{a unique set of values of \code{split}}
}
\details{
	Different columns of the \code{census1} can be defined by using
	the format: \code{census1$columnname}.  In this case DO NOT USE
	QUOTES. Do not ever provide a different parameter value for
	\code{allsplit}.  As seen from the usage, \code{allsplit} is the
	result of a function, not a user defined parameter.  The values of
	the split are sorted alphabetically and missing values, including
	\dQuote{.}  are provided.
	\cr
	\cr This function is very useful for turning any \strong{dataframe} into a 
	\strong{list} of the same records.  For instance, taking the *.full dataset 
	and creating the *.spp dataset.
	\cr
	\cr The function sends a message to the terminal to indicate the 
	progress is it making on each value of the \code{split} variable.
	}
\value{ Returns a \strong{list} of dataframes with tree records from 
one or more census1es}

\author{Rick Condit}
\seealso{\code{\link{CTFS.readwrite}}}
\examples{\dontrun{

# make the bci95.spp from bci95.full
tst.out.1 <- splitdata("tst.bci90.full")

# note this file is identical to the results of sep.data
tst.out.2 <- sep.data("tst.bci90.full")
identical(tst.out.1, tst.out.2)

# specifying a different column in the census
tst.out.codes <- splitdata("tst.bci90.full",split=tst.bci90.full$codes)

}}
\keyword{manip}
