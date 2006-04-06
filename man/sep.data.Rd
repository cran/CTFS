\name{sep.data}
\alias{sep.data}
\alias{sep.vct.data}
\concept{restructure files}
\concept{species files}
\concept{neighbors}

\title{	Re-structures a Dataframe into a List of Dataframes by Values in 
a Vector}

	\description{ 
	Creates a list from a dataframe of tree data using a vector
	with values for each tree as the names for the list.  There are as many
	dataframes in the list as there are valid values in the vector.  The
	values in the vector can be derived from a variable in the original
	dataframe or from other information about the trees.  Rows with NA for
	the separation variable can be included or excluded.  
	}

\usage{
sep.data(census1, sepcol = "sp", handle.na = NA)
sep.vct.data(census1, sepvct, handle.na = NA)
}

\arguments{
  \item{census1}{name of census datafile, must be a
  \strong{dataframe}}
  \item{sepcol}{column in datafile with values that will form the
  components of the \strong{list}.  This must be a column in
  \code{census1} and must be the column name.}
  \item{sepvct}{vector of values that will form the components of the 
	\strong{list}.  Must be of same length as \code{census1}.  The vector
	must be created before using this function.}
  \item{handle.na}{instructions on how to handle \code{NA}.  If \code{TRUE},
  missing values in the data are put last; if \code{FALSE}, they are put
  first; if \code{NA}, they are removed.}
}
\details{
  The values of the \code{sepcol} and \code{sepvct} are sorted alphabetically and
  missing values are treated according to \code{handle.na}.

  These functions are very useful for turning any \strong{dataframe}
  into a \strong{list} consisting of the same records and restructured
  so that each component of the list is a value of the variables in
  \code{sepcol} or in \code{sepvct}.
 
  For instance, taking the tst.bci90.full dataset and creating the
  tst.bci90.spp dataset can be done using the defaults of
  \code{sep.data}.
 
  \code{sep.vct.data} can make the same \strong{list} structure
  based on a variable NOT found in \code{census1}.
  \code{sep.vct.data} uses a vector of values which is the length
  as the number of rows in \code{census1}.  This vector can be created
  from any other source of information about the trees eg.  the
  quadrate in which they are located.
	
	The functions \code{\link{sep.dbh}}, \code{\link{sep.quadinfo}} and
	\code{\link{sep.sppinfo}} can be used to make vectors from the dbh of a
	tree, variables in the species information and quadrate information
	files, respectively.
	
  Both functions send a message to the terminal to indicate which
  value of the \code{sepcol} or \code{sepvct} is currently being
  processed.
}
\value{
  Returns a \strong{list} of \strong{dataframes} each of which contain
  tree records.  The values of \code{sepcol} or \code{sepvct} are the
  names of the \strong{list} components.  However, the values of
  \code{sepcol} are returned in the \strong{dataframes} as this
  variable is a part of \code{census1}.  The values for \code{sepvct}
  are NOT returned in the \strong{dataframe} because they are not in
  \code{census1}.
}

\author{Pamela Hall}
\seealso{\code{\link{CTFS.readwrite}}}

\examples{\dontrun{
1. default use
tst.out <- sep.data(tst.bci90.full)
names(tst.out)

2. separate data by status
tst.out <- sep.data(tst.bci90.full,sepcol="status")
names(tst.out)

3. separate data by status including the NAs as the last component in 
the list
tst.out <- sep.date(tst.bci90.full,sepcol="status",handle.na=TRUE)
names(tst.out)

4. separate date by quadrate (which is useful for neighborhood
analysis), leave NAs out of the list

quad.vct <- sep.quadinfo(tst.bci90.full,bciquad.info,by.col="quad")
tst.vct.out <- sep.data.vct(tst.bci90.full,quad.vct)
names(tst.vct.out)
}
}
\keyword{manip}

