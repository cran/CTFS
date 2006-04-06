\name{convert.ba.mult}
\alias{convert.ba.mult}
\concept{basal area}
\concept{mulitiple stems}

\title{Compute Basal Area for All Stems of an Individual Tree using
Census Dataset and Multi-Stem Dataset}
\description{
Computes the total basal area of an individual tagged tree by matching
tree tag numbers of census dataset to multi-stem dataset.  The basal
area of each stem is computed from its dbh and summed together to
provide a single basal area for each tagged tree.  This function is
usually called by other basal area functions, but can be called from
the command line.
}
\usage{
convert.ba.mult(census, mult, mindbh = 10)
}
\arguments{
  \item{census}{name of census datafile for first census, must be a
  \strong{dataframe}}
  \item{mult}{name of file (not in quotes) must be a dataframe 
  with trees as rows and column as DBH.  Must be of the same census
  year as \code{census1}.  Tree tag numbers can be 
  repeated up to the number of valid stems for measuement.}
  \item{mindbh}{Minimum DBH for computing population size in each 
  census.}
}
\details{
	See \code{\link{CTFS.basalarea}} for details on the computation methods 
	of basal area and associated functions.

	When \code{convert.ba.mult} is called from another function, the
	parameter \code{mindbh} should be explicitly passed to assure
	proper behavior.  Its use is determined by the function it is called from.
	
}
\value{
 \code{convert.ba.mult} returns a \strong{dataframe} with the
 following components:
 
 \item{tag}{tag number of tree}
 \item{ba}{sum total basal area (m\^2) of all multiple stems of tree}
 
 Rownames of returned dataframe are the tree tag number.  Not all
 tree tag will be in the dataframe.  Only those with at least 1
 multiple stem.
}
\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.basalarea}}
}
\examples{\dontrun{
1.  Default use of convert.ba.mult
ba.conv.out <- convert.ba.mult(tst.bci90.mult,mindbh=10)

2.  Set mindbh >= 100 mm (10 cm dbh)
ba100.conv.out <- convert.ba.mult(tst.bci90.mult,mindbh=100)


}
}
\keyword{univar}
