\name{ba.indiv}
\alias{ba.indiv}
\concept{basal area}
\concept{mulitple stems}

\title{Basal Area of Individual Trees, including Mulitiple Stems}
\description{
	 Computes the basal area (m\^2) for each tree including multiple
	 stems if there are any.  If a tree has more than one stem, the
	 other stems recorded in the multi-stem file are included in its
	 basal area.  Two datasets are used: one with the primary stem and
	 one with all multiple stems for that cenus.  The two datasets
	 contain completely non-overlapping sets of DBH. Adding the two
	 together produces the total basal area for each tagged tree.
}
\usage{
ba.indiv(census1, mult1, alivecode=c("A"), mindbh = 10)
}

\arguments{
  \item{census1}{name of census datafile for a single census, a
  \strong{dataframe}}
  \item{mult1}{name of multi-stem file from the same census as}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS". Default for \code{ba.indiv} is "A"}
  \item{mindbh}{minimum DBH to determine which trees will be included 
  in computation}
}


\details{
	See \code{\link{CTFS.basalarea}} for details on the computation methods 
	of basal area and associated functions.
	
	A tree is included for the computation based on its value for
	\code{status} and \code{dbh} for a single census. 
}
\value{
	Returns a vector of basal areas (m\^2) of the same length as
	\code{census1}
}
\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.basalarea}}}
\examples{\dontrun{
1.  Default use of ba.indiv
ba.indiv.out <- ba.indiv(tst.bci90.full,tst.bci90.mult)

2.  Basal area for trees >= 100 mm only
ba.indiv.out <- ba.indiv(tst.bci90.full,tst.bci90.mult,mindbh=100)
}
}
\keyword{univar}