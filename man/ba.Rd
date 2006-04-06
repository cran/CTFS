\name{ba}
\alias{ba}
\alias{basum}
\concept{basal area}

\title{Basal Area of Trees}
\description{
Two functions that provide calculations of basal area in m\^2 from a
dbh vector.  \code{ba} computes the basal area and returns a vector of
the same length as the input vector.  \code{basum} computes the sum of
the basal areas for all values in a vector and returns a single value.
}
\usage{
ba(dbh)
basum(dbh,mindbh=10)
}

\arguments{
  \item{dbh}{DBH in mm of tree}
  \item{mindbh}{minimum dbh for selection of value to be included in
  computation}
}
\details{
 The input is either a single number or a vector of numbers.  Valid 
 values are DBH in mm that are greater than or equal to the minimum dbh.
 If basal area cannot be computed, the value of 0 is returned.
}
\value{
 \code{ba} returns a vector of the same length as the input DBH
 vector, basal area in m\^2
 \code{basum} returns a single value (numeric), sum of the basal area 
 in m\^2
}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.basalarea}}
}
\examples{\dontrun{
ba.out <- ba(tst.bci90.full$dbh)
basum.out <- basum(tst.bci90.full$dbh)

ba100.out <- ba(tst.bci90.full[tst.bci90.full$dbh >= 100,]$dbh)
basum.out <- basum(tst.bci90.full$dbh, mindbh=100)
}}
\keyword{univar}
