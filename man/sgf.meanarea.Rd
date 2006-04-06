\name{sgf.meanarea}
\alias{sgf.meanarea}
\concept{species area}
\concept{random quadrate}
\concept{diversity}

\title{Mean Number of Species, Genera and Families for Replicates of Randomly Selected
Areas}
\description{
  Computes the mean and standard deviations of the number of species,
  genera and families of the individuals for a given number of
  relicates of randomly selected areas.  There area is defined by
  \code{quaddim} with a randomly selected initial corner.  It calls
  \code{sgf.area} to do the counting.
}
\usage{
sgf.meanarea(census1, spp.info, quaddim = c(100, 100), 
	alivecode = c("A"), mindbh = NULL, plotdim = c(1000, 500), 
	replicates = 50, unidennames = c("**", "UNID", "uniden", "UNIDEN"))
}

\arguments{
  \item{census1}{name of the \strong{datafame} of a single census}
  \item{spp.info}{name of \strong{dataframe} that contains each
  species and information about it, eg. \code{\link{bcispp.info}}}
  \item{quaddim}{dimensions (x,y) of the randomly selected area, in
  meters}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".  The default for this function is "A"}
  \item{mindbh}{Minimum DBH for inclusion in computation}
  \item{plotdim}{plot dimensions, x,y lengths, in meters}
  \item{replicates}{number of randomly selected areas to sample}
  \item{unidennames}{vector of codes for unidentified species, codes
  vary among sites and datasets}
}
\details{
  The number of samples is determined by \code{replicates}. Areas must
  be rectangles but do not need to be squares.  Only one area shape
  and size is used in each function call.  For a range of area types
  and sizes use \code{\link{sgf.meanarea.range}}.

  Pay attention to the default values of \code{alivecode} and
  \code{mindbh}.  Any tree with \code{alivecode} not equal to
  \dQuote{A} has a \code{dbh} = NA and so will not contribute to the
  results of this computation.  Setting \code{mindbh} to some value
  other than NULL will affect the trees that are included in the
  calculation.  Use \code{mindbh} to select trees larger than 10 mm
  dbh.  Use \code{alivecode} to select trees to be counted as being
  alive but have no \code{dbh} measurement.
}

\value{
  Returns a vector of the means and standard deviations of the following information.
  \item{mn.nind}{mean number of individuals}
  \item{sd.nind}{standard deviation of the number of individuals}
  \item{mn.nspp}{mean number of species}
  \item{sd.nspp}{standard deviation of the number of species}
  \item{mn.ngen}{mean number of genera}
  \item{sd.ngen}{standard deviation of the number of genera}
  \item{mn.nfam}{mean number of families}
  \item{sd.nfam}{standard deviation of the number of families}
  \item{mn.rsppgen}{mean of the ratio of number of species to number of genera}
  \item{sd.rsppgen}{standard deviation of the ratio of the number of
  species to number of genera}
  \item{mn.rsppfam}{mean of the ratio of number of species to number of 
  family}
  \item{sd.rsppfam}{standard deviation of the ratio of the number of
  species to number of family}
  \item{mn.rgenspp}{mean of the ratio of number of genera to number of 
  families}
  \item{sd.rgenspp}{standard deviation of the ratio of the number of
  genera to number of families}
}

\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.diversity}}}

\examples{\dontrun{
1. almost default use, with fewer replicates to get a result in a
shorter amount of time

sgf.out <- sgf.meanarea(bci90.full,bcispp.info,replicates=10)

}
}
\keyword{univar}

