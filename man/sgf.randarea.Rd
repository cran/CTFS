\name{sgf.randarea}
\alias{sgf.randarea}
\concept{species area}
\concept{random quadrate}
\concept{diversity}

\title{Number of Species, Genera and Families for Randomly Selected
Area}

\description{
  Counts the number of species, genera and families of the individuals
  found in a given single defined area of a plot.  There area is
  defined by \code{quaddim} with a randomly selected initial corner.
  This function only identifies a randomly selected point and the
  corners of an area of given size starting from that point. It calls
  \code{sgf.area} to do the counting.
}
\usage{
sgf.randarea(census1, spp.info, quaddim = c(100, 100), 
	alivecode = c("A"), mindbh = NULL, plotdim = c(1000, 500), 
	unidennames = c("**", "UNID", "uniden", "UNIDEN"))
}
\arguments{
  \item{census1}{name of the \strong{datafame} of a single census}
  \item{spp.info}{name of \strong{dataframe} that contains each
  species and information about it, eg. \code{\link{bcispp.info}}}
  \item{quaddim}{dimensions (x,y) of the randomly selected area in
  meters}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".  The default for this function is "A"}
  \item{mindbh}{Minimum DBH for inclusion in computation}
  \item{plotdim}{plot dimensions, x,y lengths in meters}
  \item{unidennames}{vector of codes for unidentified species, codes
  vary among sites and datasets}
}
\details{
  Only one randomly selected area is subsetted and values computed
  with each run of this function.  All areas are squares.  For other
  shapes use \code{\link{spp.randarea.range}}
  
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
  Returns a vector of the following information for the defined area.
  \item{n.ind}{number of individuals}
  \item{n.spp}{number of species}
  \item{n.gen}{number of genera}
  \item{n.fam}{number of families}
  \item{r.sppgen}{ratio of number of species to number of genera}
  \item{r.sppfam}{ratio of number of species to number of families}
  \item{r.genfam}{ratio of number of genera to number of families}
}

\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.diversity}}}

\examples{\dontrun{
1. default use gives the total number of species, genera, families for
a randomly selected area of a given size.  Run multiple times to see
that different values are computed each time.

sgf.randarea(bci90.full,bcispp.info)

}
}
\keyword{univar}
