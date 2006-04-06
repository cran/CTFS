\name{spp.randarea.range}
\alias{spp.randarea.range}
\concept{random quadrate}
\concept{species area}
\concept{diversity}

\title{Number of Species in Randomly Selected Areas of Various Shapes 
and Sizes}
\description{
}
\usage{
spp.randarea.range(census1, quadside = 100, shape = 1,
	alivecode=c("A"), mindbh = NULL, 
	plotdim = c(1000, 500), replicates = 10, 
	unidennames = c("**", "UNID", "uniden", "UNIDEN"))
}

\arguments{
  \item{census1}{name of census \strong{datafame} for a single census}
  \item{quadside}{vector of dimensions (x) of the randomly selected area
   in meters}
  \item{shape}{vector of integers, factor used to determine dimension y of quadrate.  
  for a square, shape = c(1)}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".  The default for \code{abundance} is "A"}
  \item{mindbh}{Minimum DBH for computing population size in each 
  census}
  \item{plotdim}{plot dimensions, x,y lengths in meters}
  \item{replicates}{number of randomly selected areas to sample}
  \item{unidennames}{vector of codes for unidentified species, codes
  vary among sites and datasets}
}
\details{
	\code{spp.randarea.range} produces a \code{list} and a \code{graph}.
	 The number of species and area of each replicate is graphed.  The
	 graph can be saved using the \strong{File} menu.

	\code{quadside} and \code{shape} must be vectors of equal length,
	for which each value in \code{shape} is paired with a value in
	\code{quadside} of the same location in its respective vector. 
	When these vectors are of unequal length, unexpected results
	including failure to execute, can occur.
	
 	Pay attention to the default values of \code{alivecode} and
	\code{mindbh}.  Any tree with \code{alivecode} not equal to
	\dQuote{A} has a \code{dbh} = NA and so will not contribute to
	the results of this computation.  Setting \code{mindbh} to some
	value other than NULL will affect the trees that are included in
	the calculation.  Use \code{mindbh} to select trees larger than 10
	mm dbh.  Use \code{alivecode} to select trees to be counted as
	being alive but have no \code{dbh} measurement.
}

\value{
  \code{spp.randarea.range} returns a list with 2 components
  \item{spparea }{a \code{dataframe} of the mean and standard
  deviations of the number of individuals and species of the randomly 
  selected areas. \code{spparea} contains:}
  \tabular{ll}{
	xdim \tab x dimension of selected area\cr
  ydim \tab y dimension of selected area\cr
  area \tab area, in hectares, of selected area\cr
  mn.nind \tab mean number of individuals in replicates\cr
  sd.nind \tab standard deviation of number of individuals in 
  replicates\cr
  mn.nspp \tab mean number of species in replicates\cr
  sd.nspp \tab standard deviation of number of species in
  replicates\cr
 }
  \item{full }{a \code{dataframe} of the values for each replicate. 
  \code{full} contains:}
  \tabular{ll}{
  area \tab area, in hectares, of replicate\cr
  n.ind \tab number of individuals in replicate\cr
  n.spp \tab number of species in replicate\cr
	}
}

\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.diversity}}}

\examples{\dontrun{
1.  produces a species area curve in about 2 minutes
quadside=c(1,10,20,30,40,50,60,70,80,90,100)
spp.area.out <- spp.randarea.range(bci90.full,quadside=quadside,replicates=2)
  }
}
\keyword{univar}

