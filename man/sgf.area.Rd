\name{sgf.area}
\alias{sgf.area}
\concept{species area}
\concept{diversity}

\title{Number of Species, Genera and Families for a single defined
Area}
\description{
  Counts the number of species, genera and families of the individuals
  found in a given single defined area of a plot.  There area is defined by
  x, y coordinates and can be of any size.  This function is used by
  other species area functions.
}
\usage{
sgf.area(census1, spp.info, corners = c(0, 1000, 0, 500), 
	alivecode = c("A"), mindbh = NULL, unidennames = c("**", "UNID", "uniden", "UNIDEN"))
}

\arguments{
  \item{census1}{name of the \strong{datafame} of a single census}
  \item{spp.info}{name of \strong{dataframe} that contains each
  species and information about it, eg. \code{\link{bcispp.info}}}
  \item{corners}{x, y coordinates, in pairs of the corners of the area
 of the plot to be included}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".  The default for this function is "A"}
  \item{mindbh}{Minimum DBH for inclusion in computation}
  \item{unidennames}{vector of codes for unidentified species, codes
  vary among sites and datasets}
}
\details{
	This is the base function for all of the other \code{sgf.*}
	functions.  Given any coordinates, it subsets the dataset for
	trees within those coordinates and calculates the number of
	individuals, species, genera and families.  \code{alivecode},
	\code{mindbh} and \code{unidennames} control which trees are
	included.
	
	Pay attention to the default values of \code{alivecode} and
	\code{mindbh}.  Any tree with \code{alivecode} not equal to
	\dQuote{A} has a \code{dbh} = NA and so will not contribute to
	the results of this computation.  Setting \code{mindbh} to some
	value other than NULL will affect the trees that are included in
	the calculation.  Use \code{mindbh} to select trees larger than 10
	mm dbh.  Use \code{alivecode} to select trees to be counted as
	being alive but have no \code{dbh} measurement.

	All areas are squares.  For other shapes use
	\code{\link{spp.randarea.range}}
	 
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
a given census
sgf.area(bci90.full,bcispp.info)
 
2. for the first hectare
sgf.area(bci90.full,bcispp.info,corners=c=(0,100,0,100))

3. note the differences among values for alivecode and mindbh, keep
track of the default settings for each

sgf.area(bci90.full,bcispp.info,alivecode=c("A","AB","AS"))
sgf.area(bci90.full,bcispp.info,mindbh=NULL)
sgf.area(bci90.full,bcispp.info,alivecode=c("A","AB","AS"),mindbh=NULL)
}
}
\keyword{univar}
