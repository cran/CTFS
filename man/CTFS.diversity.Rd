\name{CTFS.diversity}
\alias{CTFS.diversity}
\concept{diversity}
\concept{species area}
\concept{latin name}

\title{Analysis of the number of Species, Genera, and Families in Samples of the Plot}

\description{
  Provides an overview of the functions that count the number of
  species, genera and families for given areas of a plot.  Coordinates
  for specific areas, random samples of areas of equal size or
  different shapes (rectangular only), means for the numbers of
  species, genera or families can be computed with these functions.
  \code{\link{spp.randarea.range}} is the most versitile and general function
  for analysing species area curves and sample number and plot size
  and shape effects.
}

\details{

	LISTS OF FUNCTIONS FOR COMPUTING DIVERSITY AND SPECIES AREA CURVES
	
	\tabular{ll}{
	\code{\link{sgf.area}} \tab Number of Species, Genera and Families for a single defined
Area\cr
	\code{\link{sgf.randarea}} \tab Number of Species, Genera and Families for Randomly Selected
Area\cr
	\code{\link{sgf.meanarea}} \tab Mean Number of Species, Genera and Families for Replicates of Randomly Selected
Areas\cr
	\code{\link{sgf.meanarea.range}} \tab Mean Number of Species, Genera and Families for Replicates of Randomly Selected
Areas of a Given Range of Area Sizes and Shapes\cr
	\code{\link{spp.randarea.range}} \tab Number of Species in Randomly Selected Areas of Various Shapes 
and Sizes\cr
	\code{\link{select.randquad}} \tab Randomly Select the Coordinates of Quadrate
	}
}
\author{Lao, S, Hall, P, Condit, R}
\seealso{\code{CTFS.maps}}
\keyword{univar}
