\name{select.randquad}
\alias{select.randquad}
\concept{diversity}
\concept{random quadrate}

\title{Randomly Select the Coordinates of Quadrate}

\description{ 
	Takes a vector of quadrate sizes and randomly selects a point on
	the plot as the initial corner.  This function is called by
	\code{spp.randarea.range} and is not usually called alone.
}

\usage{
select.randquad(quadside, shape = 1, plotdim = c(1000, 500), replicates = 
100)
}

\arguments{
   \item{quadside}{vector of dimensions (x) of the randomly selected area
   in meters}
  \item{shape}{single value numeric, factor used to determine dimension y of quadrate.  
  for a square, shape = 1}
  \item{plotdim}{plot dimensions, x,y lengths in meters}
  \item{replicates}{number of randomly selected areas to sample}

}
\details{
  The dimensions of the selected quadrate are determined by the values
  in \code{quadside} and \code{shape}.  The quadrate dimensions are
  x = \code{quadside[i]} and y = \code{quadside[i]}* \code{shape[i]}. 
  If one \code{shape} value is provided then it is used for all
  quadrates.  The returned x and y
  coordinates, paired, define the four corners of the quadrates.
}
\value{
  \item{xlo }{x coordinate}
  \item{ylo }{y coordinate}
  \item{xhi }{x coordinate}
  \item{yhi }{y coordinate}
  \item{area }{area of each quadrate in hectares}
}
\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.diversity}}}

\examples{\dontrun{
1. There are no defaults for this function, parameters must be
supplied. This results in 10 sets of coordinates for 1 hectare, square
quadrates.

select.randquad(quadside=100,shape=1,replicates=10)

2. Try a rectangular quadrate and quadrates of different shapes

select.randquad(quadside=100,shape=2,replicates=10)
select.randquad(quadside=100,shape=c(2,3),replicates=10)
select.randquad(quadside=c(10,100),shape=1,replicates=10)
select.randquad(quadside=c(10,100),shape=c(2,3),replicates=10)
}
}
\keyword{univar}
