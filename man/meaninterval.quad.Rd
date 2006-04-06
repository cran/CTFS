\name{meaninterval.quad}
\alias{meaninterval.quad}
\concept{dynamics}
\concept{census interval}
\concept{quadrate}


\title{Mean Interval Between Census for Quadrates}
\description{
	Computes the mean time between censuses (in years) for any user
	specified size quadrate (area) of plot. The user can also specify any portion of the plot (or all
 of it) be considered.  The default quadrate size is 1 hectare.
}
  }
\usage{
meaninterval.quad(census1, census2, mindbh = 10, 
	plotdim = c(1000, 500), gridsize = 100,  
	xrange = c(0, plotdim[1]), yrange = c(0, plotdim[2]), 
	whichcensus = "both")
}
\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
\item{mindbh}{Minimum DBH for computing population size in each 
  census}
  \item{plotdim}{dimensions of the plot: east-west 1000m and
	north-south 500m}
  \item{gridsize}{side of the square quadrate, 100 x 100m (1 ha) by
	default}
  \item{xrange}{the x axis of the plot used for computation, standard 
  plot has 0-1000}
  \item{yrange}{the y axis of the plot used for computation, standard 
  plot has 0-500}
  \item{whichcensus}{variable to indicate which census(es) are to be
  used.  Values are \dQuote{first}, \dQuote{second}, \dQuote{both}}
}
  
\details{
  This is a generic function for any computation of dynamics.  Any
  tree with a recorded date of when a census was done is used in this 
  computation.  This function is NOT used for the computation of mean 
  time between censuses for the dynamic functions.  
	
	A tree is included in the computation depending upon the value of
	\code{whichcensus}.  It must have a dbh which equals or exceed 
	\code{mindbh} for the census specified.  For \dQuote{both} a tree is 
	included if its dbh exceeds \code{mindbh} for EITHER census.	
	} 
	
\value{ \code{meaninterval.quad} returns a
	\strong{matrix} of the mean interval, in years, for species by user
	defined quadrate.
  
  \item{$meantime}{mean time between intervals for all trees in
  quadrate}
}

\author{Rick Condit and Pamela Hall}

\seealso{}
\examples{\dontrun{
1. Default use of meaninterval.quad
meantime.both.out <- meaninterval.quad(tst.bci90.full, tst.bci95.full)
}}
\keyword{univar}
