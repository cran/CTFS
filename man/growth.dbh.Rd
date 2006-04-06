\name{growth.dbh}
\alias{growth.dbh}
\concept{dynamics}
\concept{growth}
\concept{dbh}

\title{ Annual Growth Rate for DBH Classes }
\description{
  Computes annual growth rate for all trees for DBH classes.  Two
  growth rates can be computed: simple change in dbh over time and
  relative growth rate.  Unrealistically high and low values are
  removed from the summary values.  Growth rate in mm dbh per year,
  relative growth rate in \% change in dbh per year, standard
  deviation or 95\% confidence limits, and sample size are provided.
  Two datasets must be used, one for each census.
}
\usage{
growth.dbh(census1, census2, dbhclass = c(10, 100, 300), 
	rounddown = FALSE, method = "I", stdev = FALSE, 
	err.limit = 4, maxgrow = 75) 
}

\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
  \item{dbhclass}{a vector of lower values for dbh categories}
  \item{rounddown}{logical value in caps.  When TRUE, if either of 
  census is < 55, then the \code{floor} of the dbh value / 5 is 
  provided.  When FALSE, no change in the dbh is made.}
  \item{method}{character indicating the type of growth rate to be 
  calculated.  When "I" the annual change in dbh is provided. When 
  "E" the relative growth rate is provided.  see 
  \code{\link{CTFS.growth}} for more details on \code{method}.}
  \item{stdev}{logical value in caps.  When TRUE, the 
  standard deviation of growth is provided.  When FALSE, a 95\% 
  condfidence interval is provided.  The confidence limits are not 
  computed.}
  \item{err.limit}{number of standard deviations: used for determining 
  if an individual tree growth rate is too high or low for inclusion.}
  \item{maxgrow}{maximum absolute growth rate, mm per year: used for 
  determining if an individual tree growth rate is too high for 
  inclusion.}
  }

  \details{ 
	See \code{\link{CTFS.growth}} for details on the computation of
	growth rates and associated functions.

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The vector for \code{classbreak} can be created outside of the
	function and reference by setting \code{classbreak} equal to that 
	vector.
	
	The results of \code{growth.dbh} can be organized into
	\strong{dataframes} with the use of \code{\link{assemble.demography}}.
}

\value{
  \code{growth.dbh} returns a \strong{list of arrays} with the 
  following named components.  Values are provided for each DBH class.

  \item{$N0}{the number of trees used to compute growth rate}
  \item{$rate}{the annual mean growth rate in mm per year or \% per year,
  type determined by the argument \code{method}}
 
  Only one of the following is provided determined by the argument 
  \code{$sd}{the standard deviation}
  \item{$clim}{the 95\% confidence limit based on a normal 
  distribution}
  
  \item{$dbhmean}{the mean dbh in mm for trees used to compute growth 
  rate}
  \item{$meanyrs}{mean number of years between census for trees used in 
  growth rate computation.}
  \item{$date0}{mean date of first census}
  \item{$date1}{mean date of second census}


}

\author{Rick Condit}
\seealso{\code{\link{CTFS.growth}}
}

\examples{\dontrun{
1. Default use of growth.dbh()

growth.dbh.out <- growth.dbh(tst.bci90.full, tst.bci95.full)
growth.dbh.out$rate

2.  User defined DBH classes
dbh.vct = c(10,20,30,40,50,100,300)
growth.dbh.out <- growth.dbh(tst.bci90.full, tst.bci95.full, classbreak=dbh.vct)
growth.dbh.out$rate
}}
\keyword{univar}
