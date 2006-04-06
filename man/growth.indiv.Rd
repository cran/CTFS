\name{growth.indiv}
\alias{growth.indiv}
\concept{dynamics}
\concept{growth}

\title{Annual Growth Rates for Individual Trees (3 methods of growth 
rate computation provided)}

\description{
  Computes annual growths rate for individual trees and returns a
  \strong{dataframe} of these values for all trees.  Two types of growth rates are
  computed: annual dbh increment (mm dbh per year), relative growth
  rate (\% change in dbh per year).  Growth rates can also be
  evaluated for unrealistically high and low values and removed from
  the summary values.  The dataset must have at least 2 censuses to
  compute growth.
}
\usage{
growth.indiv(census1, census2, 
	rounddown = FALSE, mindbh = 10, err.limit = 4, maxgrow = 75)
}
\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}

  \item{rounddown}{logical value in caps.  When TRUE, if either of 
  census is < 55, then the \code{floor} of the dbh value / 5 is 
  provided.  When FALSE, no change in the dbh is made.}
  \item{mindbh}{minimum DBH in \code{census1} for inclusion in computation.}
  \item{err.limit}{number of standard deviations: used for determining 
  if an individual tree growth rate is too high or low for inclusion.}
  \item{maxgrow}{maximum absolute growth rate, mm per year: used for 
  determining if an individual tree growth rate is too high for 
  inclusion.}
}
\details{
	See \code{\link{CTFS.growth}} for details on the computation methods 
	of growth rates and associated functions.

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The results of \code{growth.indiv} are growth rates for individual
	trees, \strong{NOT MEAN GROWTH VALUES}. A value is returned for
	each tree.  If growth is invalid then NA is returned.
}
\value{
	\code{growth.indiv} returns a \strong{dataframe} with the
	following components:
	
	\item{tag}{tree individual tag number}
	\item{sp}{species code}
	\item{gx}{"x" coordinate of tree within the plot}
	\item{gy}{"y" coordinate of tree within the plot}
	\item{dbh1}{DBH at the first census as defined by user}
	\item{dbh2}{DBH at the second census as defined by user}
	\item{years}{interval, in years, between censuses}
	\item{incgr}{annual growth; (dbh2 - dbh1) / time}
	\item{expgr}{relative growth; (log(dbh2)-log(dbh1))/time}
}

	
\author{Rick Condit}
\seealso{\code{\link{CTFS.growth}}
}
\examples{\dontrun{
1. Default use of growth.indiv().  Warning messages may result from
computation of exponential growth.  The values will be set to
NA in the returned \strong{dataframe}

growth.indiv.out <- growth.indiv(tst.bci90.full,tst.bci95.full)

2.  Selecting a single species for growth.indiv()

growth.alsebl.tree < -growth.indiv.out[growth.indiv.out$sp=="alsebl",]

}
}
\keyword{univar}
