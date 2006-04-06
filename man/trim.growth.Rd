\name{trim.growth}
\alias{trim.growth}
\concept{growth}
\concept{valid growth rate}
\concept{reasonable growth rate}
\concept{maximum growth}
\concept{minimum growth}
\concept{negative growth}

\title{Determines if the Growth Rate of a Tree is Within Defined Limits}
\description{
  Growth rates of trees is computed as the annual change in DBH between 
	censuses (mm / year).  Some of these rates are unreasonably high or low, due to 
	errors in measurement or data entry.  Unreasonable values are identified 
	by this function with user defined maximum and a predetermined minimum.  
	This is a core function as all valid growth rates and DBHs are assessed 
	using the parameters of this function. It is usually called within other 
	functions.  The function returns a logical vector of which trees have valid growth rates.
}
\usage{
trim.growth(dbh1,dbh2,pom1,pom2,status1,status2,timeint,
   slope=0.006214,intercept=.9036,err.limit=4,maxgrow=75)		}
\arguments{
  \item{dbh1}{DBH (mm) of tree at census 1}
  \item{dbh2}{DBH (mm) of tree at census 2}
  \item{pom1}{Point of measurement of tree at census 1}
  \item{pom2}{Point of measurement of tree at census 2}
	\item{status1}{Status of tree at census 1}
	\item{status2}{Status of tree at census 2}
  \item{timeint}{time (years) between censuses}
  \item{slope}{slope of error of measurement line}
  \item{intercept}{intercept of error of measurement line}
  \item{err.limit}{number of standard deviations for allowed for limit of 
	negative growth rate}
  \item{maxgrow}{mm per year for maximum allowable growth rate}
}
\details{
	This function is NOT normally invoked directly by a user.  It is designed 
	to be called from within other functions that compute growth rates.  
	Therefore, the other functions should provide for adjustment of the 
	acceptable limits of growth.  At this time, we recommend that only the maximum limit 
	be adjusted by the user.  The minimum limit is based on empirical data from 
	a single CTFS site (BCI, 1995 and repeated in 2000) and until this same data is collected at other 
	sites, then the BCI data will serve to set the minimum allowable 
	growth rate.
	
  See \code{\link{CTFS.datafiles}} for more information on arguments.  
	
	Only trees that have measured DBHs, a status code of \dQuote{A} in both
	censuses and no change in the point of measurement between censuses, are
	assessed for whether their growth rate is reasonable using this function.
	
	This means that trees that are missing a dbh, have a change in point of 
	measurement, have a dbh < 10 mm in either census or where the growth 
	computes to NA are not assessed for reasonable growth rates.  In these 
	cases, a growth rate cannot be calculated.
	
	The maximum allowable growth rate is just an absolute value in mm per
	year that can (usually) be set by the user, provided that the calling
	function supports it (which it should).
	
	The minimum allowable growth rate (usually negative) is based on a linear
	model estimating the standard deviation of dbh measures due to
	measurement error; the parameters slope and intercept define the linear
	relationship between this error deviation and dbh.  The lower limit of
	growth is where the second dbh measure is more than 4 standard deviations
	below the first.  The default values of slope and intercept are based on
	dbh remeasure tests done in both 1995 and 2000 at BCI. A line was fitted
	through the absolute dbh errors as a function of dbh in both years; the
	average slope and intercept is used here.  
	}
 
 \value{ The function returns a \strong{logical vector} indicating which 
 trees (rows) have valid growth rates.  No estimation of growth rates is 
 done.
}

\author{Rick Condit}
\seealso{\code{\link{CTFS.datafiles}}, \code{\link{CTFS.growth}}, 
\code{\link{CTFS.biomass}}
}
\examples{\dontrun{

dbh1 <- tst.bci90.full$dbh
dbh2 <- tst.bci95.full$dbh
pom1 <- tst.bci90.full$pom
pom2 <- tst.bci95.full$pom
status1 <- tst.bci90.full$status
status2 <- tst.bci95.full$status
time <- (tst.bci95.full$date-tst.bci90.full$date)/365.25

trim.out <- trim.growth(dbh1,dbh2,pom1,pom2,status1,status2,time,
                     slope=0.006214,intercept=.9036,err.limit=4,maxgrow=75)
										 
}}
\keyword{univar}

