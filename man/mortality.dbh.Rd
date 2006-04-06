\name{mortality.dbh}
\alias{mortality.dbh}
\concept{dynamics}
\concept{mortality}
\concept{dbh}

\title{ Annual Mortality Rates for DBH Classes}
\description{
  Computes annual mortality rate for DBH classes.It is provided as a
  shortcut for using using \code{mortality} with a \code{split1}
  vector of DBH classes.  Any two census datasets can be provided in
  order of date of census.  The annual mortality rate is computed for
  each DBH class.
}

\usage{
mortality.dbh(census1, census2, dbhclass = c(10, 100, 300), 
	alivecode = c("A", "AB", "AS"))
}

\arguments{
   \item{census1}{name of census datafile for first census, must be a
   \strong{dataframe}, must be of same length as \code{census2}}
   \item{census2}{name of census datafile for second census, must be a
   \strong{dataframe}, must be of same length as \code{census1}}
   \item{dbhclass}{vector of lower class boundaries for DBH classes, classes
   are nonoverlapping}
  \item{alivecode}{character, codes of the variable \code{status}
  that indicate the tree is alive.  The most general valid categories
  are: "A" and "AB" and "AS".}
}

\details{
	See \code{\link{CTFS.mortality}} and \code{\link{mortality}}  for details on the computation of 
	mortality rates and associated functions.  
	for 

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The results of \code{mortality.dbh} can be organized into
	\strong{dataframes} with the use of
	\code{\link{assemble.demography}}.
}

	
\value{
  \code{mortality.dbh} returns a \strong{list of arrays} with the 
  following named components:
  
  \item{$N0}{the number of living trees at the first census}
	\item{$S}{the number of surviving trees from the first to second census}
  \item{$D}{the number of trees that died between the first and 
  second census}
  \item{$rate}{the mortality in \%/year}
  \item{$lower}{the lower 95\% confidence interval}
  \item{$upper}{the upper 95\% confidence interval}
  \item{$dbhmean}{the mean dbh in mm for trees used to compute 
  mortality}
  \item{$meanyrs}{mean number of years between census for trees used in 
  mortality rate computation.}
  \item{$date0}{mean date of first census}
  \item{$date1}{mean date of second census}
}

\author{Rick Condit}
\seealso{\code{\link{CTFS.mortality}}}
}

\examples{\dontrun{
1. Default use of mortality
mort.dbh.out <- mortality.dbh(tst.bci90.full)

2.  User defined DBH classes

dbhclass=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300)
mort.dbh.out <- mortality.dbh(tst.bci90.full, classbreak = dbhclass)
}}
\keyword{univar}
