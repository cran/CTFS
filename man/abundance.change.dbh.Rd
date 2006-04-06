\name{abundance.change.dbh}
\alias{abundance.change.dbh}
\concept{dynamics}
\concept{abundance}
\concept{population size}
\concept{dbh}

\title{Annual Change in Population Size by DBH}
\description{
	Computes annual change in population size by DBH classes.  DBH class is
	determined independently for each census.  Any two census datasets can be
	provided in order of date of census.  The rate of change is computed, 
	but no probability is provided.  For categories that do not include DBH
	use \code{\link{abundance.change}}.
	}
\usage{
abundance.change.dbh(census1, census2, alivecode = c("A"), 
   dbhclass = c(10, 100, 300), split = NULL)
}

\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".}
  \item{dbhclass}{a \strong{vector} of dbh class levels - also sets minimum dbh}
  \item{split}{a \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which abundance values are computed.  This vector can be 
  composed of characters or numbers.}
}
\details{
	See \code{\link{CTFS.abundance}} for details on the computation of 
	abundance and associated functions.

	Two censuses must be used.  They do not have to be sequential, only that
	\code{census1} has to be before \code{census2}.
	
	The vector \code{split} must be of the same length as \code{census1} and
	\code{census2} but can contain \code{NA}.  The vector \code{dbhclass} 
	only has to contain the minimum dbh for each class of interest.  A 
	maximum dbh does not have to be specified.  The smallest value of 
	dbhclass is used as \code{mindbh}.
		
	The inclusion of a tree is determined independently for each census.  The
	user can only define \code{status} and \code{dbhclass} once.  These values
	are used for assessing a tree for inclusion for each census. 
}

\value{
  \code{abundance.change.dbh} returns a \strong{list of arrays} with the values
  of \code{split1} as the first dimension and the values of
  \code{split2} as the second dimension of the array.  The array
  contains the following named components:
  
  \item{$N0}{abundance in first census}
  \item{$N1}{abundance in second census}     
  \item{$rate}{annual rate of change in abundance}
  \item{$meanyrs}{mean number of years between census for trees used in 
  mortality rate computation}
  \item{$date0}{mean date of first census}
  \item{$date1}{mean date of second census}
}

\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.abundance}}
}
\examples{\dontrun{
# 1. Default use of abundance.change.dbh
abund.ch.dbh.out <- abundance.change.dbh(tst.bci90.full,tst.bci95.full)

# 2. Set different dbhclasses and also by species
abund.ch.dbh.sp.out <- 
abundance.change.dbh(tst.bci90.full,tst.bci95.full, +
  dbhclass=c(10,20,30,100,300),split=tst.bci90.full$sp)
}}
\keyword{univar}

