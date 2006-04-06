\name{abundance.change}
\alias{abundance.change}
\concept{dynamicc}
\concept{abundance}
\concept{population size}

\title{ Annual Change in Population Size by Categories (User Defined
Groups}
\description{
  Computes annual change in population size for any user defined
  populations (groups).  Any two census datasets can be provided in
  order of date of census.  The rate of change and probability that the 
	rate is = 0 is computed.  For categories of DBH, use \code{\link{abundance.change.dbh}}
	since dbh must be determined independently for each census.
  }
\usage{
abundance.change(census1, census2, alivecode=c("A"),
	mindbh=10, split1=NULL, split2=NULL)
}

\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".}
  \item{mindbh}{Minimum DBH for computing population size in each 
  census}
  \item{split1}{a \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which abundance values are computed.  This vector can be 
  composed of characters or numbers.}
  \item{split2}{a second \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which abundance values are computed.  This vector can be 
  composed of characters or numbers.}
}
\details{
	See \code{\link{CTFS.abundance}} for details on the computation of 
	abundance and associated functions.

	Two censuses must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The vectors \code{split1} and \code{split2} must be of the same
	length as \code{census1} and \code{census2} but can contain
	\code{NA}. 
	
	If dbh class is one of the spilt vectors, use 
	\code{\link{abundance.change.dbh}} instead.  The dbh class of each tree 
	is recalculated for each census in this function.  No rate of change is 
	calculated.
	
	The inclusion of a tree is determined independently for each census.  The
	user can only define \code{status} and \code{mindbh} once.  These values
	are used for assessing a tree for inclusion for each census. 
	}
	
\value{
  \code{abundance.change} returns a \strong{list of arrays} with the values
  of \code{split1} as the first dimension and the values of
  \code{split2} as the second dimension of the array.  The array
  contains the following named components:
  
  \item{$N0}{abundance in first census}
  \item{$N1}{abundance in second census}     
	\item{$S}{number of surviors from first census}
	\item{$R}{number of recruits to second census}
  \item{$rate}{annual rate of change in abundance}
	\item{$prob}{probability that N0 == N1}
  \item{$meanyrs}{mean number of years between census for trees used in 
  mortality rate computation}
  \item{$date0}{mean date of first census}
  \item{$date1}{mean date of second census}
}

\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.abundance}}
}
\examples{\dontrun{
# 1.  Default use of abundance.change
abund.ch.out <- abundance.change(tst.bci90.full,tst.bci95.full)

# 2.  Population change by species
spp.vct <- tst.bci90.full$sp
abund.ch.sp.out <-
abundance.change(tst.bci90.full,tst.bci95.full,split1=spp.vct)

}
}

\keyword{htest}
