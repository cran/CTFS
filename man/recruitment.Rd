\name{recruitment}
\alias{recruitment}
\concept{dynamics}
\concept{recruitment}

\title{Annual Recruitment Rate by Categories (User defined groups)}
\description{
  Computes annual recruitment rate for all trees or any user defined
  categorization of trees.  Any two census datasets can be provided in
  order of date of census.  The annual recruitment rate and other
  statistics are computed for user defined categories.
  }

\usage{
recruitment(census1, census2, mindbh=10, alivecode = c("A", "AB", "AS"), 
	split1 = NULL, split2 = NULL)
}
\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
  \item{mindbh}{minimum DBH in \code{census1} for inclusion in computation.}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".}
  \item{split1}{a \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which mortality values are computed.  This vector can be 
  composed of characters or numbers.}
  \item{split2}{a second \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which mortality values are computed.  This vector can be 
  composed of characters or numbers.}
  }

\details{ 
	See \code{\link{CTFS.recruitment}} for details on the computation of 
	recruitment rates and associated functions.

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The vectors \code{split1} and \code{split2} must be of the same
	length as \code{census1} and \code{census2} but can contain	NA. 
	
	Take care when creating a split vector based on dbh.  The only valid dbh
	for a recruit is in the second census of any given interval.  Recruits do
	not have a valid dbh for the first census of an interval as they are not
	large enough yet to enter the census.  If the first census is used, then
	no recruits will be identified.  So a dbh vector should be based on the
	dbh of the tree in the second census.  The inverse is true for mortality.
	Therefore the same dbh split vector CANNOT be used for computations of
	recruitment and mortality.  Use a dbh vector based on the census1 for mortality and on 
	census2 for recruitment.
	
	The results of \code{recruitment} can be organized into
	\strong{dataframes} with the use of \code{\link{assemble.demography}}.
	}
\value{
  \code{recruit} returns a \strong{list of arrays} with the values
  of \code{split1} as the first dimension and the values of
  \code{split2} as the second dimension of the array.  The array
  contains the following named components:
  
  \item{$N1}{the number of living trees at the second census}
  \item{$R}{the number of trees that were recruited between the first and 
  second census}
  \item{$rate}{the recruitment rate in \%/year}
  \item{$lower}{the lower 95\% confidence}
  \item{$upper}{the upper 95\% confidence}
  \item{$time}{mean number of years between census for trees used in 
  recruitment rate computation.}
  \item{$date0}{mean date of first census}
  \item{$date1}{mean date of second census}

  If the vector(s) \code{split1} and \code{split2} are provided by the
  user, then recruitment rates and associated statistics are computed
  for each value of the vectors.  The vectors are nested so that
  recruitment rates is computed for each category of \code{split2}
  within each category of \code{split1}.  recruitment values are
  returned for all levels of each vector and if no value can be
  computed then 0 or NA is returned as appropriate.
}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.recruitment}}}

}
\examples{\dontrun{
1. Default use of recruit()
rec.out <- recruitment(tst.bci90.full,tst.bci95.full)
rec.out

2. Create a vector of habitat for each tree based on the quadrate 
location of tree

habitat.vct <- sep.quadinfo(tst.bci90.full,bciquad.info)
rec.hab.out <- recruitment(tst.bci90.full,tst.bci95.full,split1=habitat.vct)

}
}
\keyword{univar}
