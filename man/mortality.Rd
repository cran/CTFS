\name{mortality}
\alias{mortality}
\concept{dynamics}
\concept{mortality}

\title{Annual Mortality Rates by Categories (User defined groups)}
\description{
  Computes annual mortality rate for all trees or any user defined
  categorization of trees.  Any two census datasets can be provided in
  order of date of census.  The annual mortality rate and other
  statistics are computed for user defined categories.
  }
\usage{
mortality(census1, census2, alivecode = c("A", "AB", "AS"), 
	split1 = NULL, split2 = NULL)
}
\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
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
	See \code{\link{CTFS.mortality}} for details on the computation of 
	mortality rates and associated functions.

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The vectors \code{split1} and \code{split2} must be of the same
	length as \code{census1} and \code{census2} but can contain NA. 
	
	Take care when creating a split vector based on dbh.  Mortality is
	computed on the surviving trees which have valid dbhs in both censuses.
	However, dead trees only have a valid dbh for the first census of any
	interval.  If the second census is used then no dead trees will be
	identified.  So a dbh vector must be based on the first census of the
	interval.  The inverse if true for recruitment.  Therefore the same dbh
	split vector CANNOT be used for computations of recruitment and
	mortality.  Use a dbh vector based on the census1 for mortality and on 
	census2 for recruitment.
	
	The results of \code{mortality} can be organized into
	\strong{dataframes} with the use of \code{\link{assemble.demography}}.
	}

\value{
  \code{mortality} returns a \strong{list of arrays} with the values
  of \code{split1} as the first dimension and the values of
  \code{split2} as the second dimension of the array.  The array
  contains the following named components:
  
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
  \item{$date1}{mean date of first census}
  \item{$date2}{mean date of second census}
 
  If the vector(s) \code{split1} and \code{split2} are provided by the
  user, then mortality rates and associated statistics are computed
  for each value of the vectors.  The vectors are nested so that
  mortality rates is computed for each category of \code{split2}
  within each category of \code{split1}.  Mortality values are
  returned for all levels of each vector and if no value can be
  computed then 0 or NA is returned as appropriate.
}

\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.mortality}}}

\examples{\dontrun{

1. Default use of mortality()
mort.out <- mortality(tst.bci90.full, tst.bci95.full)
mort.out

2. Create a vector of species names for each tree to compute mortality 
rates for each species.

spp.vct <- tst.bci90.full$sp
mort.spp.out <- mortality(tst.bci90.full, tst.bci95.full, split1 = spp.vct)

3. Create 2 vectors: the first of habitats based on the quad location 
of each tree and the second of tree species names.

habitat.vct <- sep.quadinfo(tst.bci90.full, bciquad.info, by.col = "hab")
spp.vct <- tst.bci90.full$sp

mort.spp.hab.out <- mortality(tst.bci90.full, tst.bci95.full, split1=spp.vct, split2=habitat.vct)
}
}
\keyword{univar}
