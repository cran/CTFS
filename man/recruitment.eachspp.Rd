\name{recruitment.eachspp}
\alias{recruitment.eachspp}
\concept{dynamics}
\concept{recruitment}
\concept{species}

\title{Annual Recruitment Rate for Species }
\description{
  Computes annual recruitment rate for each species.  Any two census
  datasets can be provided in order of date of census.  The annual
  recruitment rate and other statistics are computed for each species.
  }
  
\usage{
recruitment.eachspp(census1, census2, mindbh=10, 
	alivecode = c("A", "AB","AS"))
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
  }

\details{ 
	See \code{\link{CTFS.recruitment}} for details on the computation of 
	recruitment rates and associated functions.

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The results of \code{recruitment} can be organized into
	\strong{dataframes} with the use of \code{\link{assemble.demography}}.
	}
\value{
  \code{recruit} returns a \strong{list of arrays} with the 
  following named component.  Values are reported for each species. If no value can
	be computed then 0 or NA is returned as appropriate.
  
  \item{$N1}{the number of living trees at the second census}
  \item{$R}{the number of trees that were recruited between the first and 
  second census}
  \item{$rate}{the recruitment rate in \%/year}
  \item{$lower}{the lower 95\% confidence interval}
  \item{$upper}{the upper 95\% confidence interval}
  \item{$time}{mean number of years between census for trees used in 
  recruitment rate computation.}
  \item{$date0}{mean date of first census}
  \item{$date1}{mean date of second census}
}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.recruitment}}}
\examples{\dontrun{
rec.spp.out <- recruitment.eachspp(tst.bci90.full,tst.bci95.full)
rec.spp.out
}
}
\keyword{univar}
