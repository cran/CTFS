\name{growth}
\alias{growth}
\concept{dynamics}
\concept{growth}

\title{ Annual Growth Rates by Categories (User defined groups) }
\description{
  Computes annual growth rate for all trees or any user defined
  categorization of trees.  Two growth rates can be computed: simple
  change in dbh over time and relative growth rate.  Unrealistically
  high and low values are removed from the summary values.  Growth
  rate in mm dbh per year, relative growth rate in \% change in dbh
  per year, standard deviation or 95\% confidence limits, and sample
  size are provided.  Two datasets must be used, one for each census.
}
\usage{
growth(census1, census2, 
	rounddown = FALSE, method = "I", stdev = FALSE, 
	mindbh = 10, err.limit = 4, maxgrow = 75, 
	split1 = NULL, split2 = NULL)
}

\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}

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
  \item{mindbh}{minimum DBH in \code{census1} for inclusion in computation.}
  \item{err.limit}{number of standard deviations: used for determining 
  if an individual tree growth rate is too high or low for inclusion.}
  \item{maxgrow}{maximum absolute growth rate, mm per year: used for 
  determining if an individual tree growth rate is too high for 
  inclusion.}
  \item{split1}{a \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which growth values are computed.  This vector can be 
  composed of characters or numbers.}
  \item{split2}{a second \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which growth values are computed.  This vector can be 
  composed of characters or numbers.}
}
\details{ 
	See \code{\link{CTFS.growth}} for details on the computation of 
	growth rates and associated functions.

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The vectors \code{split1} and \code{split2} must be of the same
	length as \code{census1} and \code{census2} but can contain	NA. 
	
	The results of \code{growth} can be organized into
	\strong{dataframes} with the use of \code{\link{assemble.demography}}.
}

\value{
  \code{growth} returns a \strong{list of arrays} with the values
  of \code{split1} as the first dimension and the values of
  \code{split2} as the second dimension of the array.  The array
  contains the following named components:
 
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

  
  If the vector(s) \code{split1} and \code{split2} are provided by the
  user, then growth rates and associated statistics are computed for
  each value of the vectors.  The vectors are nested so that growth
  rates is computed for each category of \code{split2} within each
  category of \code{split1}
}

\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.growth}}
}
\examples{\dontrun{
1. Default use of growth()

growth.all.out <- growth(tst.bci90.full,tst.bci95.full)
growth.all.out

2. Create a vector of species names for each tree to compute growth 
rates for each species.

bci.spp.vct <- tst.bci90.full$sp
names(growth.spp.out)
growth.spp.out$rate

3.  Create 2 vectors: the first of habitats based on the quad location 
of each tree and the second of tree species names.

bci.hab.vct <- sep.quadinfo(tst.bci90.full,bciquad.info,by.col="hab")
growth.spp.hab.out <- growth(tst.bci90.full,tst.bci95.full,split1=bci.spp.vct,split2=bci.hab.vct)
growth.spp.hab.out$rate

}
}
\keyword{univar}

