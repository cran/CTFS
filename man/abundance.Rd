\name{abundance}
\alias{abundance}
\concept{population size}
\concept{tree density}

\title{Abundance of Trees by Categories (User Defined Groups)}
\description{
  Computes the population size of trees (number alive) for user
  defined groups of trees.  A single dataset is used.  Abundance and other statistics
  are computed for user defined categories.
}

\usage{
abundance(census1, alivecode = c("A"), mindbh = 10,
                    split1 = NULL, split2 = NULL)
}
\arguments{
  \item{census1}{name of census datafile for a single census}
  
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".}
  \item{mindbh}{Minimum DBH for computing population size in each 
  census}
  \item{split1}{a \strong{vector} of categorical values of the same 
  length as \code{datafile} which groups trees into classes of 
  interest for which abundance values are computed.  This vector can be 
  composed of charcters or numbers.}
  \item{split2}{a second \strong{vector} of categorical values of the same 
  length as \code{datafile} which groups trees into classes of 
  interest for which abundance values are computed.  This vector can be 
  composed of charcters or numbers.}
}
\details{
	See \code{\link{CTFS.abundance}} for details on the computation methods
	of abundance and associated functions.  Values for only 1 census can be
	computed at a time.  A tree is included for the computation based on its
	value for \code{status} and \code{dbh} for a single census.

	The vectors \code{split1} and \code{split2} must be of the same
	length as \code{census1} and \code{census2} but can contain
	NAs. 
	
	The results of \code{abundance} can be organized into
	\strong{dataframes} with the use of \code{\link{assemble.demography}}.
}


\value{
  \code{abundance}returns a \strong{list of arrays} with the values
  of \code{split1} as the first dimension and the values of
  \code{split2} as the second dimension of the array.  The array
  contains the following named components:
  
  \item{$N}{the abundance of each category, population size}
  \item{$meandbh}{the mean dbh in mm for trees used to compute 
  abundance}
  \item{$meandate}{mean date of census for trees used to compute abundance}
  
  If the vector(s) \code{split1} and \code{split2} are provided by the
  user, then abundance and associated statistics are computed
  for each value of the vectors.  The vectors are nested so that
  abundance rates is computed for each category of \code{split2}
  within each category of \code{split1}.  Abundance values are
  returned for all levels of each vector and if no value can be
  computed then 0 or \code{NA} is returned as appropriate.

}
\author{Rick Condit, Suzanne Lao and Pamela Hall}
\seealso{\code{\link{CTFS.abundance}}
}
\examples{\dontrun{
1. Default use of abundance
abund.out <- abundance(tst.bci90.full)
names(abund.out)
abund.out$N
abund.out[[1]]

2. Create vectors that categorize each tree by species name 
and habitat.  This is the same process are used for other dynamic 
functions (eg growth, abundance) but the vectors made MUST be put 
together as a list.

spp.vct <- tst.bci9095.full$sp
hab.vct <- sep.quadinfo(tst.bci9095.full,bciquad.info,by.col="hab")->hab.vct
abund.hab <- abundance("tst.bci9095.full",split1=hab.vct)

3. Create vectors that define the quadrate of tree and species to run 
with abundance().  Compare to result of abundance.quad() 

NOTE: the number of dimensions in the returned arrays may NOT be equal to the
number of quadrates in the plot if the census data used did not have
at least 1 tree in each quadrate.  Note the results of these test
files.

spp.vct <- tst.bci90.full$sp

create a vector of quadrates, 2 ways to do it:

quad.vct <- gxgy.to.index(tst.bci90.full$gx,tst.bci90.full$gy,gridsize=20,plotdim=c(1000,500))
quad.vct <- sep.quadinfo(tst.bci90.full,bciquad.info,by.col="quad")

run abundance() and abundance.quad() and compare results

abund.sp.quad.out <- abundance(tst.bci90.full,split1=spp.vct,split2=quad.vct)
abund.quad.out <- abundance.quad(tst.bci90.full,gridsize=20)
identical(abund.sp.quad.out[[1]],abund.quad.out[[1]])
}
}
\keyword{univar}

