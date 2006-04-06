\name{ba.change}
\alias{ba.change}
\concept{basal area}
\concept{dynamics}
\concept{multiple stems}

\title{Annual Change in Basal Area by Categories (User Defined Group)}
\description{
  Computes annual change in basal area (m\^2 per year) for any user
  defined groups.  If a tree has more than one stem, the other stems
  recorded in the multi-stem file are included in its basal area, if
  the files exist.  Two datasets are used: one with the primary stem
  and one with all multiple stems for that cenus.  The two datasets
  contain completely non-overlapping sets of DBH. Adding the two
  together produces the total basal area for each tagged tree.  Two
  censuses are required to compute annual change.
}
\usage{
ba.change(census1, census2, mult1 = NULL, mult2 = NULL, 
	alivecode = c("A"), mindbh = 10, 
	split1 = NULL, split2 = NULL)
}
\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
  \item{mult1}{name of file (not in quotes) must be a dataframe 
  with trees as rows and column as DBH.  Tree tag numbers can appear
  multiple times in this file if there is more than two stems for a
  tree.}
  \item{mult2}{name of file (not in quotes) must be a dataframe 
  with trees as rows and column as DBH.  Tree tag numbers can appear
  multiple times in this file if there is more than two stems for a
  tree.}
  \item{alivecode}{character, codes of the variable \code{status} 
  that indicate the tree is alive.  The most general valid categories are: "A" and 
  "AB" and "AS".}
  \item{mindbh}{Minimum DBH for computing population size in each 
  census}
  \item{split1}{a \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which basal area values are computed.  This vector can be 
  composed of characters or numbers.}
  \item{split2}{a second \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which basal area values are computed.  This vector can be 
  composed of characters or numbers.}
}
\details{
	See \code{\link{CTFS.basalarea}} for details on the computation methods 
	of basal area and associated functions.

	The inclusion of a tree is determined independently for each census.  The
	user can only define \code{status} and \code{mindbh} once.  These values
	are used for assessing a tree for inclusion for each census. 

}

\value{
  \code{ba.change} returns a \strong{list of arrays} with the values
  of \code{split1} as the first dimension and the values of
  \code{split2} as the second dimension of the array.  The array
  contains the following named components:
  
  \item{$B0}{abundance of species in census 1}
  \item{$B1}{abundance of species in census 2}
  \item{$rate}{annual rate of change in abundance}
  \item{$time}{mean number of years between census for trees used in 
  change in basal area computation.}
  \item{$date1}{mean date of first census}
  \item{$date2}{mean date of second census}
  }
\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.basalarea}}}
\examples{\dontrun{
1.  Default use of ba.change
ba.change.out <- ba.change(tst.bci90.full, tst.bci95.full, tst.bci90.mult,
tst.bci95.mult)

2.  Change in basal area for species
spp.vct <- tst.bci90.full$sp
ba.change.out <- ba.change(tst.bci90.full, tst.bci95.full, tst.bci90.mult,
tst.bci95.mult, split1=spp.vct)

3.  Basal area dynamics for larger trees
ba.change.large.trees.out <- ba.change(tst.bci90.full, tst.bci95.full,
tst.bci90.mult, tst.bci95.mult, mindbh=100)
  }
}
\keyword{univar}