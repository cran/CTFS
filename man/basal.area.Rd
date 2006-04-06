\name{basal.area}
\alias{basal.area}
\concept{basal area}
\concept{multiple stems}

\title{ Basal Area by Categories (User defined groups) }
\description{
	Computes basal area (m\^2) for all trees or any user defined
	categorization of trees.  If a tree has more than one stem,
	the other stems recorded in the multi-stem file are included in
	its basal area, if the file exists.  Two datasets are used: one with the primary stem
	and one with all multiple stems for that cenus.  The two datasets
	contain completely non-overlapping sets of DBH. Adding the two
	together produces the total basal area for each tagged tree.
}
\usage{
basal.area(census1, mult1, alivecode = c("A"), mindbh = 10, 
	split1 = NULL, split2 = NULL)
	}

\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}}
  \item{mult1}{name of file (not in quotes) must be a dataframe 
  with trees as rows and column as DBH.  Tree tag numbers can appear
  multiple times in this file if there is more than two stems for a
  tree.}
  \item{alivecode}{character, codes of the variable \code{status}
  that indicate the tree is alive and has a valid DBH. The default is:
  "A" for trees >= 10 mm DBH}
  \item{mindbh}{minimum DBH to determine which trees will be included 
  in computation}
  \item{split1}{a \strong{vector} of categorical values of the same 
  length as \code{datafile} which groups trees into classes of 
  interest for which growth values are computed.  This vector can be 
  composed of charcters or numbers.}
  \item{split2}{a second \strong{vector} of categorical values of the same 
  length as \code{datafile} which groups trees into classes of 
  interest for which growth values are computed.  This vector can be 
  composed of charcters or numbers.}

}
\details{
	See \code{\link{CTFS.basalarea}} for details on the computation methods 
	of basal area and associated functions.

		A tree is included for the computation based on its value for
	\code{status} and \code{dbh} for a single census. 

}

\value{
  \code{basal.area} returns a \strong{list of arrays} with the 
  following named components.  Values for each category of the split 
  vectors are provided.

  \item{$basalarea}{sum of basal area m\^2}
  \item{$meandate}{mean date of enumeration}

	If the vector(s) \code{split1} and \code{split2} are provided by the
	user, then basal area and associated statistics are computed for each
	value of the vectors.  The vectors are nested so that basal area is
	computed for each category of \code{split2} within each category of
	\code{split1}

}
\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.basalarea}}}
\examples{\dontrun{
1. Default use of basal.area
basal.out <- basal.area(tst.bci90.full, tst.bci90.mult)

2. Basal area for trees >= 100 mm (10 cm)
basal100.out <- basal.area(tst.bci90.full, tst.bci90.mult, mindbh=100)

3.  Basal area for each species
spp.vct <- tst.bci90.full$sp
basal.spp.out <- basal.area(tst.bci90.full, tst.bci90.mult,
split1=spp.vct)
}
}
\keyword{univar}
