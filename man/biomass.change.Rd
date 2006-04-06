\name{biomass.change}
\alias{biomass.change}
\alias{biomass.change.dbh}
\concept{biomass}
\concept{dynamics}

\title{Biomass Dynamics}
\description{
 Summarizes the dynamics of above ground biomass (AGB) by computing Mg of
 biomass in each census including the amount loss due to mortality, amount
 gained due to recruitment and growth of persistent trees.  Second function
 is used for dividing these quantities into dbh classes.
}
\usage{
biomass.change(biomass.out, split1 = NULL, split2 = NULL)
biomass.change.dbh(census1, census2, biomass.out, 
   dbhclass = c(10, 100, 300), split = NULL)
}

\arguments{
  \item{biomass.out}{result from \code{\link{biomass}} which is a 
	\strong{dataframe} of dimensions \code{census}}
  \item{split1}{a \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which growth values are computed.  This vector can be 
  composed of characters or numbers.}
  \item{split2}{a second \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which growth values are computed.  This vector can be 
  composed of characters or numbers.}
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
  \item{dbhclass}{a vector of lower values for dbh categories}
  \item{split}{a second \strong{vector} of categorical values of the same 
  length as \code{census} which groups trees into classes of 
  interest for which growth values are computed.  This vector can be 
  composed of characters or numbers. This is the second split value for 
	\code{biomass.change.dbh}}
}
\details{
The output from biomass is a dataframe of the same dimensions as census1 or
census2 with AGB for each tree where a DBH value was available.  The
presence or absence of a DBH value in a census indicates whether the tree
died during the census interval (AGB census 2 is missing), was recruited
(AGB census 1 is missing) or persisted in both censuss (AGB not missing).
AGB may be missing in both censuses is the data come from a site with more
than the two censuses used to run biomass.

The first function, biomass.change can divide AGB into 2 sets of mutually
exclusive categorical levels.  Do NOT use DBH classes with this function.
For DBH classes use the second function, biomass.change.dbh.  Here the
trees are categorized into dbh classes in each census independently of each
other.  Therefore, trees alive in the first census are categorized by their
census 1 dbh to determine their dbh class, etc.  Trees that died are
categorized by their census 1 dbh, recruits by their census 2 dbh.  Growth
in AGB is categorized into DBH classes based on census 1.

To obtain AGB Mg per year, compute the mean number of years between
censuses and apply to all values in the returned list.  For area, compute
total area in ha for plot and divide list contents by this also.}
\value{
  Returns a list of arrays:
	\item{alive1 }{AGB Mg of trees alive in census 1}
	\item{alive2 }{AGB Mg of trees alive in census 2}
	\item{dead }{AGB Mg of trees that died between census 1 and 2}
	\item{rec }{AGB Mg of trees that recruited between census 1 and 2}
	\item{aliveboth1 }{AGB Mg of census 1 for trees that were alive in both
	censuses}
	\item{aliveboth2 }{AGB Mg of census 2 for trees that were alive in both
	censuses}
	\item{alivePC1 }{AGB Mg of census 1 for trees that were alive in both
	censuses but for which the point of measurement was changed (POM)}
	\item{alivePC2 }{AGB Mg of census 2 for trees that were alive in both
	censuses but for which the point of measurement was changed (POM)}
	\item{nalive1 }{Number of trees alive in census 1}
	\item{nalive2 }{Number of trees alive in census 2}
	\item{ndead }{Number of trees that died between census 1 and 2}
	\item{nrec }{Number of trees that recruited between census 1 and 2}
	\item{naliveboth1 }{Number of census 1 for trees that were alive in both
	censuses}
	\item{naliveboth2 }{Number of census 2 for trees that were alive in both
	censuses}
	\item{nalivePC1 }{Number of census 1 for trees that were alive in both
	censuses but for which the point of measurement was changed (POM)}
	\item{nalivePC2 }{Number of census 2 for trees that were alive in both
	censuses but for which the point of measurement was changed (POM)}
	\item{alivegr }{Change in AGB for trees that were alive in both censuses.
	Does not include trees with a change in POM}
	\item{alivegrPC }{Change in AGB for trees that were alive in both
	censuses and had a change in POM. Since POM was changed, this AGB growth
	is not valid but is provided for completeness}
}
\references{\cite{Chave, J. et.al.(2005). Tree allometry 
	and improved estimation of carbon stocks and balance in tropical 
	forests. Oecologia. 145. pg 87}}
\author{Pamela Hall and Jerome Chave}
\seealso{\code{\link{CTFS.biomass}}}
\examples{\dontrun{
1.  Use result from example for biomass:

bci9095.biomass <- biomass(bci90.full,bci95.full,mult1=tst.bci90.mult,mult2=tst.bci95.mult,wsgfile=wsg.ctfs,
        forest="moist",badgrowth="meangrowth",dbhuse="second",trim=c(-5,75),recsize=110,rounddown=FALSE)
        
biomass.change.out <- biomass.change(biomass.out,split1=tst.bci90.full$sp)
biomass.dbh.out <- biomass.change.dbh(tst.bci90.full,tst.bci95.full,biomass.out,split=tst.bci90.full$sp,dbhclass=c(10,50,100,300))
}}
\keyword{univar}
