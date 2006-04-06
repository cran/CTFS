\name{badgrowth.tree}
\alias{badgrowth.tree}
\concept{biomass}
\concept{wood specific gravity}
\concept{wood density}

\title{Identify Trees for which DBH was Estimated for Biomass
Calculation}
\description{
 Provides information for each tree for which a DBH was estimated when
 computing biomass.  The DBHs and growth which was measured and
 estimated are provided in a dataframe.
}
\usage{
badgrowth.tree(census1, census2, badgrowth = "meangrowth", 
dbhuse = "second", trim = c(-5, 75), recsize = 110, rounddown = FALSE)
}
\arguments{
  \item{census1}{name \strong{dataframe} of first census, must be of
  same length as \code{census2}}
  \item{census2}{name \strong{dataframe} of second census, must be of
  same length as \code{census1}}
  \item{badgrowth}{parameter indicating how to handle trees with
  \dQuote{bad} growth. \dQuote{bad} growth is defined by \code{trim}.
  Valid values are: \dQuote{keep} for retaining 
  all bad growth rates without alteration, \dQuote{nogrowth} for
  setting bad growths to 0, \dQuote{meangrowth} for setting bad growth
  to the mean value of the dbh size class of the tree, \dQuote{trim}
  for using the values of \code{trim} to set growth to maximum and
  minium as appropriate.}
  \item{dbhuse}{indicate which \code{dbh} to use for computing growth 
  after implementing value of \code{badgrowth}.  Valid values are:
  \dQuote{second} for the dbh in the first census, \dQuote{second} for 
  dbh in the second census.}
  \item{trim}{minimum and maximum value of growth that is considered
  acceptable when \code{badgrowth} = \dQuote{trim}}
	\item{recsize}{the maximum DBH (mmm) a recruit to the second census that is 
	considered appropriate, a very site and species specific value that 
	should be carefully determined by the user.}
  \item{rounddown}{how to round down values of growth for trees enumerated 
  at Pasoh in 1990 due to a measurement standard. See
  \code{\link{CTFS.datafiles}}}
}


\details{
  This function computes the growth rate of each tree and determines
  if that growth rate is too high or too low based on user provided 
  values, \code{trim}. If the growth rate is out of these bounds, a
  dbh is estimated, based on a user provided choice, \code{dbhuse}. 
  These trees are returned in a dataframe for inspection to verify how
  well \code{biomass} is performing.
}
\value{
  Returns a \strong{dataframe} with a row for each tree with an
  estimated dbh.  The dbh that is estimated is determined by the user.
    \item{$tag}{an integer, tree identification number}
    \item{$sp}{character, species name code}
    \item{$gx}{a real number, location on East-West axis, in meters}
    \item{$gy}{a real number, location on North-South axis, in meters}
	\item{$dbh1}{a numeric vector, DBH in mm measured in year indicated}
 	\item{$dbh2}{a numeric vector, DBH in mm measured in year indicated}
	\item{$realGR}{a real number, growth rate of tree without any
	elimination of extreme values}
	\item{$estdbh1}{a numeric vector, estimated DBH in mm measured in year indicated}
 	\item{$estdbh2}{a numeric vector, estimated DBH in mm measured in year indicated}
	\item{$estlGR}{a real number, estimated growth rate of tree based
	on value of \code{badgrowth}}
	\item{$pom1}{a numeric vector, point of measurement, in m measured 
    in year indicated}
	\item{$pom1}{a numeric vector, point of measurement, in m measured 
    in year indicated}
	\item{$status1}{character, value indicating whether the tree is
	alive, dead or a recruit of that census.}
	\item{$status2}{character, value indicating whether the tree is
	alive, dead or a recruit of that census}
}

\author{Pamela Hall}
\note{
	The \code{estdbh} will depend on the value of \code{dbhuse}.
	\code{dbhuse} is the census dbh that is retained and used to compute \code{estdbh}.
}

\seealso{\code{\link{CTFS.biomass}}}
\examples{\dontrun{
bci.badtrees <- badgrowth.tree(tst.bci90.full, tst.bci95.full,
	badgrowth="meangrowth",dbhuse="second")
}
}
\keyword{univar}