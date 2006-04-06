\name{biomass}
\alias{biomass}
\concept{biomass}
\concept{wood specific gravity}
\concept{wood density}

\title{Compute Above Ground Biomass for Trees for 2 Censuses using Wood Specific Gravity for
Species}
\description{
  Computes the above ground biomass for individual trees for 2
  censuses including biomass of main stem, any and all multiple stems 
  and total for tree.  Growth rate used for the computation is also
  provided.  The biomass is computed using the models in Chave, J.
  et.al. (2005).
}
\usage{
biomass(census1, census2, mult1 = NULL, mult2 = NULL, wsgfile = NULL, 
	forest = NULL, badgrowth = "meangrowth", dbhuse = "second", 
	trim = c(-5, 75), recsize = 110, rounddown = FALSE)
}

\arguments{
  \item{census1}{name \strong{dataframe} of first census, must be of
  same length as \code{census2}}
  \item{census2}{name \strong{dataframe} of second census, must be of
  same length as \code{census1}}
  \item{mult1}{name of \strong{dataframe} with dbh measurements of
  secondary stems of census1.  Tree tag numbers can appear
  multiple times in this file if there is more than two stems for a
  tree.}
  \item{mult2}{name of \strong{dataframe} with dbh measurements of
  secondary stems of census2.  Tree tag numbers can appear
  multiple times in this file if there is more than two stems for a
  tree.}
  \item{wsgfile}{name of the \strong{dataframe} that contains wood
  density values for selected species. See \code{\link{tst.ctfs.wsg}}
}
  \item{forest}{parameter indicating forest type which affects the
  coeficients used to compute above ground biomass.  Valid values are 
  \dQuote{wet}, \dQuote{moist}, \dQuote{dry} and NULL.  The last
  results in a computation based on coefficients based on all forest
  types, a \dQuote{mean} value.}
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
  \dQuote{first} for the dbh in the first census, \dQuote{second} for 
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
  This function computes above ground biomass relying upon
  coefficients and models presented in Chave, J. et.al. (2005) (see
  references below). The function can be used to explore the impact of
  different choices for handling mismeasured dbh which affect growth
  rates.  See \code{\link{CTFS.biomass}} for further details.
	
	This function makes NO selection of records.  The returned object is of the
	same length and structure as the input file.  Use this object for further
	analysis of biomass for instance with \code{\link{biomass.change}}.
	
}

\value{
	Returns a \strong{dataframe} with row of \code{tag} numbers and
	the following columns:
	\item{tag}{tree tag number, same as in CTFS main dataset}
	\item{errtype}{indication of whether a DBH was estimated for AGB, values 
	are:}
	\tabular{ll}{
	\tab \dquote{NE} \tab no error \cr
  \tab \dquote{FG} \tab growth too fast, > trim[2] \cr
	\tab \dquote{SG} \tab growth too slow, < trim[1] \cr
  \tab \dquote{M1} \tab missing dbh1 even though status1 == A \cr
  \tab \dquote{M2} missing dbh2 even though status2 == A \cr
  \tab \dquote{LR} too large for a recruit \cr
	\tab \dquote{PC} POM changed - not technically an error, but tree should 
	be used with caution to compute AGB growth for trees that were alive in 
	both censuses\cr
  }
	\item{AGBm1}{above ground biomass of main stems in census1}
	\item{AGBm2}{above ground biomass of main stems in census2}
	\item{AGBs1}{summed above ground biomass of all secondary stems in census1}
	\item{AGBs2}{summed above ground biomass of all secondary stems in
	census2}
	\item{AGBt1}{summed above ground biomass of main and secondary stems in census1}
	\item{AGBt2}{summed above ground biomass of main and secondary stems in census1}
}
\references{\cite{Chave, J. et.al.(2005). Tree allometry 
	and improved estimation of carbon stocks and balance in tropical 
	forests. Oecologia. 145. pg 87}}
\author{Pamela Hall}
\note{This function calls the following functions which are specific
for \code{biomass}:
\code{\link{pick.dbh}}, \code{\link{getwsg.vct}}, \code{\link{agb.ind}},
\code{\link{agb.mult}}, \code{\link{gr.raw}}, \code{\link{growth.dbh}}, 
\code{\link{trim.growth}} and \code{\link{badgrowth.tree}}
}



\seealso{\code{\link{CTFS.biomass}}}
\examples{\dontrun{
1.  AGB in 1990 and 1995 for BCI which is a moist forest, use mean growth
for dbh class as replacement for bad growth rates and rely upon the dbh in
1995 (second).

bci9095.biomass <- biomass(bci90.full,bci95.full,mult1=tst.bci90.mult,mult2=tst.bci95.mult,wsgfile=wsg.ctfs,
	forest="moist",badgrowth="meangrowth",dbhuse="second",trim=c(-5,75),recsize=110,rounddown=FALSE)
  }
}
\keyword{univar}
