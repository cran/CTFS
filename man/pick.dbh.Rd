\name{pick.dbh}
\alias{pick.dbh}
\concept{biomass}
\concept{wood specific gravity}
\concept{wood density}
\concept{valid dbh}
\concept{unreasonable growth}

\title{Determine the Appropriate DBH to Use for Computation of Above
Ground Biomass}
\description{
   These functions are called from \code{biomass}.  Given the parameter
   values provided by \code{biomass} the dbh of tree with a suspect
   growth rate is set in order to have a reasonable estimate of the
   growth rate.  The growth rate estimations can have a very large
   effect on the final above ground biomass for the entire site. Use 
	 \code{\link{badgrowth.tree}} to generate a file of the trees with 
	 estimated dbh as determined by \code{pick.dbh}.
}
\usage{
pick.dbh(census1, census2, badgrowth = "keep", 
	dbhuse = "second", trim = c(-5, 75), recsize = 110, rounddown = FALSE)	
}

\arguments{
  \item{census1}{name \strong{dataframe} of first census, must be of
  same length as \code{census2}}
  \item{census2}{name \strong{dataframe} of second census, must be of
  same length as \code{census1}}
  \item{badgrowth}{how to treat growth rates that are determined as
  too high or low by \code{trim} values.  Valid values are: \code{keep}
  for retaining all growth rates regardless of value, i.e. \code{trim}
  limits are not applied; \code{nogrowth} set growth rates to 0;
  \code{meangrowth} set growth rate to the mean growth for the dbh
  category of the tree; \code{trim} set growth values to \code{trim}
  maximum or minimum if it exceeds or is less than the values for
  \code{trim}.} 
  \item{dbhuse}{for trees with \dQuote{bad} growth, identify which dbh
  to trust, then use this dbh to apply the rules selected for
  \code{\link{badgrowth.tree}}.}
  \item{trim}{limits of acceptable growth minimum and maximum}
	\item{recsize}{the maximum DBH (mmm) a recruit to the second census that is 
	considered appropriate, a very site and species specific value that 
	should be carefully determined by the user.}
  \item{rounddown}{rounding rules for the first and second censuses at
  Pasoh}
  }
\details{
 This function is only called by \code{calc.neighborhood.sp} and is
 used internally.  
}
\value{
	\code{pick.dbh} returns a \strong{list} of:
    \item{dbh1}{dbh in census1 as determined by \code{pickdbh}}
    \item{dbh2}{dbh in census2 as determined by \code{pickdbh}}
    \item{GRm}{growth (mm/year) of the main stem}

	\code{gr.raw} returns a vector of growth rates without adjustment.
	}

\author{Pamela Hall}
\seealso{\code{\link{CTFS.biomass}}, \code{\link{biomass}}}
\keyword{univar}

