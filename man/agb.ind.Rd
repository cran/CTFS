\name{agb.ind}
\alias{agb.ind}
\alias{agb.mult}
\alias{getwsg.vct}
\concept{biomass}
\concept{wood specific gravity}
\concept{wood density}

\title{Computes Above Growth Biomass for Individual Tree}
\description{
 This function is called from \code{biomass}.  This is the core function
 that uses the wood specific gravity data file and the coefficients from
 the model in Chave, J. et.al.  to compute the above ground biomass for
 each tree and multiple stem.  This package contains only wood specific
 gravity values from species found at the CTFS site, BCI. For access to the
 full data set of wood specific gravity, contact J. Chave at: chave@cict.fr
}
\usage{
agb.ind(dbh, wsg.vct, forest = forest)
agb.mult(census, multfile, wsgfile, forest = forest)
getwsg.vct(census, wsgfile, by.col = "wsg")

}
\arguments{
  \item{dbh}{dbh in mm of tree, determined through logic and
  parameters of \code{biomass}}
  \item{wsg.vct}{wood specific gravity for all species}
  \item{forest}{forest type which determines which set of coefficients
  are used for calculation}
  
  \item{census}{name \strong{dataframe} of a census}
  \item{multfile}{name of \strong{dataframe} with dbh measurements of
  secondary stems of census.  Tree tag numbers can appear
  multiple times in this file if there is more than two stems for a
  tree.}
  \item{wsgfile}{name of the \strong{dataframe} that contains wood
  density values for selected species. See \code{\link{tst.ctfs.wsg}}}
  
  \item{by.col}{column name of the wood specific gravity values
  in \code{tst.ctfs.wsg.Rdata}, usually \dQuote{wsg}}

}
\details{
  These functions are only used internally but can, with the correct
  parameter values, be used directly. \code{agb.ind} is the core
  function with the coefficients for the estimation of above ground
  biomass.  \code{agb.mult} calls \code{agb.ind} to compute the sum of
  the above ground biomass of all secondary stems of the same
  individual tree as defined by \code{tag}.  
  
  The current set of coefficients used are given below anda are taken
  from Table 4, II.3 and II.4.  \code{forest = NULL} is the same as
  \dQuote{all types} as described in Chave, J. et.al.
  
	\deqn{AGB = p * exp(aprime + b*lndbh + c*lndbh^2 + d*lndbh^3)}
	\deqn{p = wood specific gravity for species}
	\deqn{aprime = \frac{(a*RSE^2)}{2}}{%
	aprime = (a*RSE^2) / 2}
	
	\tabular{lrrrrr}{
	Forest \tab a \tab b \tab c \tab d \tab RSE\cr
	NULL \tab -1.589 \tab 2.284 \tab 0.129 \tab -0.0197 \tab 0.377\cr
	wet \tab -1.302 \tab 1.98 \tab 0.207 \tab -0.0281 \tab 0.356\cr
	moist \tab -1.562 \tab 2.148 \tab 0.207 \tab -0.0281 \tab 0.356\cr
	dry \tab -0.730 \tab 7.784 \tab 0.207 \tab -0.0281 \tab 0.356\cr
	}

  \code{getwsg.vct} matches the species name to the wood specific
  gravity for that species in \code{wsg.ctfs.RData}.
  }
	
\value{
 Returns a number, \code{AGB}, the above ground biomass for the dbh
 submitted.
}
\seealso{\code{\link{CTFS.biomass}}
}
\keyword{univar}