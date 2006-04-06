\name{tst.ctfs.wsg}
\alias{tst.ctfs.wsg}
\concept{wood specific gravity}
\concept{wood density}
\concept{biomass}

\docType{data}
\title{Wood Specific Gravity for Tropical Tree Species}
\description{
	This dataset contains information a subset of the full data set of
	the wood specific gravity (\code{wsg}) for tropical tree species for
	calculation of above ground biomass.  Only species found at the CTFS
	site, BCI, are included in the file distributed with this package.
	Species name codes as found in CTFS datasets, full latin name for
	Family, Genus and Species, the CTFS site the species has been
	recorded at and the level of identification used to assign a wood
	specific gravity.
}
\usage{data(wsg.ctfs)}
\format{
  A data frame with 5761 observations on the following 7 variables.
  \describe{
    \item{\code{wsg}}{a numeric vector, wood specific gravity}
    \item{\code{idlevel}}{a character vector, level of matching wsg 
		with species name}
    \item{\code{site}}{a character vector, CTFS site for species}
    \item{\code{sp}}{a character vector, CTFS species name code}
    \item{\code{genus}}{a character vector, latin genus}
    \item{\code{species}}{a character vector, latin species}
    \item{\code{fam}}{a character vector, latin family}
}
}
\details{
  This data file is used for the analysis of above ground biomass in
  particular with the function \code{biomass}.  It is the parameter
  \code{wsgfile}.
}
\source{ \cite{Chave, J. et.al.(2005). Tree allometry 
	and improved estimation of carbon stocks and balance in tropical 
	forests. Oecologia. 145. pg 87}}
}
}
\examples{\dontrun{
str(tst.ctfs.wsg)
}}
\keyword{datasets}

