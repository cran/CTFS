\name{socrex90.neigh}
\alias{socrex90.neigh}
\concept{CTFS example dataset}
\concept{neighborhood}
\concept{density dependence}
\docType{data}
\title{Structure of Datasets for Density Dependence Analysis}

\description{ This file can be used to investigate the effect of
	neighborhood attributes (density and basal area) on the population
	structure and dynamics of focal species. The file contains the number and
	basal area of conspecific and heterospecific trees for given neighborhood
	distances (annuli) around each tree of a single species.  }

\usage{data(socrex90.neigh)}
\format{
  
A \strong{data frame} for the species \dQuote{socrex} with 1133
observations on the following 17 variables.

	\describe{
    \item{\code{tag}}{a numeric vector, tree identification number for a 
		tree of \dQuote{socrex}}
    \item{\code{Ncon.0.5}}{a numeric vector, number of conspecific trees 
		>=0 and < 5 m from the tag tree}
    \item{\code{Ncon.5.10}}{a numeric vector, number of conspecific trees 
		>=05 and < 10 m from the tag tree}
    \item{\code{Ncon.10.20}}{a numeric vector, number of conspecific trees 
		>=10 and < 20 m from the tag tree}
    \item{\code{Ncon.20.30}}{a numeric vector, number of conspecific trees 
		>=20 and < 30 m from the tag tree}
    \item{\code{BAcon.0.5}}{a numeric vector, total basal area of 
		conspecific trees >=0 and < 5 m from the tag tree}
    \item{\code{BAcon.5.10}}{a numeric vector, total basal area of 
		conspecific trees >=5 and < 10 m from the tag tree}
    \item{\code{BAcon.10.20}}{a numeric vector, total basal area of 
		conspecific trees >=10 and < 20 m from the tag tree}
    \item{\code{BAcon.20.30}}{a numeric vector, total basal area of 
		conspecific trees >=20 and < 30 m from the tag tree}
    \item{\code{Nhet.0.5}}{a numeric vector, number of heterospecific trees 
		>=0 and < 5 m from the tag tree}
    \item{\code{Nhet.5.10}}{a numeric vector, number of heterospecific trees 
		>=05 and < 10 m from the tag tree}
    \item{\code{Nhet.10.20}}{a numeric vector, number of heterospecific trees 
		>=10 and < 20 m from the tag tree}
    \item{\code{Nhet.20.30}}{a numeric vector, number of heterospecific trees 
		>=20 and < 30 m from the tag tree}
    \item{\code{BAhet.0.5}}{a numeric vector, total basal area of 
		heterospecific trees >=0 and < 5 m from the tag tree}
    \item{\code{BAhet.5.10}}{a numeric vector, total basal area of 
		heterospecific trees >=5 and < 10 m from the tag tree}
    \item{\code{BAhet.10.20}}{a numeric vector, total basal area of 
		heterospecific trees >=10 and < 20 m from the tag tree}
    \item{\code{BAhet.20.30}}{a numeric vector, total basal area of 
		heterospecific trees >=20 and < 30 m from the tag tree}
  }
}
\details{
 This data file is the result of \code{\link{neighbordens.1sp}}.  It was 
 created using the default neighborhood distances, gridsize and 
 plot dimensions:
 
 \code{socrex90.neigh <- neighbordens.1sp("socrex",tst.bci90.quad, tst.bci90.spp)}}

 \source{
 The Center for Tropical Forest Science\cr
\url{http://www.ctfs.si.edu}
}

\examples{\dontrun{
data(socrex90.neigh)
str(socrex90.neigh)
}}
\keyword{datasets}
