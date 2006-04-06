\name{sep.species}
\alias{sep.species}
\alias{load.species}
\concept{creating datasets}
\concept{species dataframes, separate}

\title{Separate One or More Dataframes of Tree Records from a List}
\description{
	Provides access to a single or groups of \strong{dataframes} in a
	\strong{list} where each dataframe is the tree records of a single
	species.  This list data file is generated with \code{\link{sep.data}}.  The
	example data file is \code{\link{tst.bci90.spp}}.
}
\usage{
sep.species(spname, census.spp)
load.species(spname, census.spp)
}

\arguments{
  \item{spname}{species names in quotes, either a single name or a vector 
	of names}
  \item{census.spp}{the list structured dataset where each species is a 
  dataframe, eg. \code{tst.bci90.spp} }
}
\details{
 These functions are similar to each other.  \code{load.species} only 
 works with one species name at a time and returns a data frame.  
 \code{sep.species} returns a list that is a subset of data frames in \code{census.spp}.
}
\value{
  For \code{sep.species}, if \code{spname} is a single name then a data frame of tree records is 
	returned.  If \code{spname} is a vector of names, then a list of 
	data frames is returned, but only those data frames of the species names 
	in \code{spname}.  

	For \code{load.species}, \code{spname} can only be a single value and a 
	data frame is always returned.
}

\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.readwrite}}}

\examples{\dontrun{
socrex <- load.species("socrex",tst.bci95.spp)
socrex[1:4,]

socrex <- sep.species("socrex",tst.bci95.spp)
socrex[1:4,]

Using sep.species or load.species with one species name is the same as
socrex <- tst.bci90.spp$socrex
socrex[1:4,]

spnames <- c("alsebl","socrex")
alsebl.socrex <- sep.species(spnames,tst.bci95.spp)
is.list(alsebl.socrex)
str(alsebl.socrex)
}}
\keyword{manip}
