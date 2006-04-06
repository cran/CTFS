\name{bcispp.info}
\alias{bcispp.info}
\concept{species}
\concept{latin name}
\concept{growth form}
\concept{breeding system}
\concept{reproductive size}
\concept{maximum height}

\docType{data}
\title{Dataset of Information on each Species found in the BCI 50 HA plot}
\description{
 A tabulation of information on each species found during all cenuses 
 of the BCI 50 HA plot.  Information currently includes taxonomic 
 levels, growth form, size at reproductions, breed system and maximum 
 height.  The information for each species is derived from various 
 sources.
}
\usage{data(bcispp.info)}
\format{
  A \strong{dataframe} with 318 rows and 8 columns
  \describe{
		\item{rowname}{a character vector, 6 letter code for a species name: 
    genus is 2 letters and species is 4.  Names with numbers indicate 
    species that have yet to be formally named.  Names may be upper or 
		lower case but \code{rowname} and \code{sp} are of the same case.}
    \item{sp}{a character vector, identical to rowname}
    \item{genus}{a character vector, full latin Genus name}
    \item{species}{a character vector, full latin Species name}
    \item{family}{a character vector, full latin Family name}
    \item{grform}{a character vector}
    \item{repsize}{a numeric vector}
    \item{breedsys}{a character vector}
    \item{maxht}{a numeric vector}
  }
}
\details{
 \code{bcispp.info} variables have the following values:
 \describe{
	\item{grform}{the growth form of a species, derived from
	knowledge of the species at this site and elsewhere.}\cr
 	\tabular{lcr}{
	Code  \tab Growth Form  \tab Number of Species \cr
	H  \tab herb \tab 2\cr
	M  \tab main canopy \tab 75 \cr
	S  \tab shrub \tab 59 \cr
	T  \tab treelet \tab 108 \cr
	U  \tab understorey \tab 64 \cr
	Z  \tab unknown \tab 4 \cr
	\dQuote{.} \tab other missing \tab 6 \cr
	}
	\item{breedsys}{the breeding system of the species,
	derived from knowledge of the species at this site an elsewhere.}
	
	\tabular{lcr}{
	Code  \tab Breeding System \tab Number of Species \cr
	B  \tab bisexual \tab 202 \cr
	D  \tab dioecious \tab 64 \cr
	M  \tab monecious \tab 35 \cr
	P  \tab polygamous \tab 11 \cr
	\dQuote{.} \tab other missing \tab 6 \cr
	}
	
	\item{repsize}{the size of the tree upon reproduction, DBH
	in cm, derived from knowledge of the species at this site and
	elsewhere.}
	\item{maxht}{the maximum height achieved by a species,
	derived from knowledge of the speices at this site and
	elsewhere.}
	}
}
\source{
The Center for Tropical Forest Science\cr
\url{http://www.ctfs.si.edu}
}
\examples{\dontrun{
str(bcispp.info)
`data.frame':	318 obs. of  8 variables:
 $ sp      : chr  "acacme" "acaldi" "acalma" "ade1tr" ...
 $ genus   : chr  "Acacia" "Acalypha" "Acalypha" "Adelia" ...
 $ species : chr  "melanoceras" "diversifolia" "macrostachya" "triloba" ...
 $ family  : chr  "Fabaceae:Mimos." "Euphorbiaceae" "Euphorbiaceae" "Euphorbiaceae" ...
 $ grform  : chr  "U" "S" "U" "U" ...
 $ repsize : int  4 2 2 10 4 20 16 4 4 20 ...
 $ breedsys: chr  "B" "M" "M" "D" ...
 $ maxht   : chr  "6" "6" "5" "5" ...
}
}
\keyword{datasets}
