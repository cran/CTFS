\name{bciquad.info}
\alias{bciquad.info}
\concept{quadrate}
\concept{habitat}
\concept{elevation}
\concept{convexity}
\docType{data}

\title{Dataset of Information for each 20 by 20 m quadrate in the BCI 50 ha plot}
\description{
  A tabulation of information on each quadrate that is dervied from 
  the initial establishment of the plot (elevation, convexit and 
  slope) and from further analyses of these and other topographic or 
  local environmental features.  These characteristics are 
  independent of the trees that occupy the quadrates and their 
  species identification.
}
\usage{data(bciquad.info)}
\format{
		A data frame with 1250 observations on the following 5 variables.
		\describe{ 
		\item{elev}{a numeric vector, elevation, msl}
		\item{convex}{a numeric vector} \item{slope}{a numeric vector}
		\item{hab}{a numeric vector, a value from 1 to n, assigned by
		investigator.  Combinations of altitude, slope and convexity, mapped 
		location of streams, swamps, etc. can be used to define habitats}
		\item{quad}{a factor with levels 1250 in the standard 50 Ha Plots}
}}

\details{ 
	\code{bciquad.info} variables have the following values:
	\describe{
		\item{hab}{the habitat assigned to a quadrate by topographic values.
		These values have been assigned by the habitat analysis function
		associated with the "torus" analysis.} 
\tabular{llr}{
	Code  \tab Habitat\_Type \tab Number\_Quadrates \cr
	1  \tab low\_plateau \tab 168 \cr
	2  \tab hi\_plateau \tab 81 \cr
	3  \tab low\_slope \tab 277 \cr
	4  \tab hi\_slope \tab 99 \cr
	5 \tab low\_hill \tab 173 \cr
	6 \tab hi\_hill \tab 203 \cr
	7 \tab low\_hilly\_slope \tab 129 \cr
	8 \tab hi\_hilly\_slope \tab 120 \cr
	\dQuote{.} \tab missing \tab 0 \cr
	}
\item{quad}{an "index" number for each quadrates. For this 
  BCI dataset a quadrate is 20 by 20 m and there are 1250 of them in a 
 1000 by 500 m plot.  The convention is to number the quadrates starting in 
 the southwest corner and going northward up a "column" of 25, and 
 returning to the southern border for 26-50, etc.}
 
\tabular{rrrr}{
         25\tab 50 \tab 75 \ldots \tab 1250\cr
         24\tab 49 \tab 74  \ldots \tab 1249\cr
         . \tab . \tab .  \ldots \tab .\cr
         . \tab . \tab .  \ldots \tab . \cr
         . \tab . \tab .  \ldots \tab . \cr
         3\tab 28 \tab 53  \ldots \tab 1228\cr
         2\tab 27 \tab 52  \ldots \tab 1227\cr
         1\tab 26 \tab 51  \ldots \tab 1226\cr
 }}}
	 
\source{
 The Center for Tropical Forest Science\cr
\url{http://www.ctfs.si.edu}
}
%\references{  articles?}
\examples{\dontrun{
str(bciquad.info)
}
}
\keyword{datasets}
