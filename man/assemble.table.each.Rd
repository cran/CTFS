\name{assemble.table.each}
\alias{assemble.table.each}
\concept{demography}
\concept{reformat}

\title{Reformat the Output from Demographic Functions from List to 
Dataframe}
\description{
	This function reformats the \strong{list of lists} output from the 
	demographic functions \code{mortality}, \code{recruitment}, 
	\code{growth} and \code{abundance.change}.  This function is designed to 
	provide separate tables for each dbh class and only returns 
	only rate, confidence limits and Ns. It is currently only implemented 
	for growth and mortality.  An alternative function that has similar 
	results is \code{\link{assemble.demography}}
}
\usage{
assemble.table.each(output, type = "g", whichdbhcat = 1:2)
}
\arguments{
  \item{output}{file containing the results from one of the 
  demographic functions}
  \item{type}{type of demographic function with values:
	\tabular{ll}{
  	\dQuote{g} \tab output from \code{growth}\cr
		\dQuote{a} \tab output from \code{abundance.change}\cr
		\dQuote{m} \tab output from \code{mortality}\cr
		\dQuote{r} \tab  output from \code{recruitment}\cr
}}
  \item{whichdbhcat}{number indicating how many dbh categories to 
  include, cannot exceed the number of categories in 
  \code{output}}
}

\details{
 The results of the demographic functions (\code{growth, mortality, 
 recruitment, totalabundance} and related functions are 
 \strong{lists}.  This is a very convenient format for further 
 analysis but is not the easiest way to view the results of these 
 functions.  \code{assemble.table.each} reformats the list output 
 into a dataframe so that is resembles a table.  
 
 This function only works for one category for computation of growth
 rates.  It is most often used for DBH classes.\cr
 
 Returned Values
 
  \code{assemble.table.each} returns a \strong{dataframe}.\cr
  
  For \code{growth} the following values are returned:
  \describe{
	\item{\code{$rate}}{annual growth rate in mm per year or \% per year,
  type determined by the argument \code{method} are provided for \code{growth}}
  \item{\code{$clim}}{95\% confidence limit based on a normal 
  distribution}
  \item{\code{$N0}}{number of trees used to compute growth rate}\cr}

  For \code{mortality} the following values are returned:
	\describe{
  \item{\code{$rate}}{mortality rate in \%/year}
  \item{\code{$upper}}{upper 95\% confidence interval based on a normal
  distribution}
  \item{\code{$lower}}{lower 95\% confidence interval based
  on a normal distribution}
  \item{\code{$N0}}{number of living trees at the first census}
  \item{\code{$S}}{number of surviving trees at the second census}
  \item{\code{$meanyrs}}{mean number of years between censuses}\cr}
}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{assemble.demography}}}
\examples{\dontrun{
dbh0.vct=sep.dbh(tst.bci90.full)
growth.dbh.out <- growth(tst.bci90.full,tst.bci95.full,split1=dbh0.vct)
assemble.demography(growth.dbh.out)

growth.sp.dbh.out <- 
growth(tst.bci90.full,tst.bci95.full,split1=tst.bci90.full$sp,split2=dbh0.vct)
assemble.table.each(growth.sp.dbh.out,whichdbhcat=1:3)

}}
\keyword{manip}
