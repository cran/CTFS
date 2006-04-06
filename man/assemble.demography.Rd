\name{assemble.demography}
\alias{assemble.demography}
\concept{demography}
\concept{dynamics}
\concept{reformat}

\title{Reformat the Output from Demographic Functions from List to 
Dataframe}
\description{
	This function reformats the \strong{list of lists} output from the 
	demographic functions \code{mortality}, \code{recruitment}, 
	\code{growth} and \code{abundance.change}.  	
}
\usage{
assemble.demography(output, type = "g", whichdbhcat = 1, 
    datefirst = "1980-1-1")
}
\arguments{
  \item{output}{file containing the results from one of the 
  demographic functions}
  \item{type}{type of demographic function with values:
	\tabular{ll}{
  	\dQuote{g} \tab output from \code{\link{growth}}\cr
		\dQuote{a} \tab output from \code{\link{abundance.change}}\cr
		\dQuote{m} \tab output from \code{\link{mortality}}\cr
		\dQuote{r} \tab  output from \code{\link{recruitment}}\cr
}}
  \item{whichdbhcat}{number indicating how many dbh categories to 
  include, cannot exceed the number of categories in 
  \code{output}}
	\item{datefirst}{date of first enumeration at site, format: yyyy-mm-dd}
}


\details{
 The results of the demographic functions (\code{growth, mortality, 
 recruitment, abundance} and related functions are 
 \strong{lists}.  This is a very convenient format for further 
 analysis but is not the easiest way to view the results of these 
 functions.  \code{assemble.demography} reformats the list output 
 into a dataframe so that is resembles a table.  
 
 The columns of the dataframe are always composed of the values
 computed by the demographic function.  The rows are the categories 
 provided by the first split variable used in the demographic 
 function.  If a second split varible is used, then the argument 
 \code{whichdbhcat} determines which category of the second spit 
 variable is put into the dataframe.

 Returned Values
 
  \code{assemble.demography} returns a \strong{dataframe}.\cr
  
 	For \code{growth} the following values are returned:
	\describe{
	\item{\code{$rate}}{the annual growth rate in mm per year or \% per year,
  type determined by the argument \code{method} are provided for \code{growth}}
  \item{\code{$clim}}{the 95\% confidence limit based on a normal 
  distribution}
  \item{\code{$N0}}{the number of trees used to compute growth rate}
  \item{\code{$dbhmean}}{the mean dbh in mm for trees used to compute growth 
  rate}
  \item{\code{$meanyrs}}{mean number of years between census for trees used in 
  growth rate computation.}
  \item{\code{$date0}}{mean date of first census}
  \item{\code{$date1}}{mean date of second census}\cr}
	
  For \code{mortality} the following values are returned:
  \describe{
	\item{\code{$rate}}{the mortality in \%/year}
  \item{\code{$upper}}{the upper 95\% confidence interval based on a normal
  distribution}
  \item{\code{$lower}}{the lower 95\% confidence interval based
  on a normal distribution}
  \item{\code{$N0}}{the number of living trees at the first census}
	\item{\code{$S}}{the number of survivors from the first to second census}
	\item{\code{$dbhmean}}{the mean dbh in mm for trees used to compute growth 
  rate}
  \item{\code{$meanyrs}}{mean number of years between census for trees used in 
  growth rate computation.}
  \item{\code{$date0}}{mean date of first census}
  \item{\code{$date1}}{mean date of second census}\cr}
	
	For \code{recruitment} the following values are returned:
	\describe{
	\item{\code{$rate}}{the recruitment in \%/year}
  \item{\code{$upper}}{the upper 95\% confidence interval based on a normal
  distribution}
  \item{\code{$lower}}{the lower 95\% confidence interval based
  on a normal distribution}
  \item{\code{$N1}}{the number of living trees at the second census}
	\item{\code{$R}}{the number of recruits from the first to second census}
  \item{\code{$meanyrs}}{mean number of years between census for trees used in 
  growth rate computation.}
  \item{\code{$date0}}{mean date of first census}
  \item{\code{$date1}}{mean date of second census}\cr}

	For \code{abundance.change} the following values are returned:
  \describe{
	\item{\code{$rate}}{the change in abundance in \%/year}
  \item{\code{$prob}}{the prob \code{N0} == \code{N1}}
  \item{\code{$N0}}{the number of living trees at the first census}	
  \item{\code{$N1}}{the number of living trees at the second census}
	\item{\code{$D}}{the number of trees that died from the first to second census}
	\item{\code{$R}}{the number of recruits from the first to second census}
  \item{\code{$meanyrs}}{mean number of years between census for trees used in 
  growth rate computation.}
  \item{\code{$date0}}{mean date of first census}
  \item{\code{$date1}}{mean date of second census}\cr}
	}
	
\author{Rick Condit and Pamela Hall}
\seealso{ \code{\link{assemble.table.each}}}
\examples{\dontrun{
1. Run the demographic function

dbh.vct <- sep.dbh(tst.bci90.full,dbhclass=c(10,50,100,200,300,400,500))
growth.dbh <- growth(tst.bci90.full,tst.bcio95.full)
growth.dbh.out

2. Reformat the output into a dataframe.  Note there is only one row in 
this dataframe because no split variable was used in growth()

assemble.demography(growth.out)

3.  Run the demographic functions with two split variables, species 
names and dbh classes. Only the first list of the output is shown.
}
}
\keyword{manip}
