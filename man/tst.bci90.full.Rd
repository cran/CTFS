\name{tst.bci90.full}
\alias{tst.bci90.full}
\alias{tst.bci95.full}
\alias{tst.bci90.spp}
\alias{tst.bci95.spp}
\alias{tst.bci90.quad}
\alias{tst.bci95.quad}
\concept{CTFS example dataset}
\docType{data}
\title{Structure of CTFS Long Term Forest Data Set for a Single Census}
\description{
  The structure of the CTFS Long Term Forest data set for all trees
  with census information from a single census including variable
  names, values for categorical variable and the contents of a test
  file for demonstration of functions.
}
\usage{data(tst.bci90.full)}
\format{
		A \strong{data frame} with 12421 observations on the following 9 
		variables.
  \describe{
  \item{tag}{a numeric vector, tree identification number}
  \item{sp}{character, species name code}
	\item{gx}{a numeric vector, location on East-West axis, in meters}
  \item{gy}{a numeric vector, location on North-South axis, in 
    meters}
	\item{dbh}{a numeric vector, DBH in mm measured in year indicated}
  \item{pom}{a numeric vector, Point of Measurement, in m measured 
    in year indicated}
	\item{date}{a numeric vector, Date of measurement of quadrat in 
    days since initial measurement of plot since year indicated}
	\item{codes}{character, various values which indicate condition of 
    tree such as broken stems, etc. in year indicated}
	\item{status}{character, value indicating whether the tree is
	alive, dead or a recruit of that census}
}}

\details{
		
	\code{tst.bci90.full} and \code{tst.bci95.full} are subsets of the
	complete data files for BCI census 1990 and 1995, respectively. 

  \code{tst.bci90.spp}, \code{tst.bci90.quad} and \code{tst.bci95.spp},
	\code{tst.bci95.quad} are the same data as the \code{*.full} files
	restructured as \strong{a list of dataframes}.  For \code{*.spp} each
	data frame consists of the records of a single species.  For
	\code{*.quad} each data frame consists of the record for a single
	quadrate (20 x 20).  The names on the list are the species names and
	quadrate number, respectively.  The \code{*.quad} file is used for the
	density dependent analysis.  They are created using
	\code{\link{sep.data}} for \code{*.spp} and
	\code{\link{sep.vct.data}} for for \code{*.quad} files.
	
	The variable \code{codes} can have many different values which have 
	different meaning at each site.  It is best not to use this variable 
	because of its inconsistency across sites.
	}
	
\section{Contents ot Test Files}{	

	Only three species are included in these test files: \strong{alsebl},
	\strong{psycde}, and \strong{socrex}.\cr 
  
  the species names and number of trees for each are:\cr
	alsebl:`data.frame':	11128 obs. of  9 variables:\cr
	psycde:`data.frame':	160 obs. of  9 variables:\cr
	socrex:`data.frame':	1133 obs. of  9 variables:\cr
	}

\section{Variable Values}{
  \tabular{ll}{
  POM \tab Point of Measurement\cr
	0  \tab at initial position\cr
	1  \tab first change of POM\cr
	2 \ldots \tab second to n change of POM\cr
	
  STATUS \tab Values\cr
	\dQuote{A} \tab alive\cr
	\dQuote{AB} \tab alive but broken\cr
	\dQuote{AS} \tab alive but broken, used at the site Pasoh\cr
	\dQuote{D} \tab dead\cr
	\dQuote{P} \tab not yet large enough for measurement, future recruit\cr
	NA \tab missing data, unknown cause\cr
	}
	Only trees with STATUS == \dQuote{A} have valid dbh measurements.
	All other trees have a dbh = NA.
}

\source{
 The Center for Tropical Forest Science\cr
\url{http://www.ctfs.si.edu}
}

\examples{\dontrun{
data(tst.bci90.full)
str(tst.bci90.full)
data(tst.bci90.spp)
str(tst.bci90.spp)
}
}
\keyword{datasets}
