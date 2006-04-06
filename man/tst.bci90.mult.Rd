\name{tst.bci90.mult}
\alias{tst.bci95.mult}
\alias{tst.bci90.mult}
\concept{BCI}
\concept{CTFS dataset}
\concept{multiple stems}
\docType{data}

\title{Structure of the CTFS Long Term Forest Data Set Multiple Stem
dataset for a single census}
\description{
	The structure of the CTFS Long Term Forest data set of
	measurements of multiple stems.  In addition to the measurement of
	a main stem, if a tree has other stems that are greater than 10 mm
	in dbh, they are also measured for diameter (DBH) only.  These
	measurements are recorded in this file.  These measurements are
	useful for computations of basal area and biomass but are not
	normally used for computations of indiviual tree dynamics.
}
\usage{data(tst.bci90.mult)}
\format{
	\describe{
		\item{tag }{a numeric vector, tree identification number}
		\item{dbh }{a numeric vector, DBH in mm measured in year indicated}
		}
}		
\details{
	\code{tst.bci90.mult} and \code{tst.bci95.mult} are subsets of the
	complete multi stem data files for BCI census 1990 and 1995,
	respectively.  \code{tst.bci90.mult} contains 584 records.
	\code{tst.bci95.mult} contains 548 records.  Only the the multiple
	stem measures of trees in the \code{tst.bci90.full} and
	\code{tst.bci95.full}, respectively are included.  The \code{tag}
	number is repeated for as many multiple stems as were measured for
	the main stem.\cr

}

\source{
 The Center for Tropical Forest Science\cr
\url{http://www.ctfs.si.edu}
}

\examples{\dontrun{
str(tst.bci90.mult)
}
}
\keyword{datasets}
