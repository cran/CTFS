\name{create.status}
\alias{create.status}
\concept{census}
\concept{status}
\concept{create data file}

\title{Create STATUS variable for CTFS Datasets }
\description{Creates the \code{status} variable based on the value of
\code{dbh} in CTFS datasets of a single census.}

\usage{
create.status(census, pasoh = FALSE)
}

\arguments{
  \item{census}{CTFS datafile of a single census which does not
  contain \code{status}}
  \item{pasoh}{logical, FALSE if datafile is not from Pasoh site}
}
\details{
	\code{status} indicates whether a tree is alive, dead, or a
	recruit for a single census.  Trees that were broken below the
	point of measurement or shrunk to less than 10 mm are indicated as
	alive but broken and alive but too small, respectively.  At Pasoh
	a different dbh value for alive but broken was used than at the
	other sites.  \code{dbh} is set to NA if \code{status} does not
	equal \dQuote{A}.

\tabular{llll}{
\code{status} \tab \tab \code{dbh} \tab condition\cr
\dQuote{A} \tab \tab dbh >= 10 mm \tab alive, dbh measured and value
retained\cr
\dQuote{AS} \tab \tab dbh > 0 & dbh < 10 \tab alive, shrunk, dbh = NA\cr
\dQuote{AB} \tab \tab dbh == 0 \tab alive but broken, do not alter recorded dbh\cr
\dQuote{AB} \tab \tab dbh == 5 \tab for Pasoh: alive, but dbh = NA\cr
\dQuote{P} \tab \tab dbh == -2 \tab recruit in future census, dbh = NA\cr
\dQuote{D} \tab \tab dbh == -1 \tab dead in this and/or past census, dbh=NA\cr
NA \tab \tab dbh == -9 \tab no information, dbh = NA\cr
}
}
\value{
	returns a data frame of the same number of rows as the input files 
	and with all of the columns of input file plus \code{status} 
	}
\author{Suzanne Lao and Pamela Hall}

\seealso{\code{\link{CTFS.readwrite}}}
\examples{
}
\keyword{univar}

