\name{tojulian}
\alias{tojulian}
\alias{fromjulian}
\concept{dates}
\concept{calendar dates}

\title{Converts Calendar Dates to Julian Dates and Back}
\description{
  These two functions can be used to convert dates to and from julian 
  formats for computations.  There are many functions that do this but
  these are useful for CTFS.  The functions rely upon functions found 
  in the R package \bold{survival}.
}
\usage{
tojulian(x, dateform = "\%m/\%d/\%Y", centurycut = 2020)
fromjulian(j, dateform = "\%m/\%d/\%Y")
}

\arguments{
  \item{x}{calendate date in the default format \dQuote{m/d/yyyy}.  The year
  should be specified with the century included, eg. \dQuote{12/4/1999}}
  \item{j}{the number of days since 1/1/1960}
  \item{dateform}{format for calendar date, see \code{strptime} for
  format information}
  \item{centurycut}{years after this date are assumed to be in the
  21st century}
}
\details{
 These functions are dependent upon the package \bold{survival}.  They
 use the functions \code{mdy.date} and \code{date.mmddyyy}.
 
 The results of \code{tojulian} is a number of days, so it is an integer and can be used for
 computation.
}
\value{
 \item{\code{tojulian}}{returns an integer, the number of days since 1/1/1960
 formatted as d/m/y where d = day, m = 3 letter abreviation for month 
 eg Jan, Feb, Mar\ldots, and y = year.}
 \item{\code{fromjulian}}{returns a calendar date in mm/dd/yyyy format}
}
\references{depends on the package \bold{survival}}
\author{Rick Condit}
\seealso{\code{\link{CTFS.utilities}}
}

\examples{\dontrun{
1.  tojulian()

x1="3/5/1987"
x2="3/5/2001"
y1 <- tojulian(x1)
y2 <- tojulian(x2)
str(y1)
str(y2)
y2-y1

2.  fromjulian()
fromjulian(9925)
fromjulian(15039)
fromjulian(y1+(y2-y1))
  }
}
\keyword{manip}
