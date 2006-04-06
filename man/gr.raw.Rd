\name{gr.raw}
\alias{gr.raw}
\concept{growth}
\concept{dynamics}

\title{Unadjusted, Untrimmed Annual DBH Increment for Individual Trees}

\description{
		Computes the annual dbh increment for each tree and compares it to an
		expected range of values based on the function
		\code{\link{trim.growth}}.  The range of values is set so high and low
		as to include all possible growth rates (or so it is believed).  This
		function is intended for use by \code{\link{biomass}} but can be used
		by itself if unadjusted, raw DBH increments for individuals trees is
		desired.
}
\usage{
gr.raw(census1, census2, rounddown = FALSE)
}

\arguments{
  \item{census1}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census2}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}

  \item{rounddown}{logical value in caps.  When TRUE, if either of 
  census is < 55, then the \code{floor} of the dbh value / 5 is 
  provided.  When FALSE, no change in the dbh is made.}
}
\details{
	See \code{\link{CTFS.growth}} for details on the computation methods 
	of growth rates and associated functions.

	Any two censuses on a datafile must be used. They do not have to be 
	sequential, only that \code{census1} has to be before \code{census2}.
	
	The results of \code{gr.raw} are annual DBH increments for individual
	trees, \strong{NOT MEAN GROWTH VALUES}. A value is returned for
	each tree.  No invalid values should result but this means that no
	reasonable trimming of growth rates has been done.
}
\value{
		\code{gr.raw} returns a vector of annual dbh increments.  It is not
		intended for use outside of \code{\link{biomass}} but it will function
		on its own.
		
		If you want the annual dbh increments for individual trees with the 
		full compliment of computed values, use \code{\link{growth.indiv}}.
}

\author{Pamela Hall}
\seealso{\code{\link{growth.indiv}}, \code{\link{biomass}},
		\code{\link{CTFS.growth}}}
\examples{\dontrun{

gr.raw.out <- gr.raw(tst.bci90.full,tst.bci95.full)

}
}

\keyword{univar}
