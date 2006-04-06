\name{abund.ch.prob}
\alias{abund.ch.prob}
\concept{dynamics}
\concept{abundance}
\concept{population size}

\title{Calculation of Probability of Change in Population size}
\description{
  The probability that N0 = N1 is computed.
	}
\usage{
abund.ch.prob(N0, N1, recrate, survrate)
}
\arguments{
  \item{N0}{abundance in first census}
  \item{N1}{abundance in second census}     
  \item{recrate}{number of recruits / N1}
  \item{survrate}{number of survivors / N0}
}
\details{
	See \code{\link{CTFS.abundance}} for details on the computation of 
	abundance and associated functions.
	
	This function is designed to be called from
	\code{\link{abundance.change.calculation}}.  It can be used independently
	if the number of recruits and survivors are known.  The probability
	returned is the sum of a biniomial (survivorship) and poisson
	(recruitment) distribution.  The arguments can be vectors of any length.
	This function does not work with multi-dimensional arrays.  To use if for
	these, use \code{\link{abundance.change.calculation}}
	
}
\value{
  \item{prob }{\strong{vector} of probabilities of same length as the arguments}
}

\author{Rick Condit}
\seealso{\code{\link{CTFS.abundance}}, \code{\link{abundance.change}}, 
\code{\link{abundance.change.calculation}}
}
\keyword{htest}
