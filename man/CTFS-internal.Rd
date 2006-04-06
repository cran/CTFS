\name{CTFS.internal}
\alias{borderdist}
\alias{fill.dimension}
\alias{fill.1dimension}
\alias{fill.distclasses}
\alias{unidentified.species}
\alias{unwind.matrix}
\alias{rndown5}
\alias{xydist}
\alias{xydistvect}
\alias{setsize}
\alias{quadslope}
\alias{calcslope}
%\alias{showci}

\title{internal CTFS functions}
\description{
	Internal CTFS functions that are called from within other functions.
	}

\usage{
borderdist(coord,plotdim)
fill.dimension(dataarray,class1,class2,fill=0)
fill.1dimension(dataarray,class1,fill=0)
fill.distclasses(dataarray,distclasses)
unidentified.species(spnames,unidenstr=c("**","UNID","unid","UNIDEN","uniden"))
unwind.matrix(matrix.in)
rndown5(s)
xydist(x1,y1,x2,y2)
xydistvect(pt1,pt2)
setsize(n, s) 
quadslope(cornerelev,gridsize=20)
calcslope(z, gridsize=20)
%showci(vector,number,number,character vector)
}

\details{
 These are not usually called by the user and are not documented
 individually.\cr  
 \item{\code{borderdist}}{computes the distance to the nearest plot border
 for a given pair of coordinates}\cr
 \item{\code{fill.dimension}}{check values and/or fills a 2 dimensional
 array after the use of \code{tapply}which eliminates dimension values if
 there is no content for them}\cr
 \item{\code{fill.1dimension}}{check values and/or fills a 1 dimensional
 array after the use \code{tapply} which eliminates dimension values if
 there is no content for them}\cr
 \item{\code{fill.distclases}}{check values and/or fills a 1 dimensional
 array after the use \code{tapply} which eliminates dimension values if
 there is no content for them.  This is specific for neighborhood 
 calculations.}\cr
 \item{\code{unidentified.species}}{creates a vector that indicates whether
 a tree belongs to an unidentifed species.  \code{unidenstr} contains the
 possible codes for unidentified species at the sites.}\cr
 \item{\code{unwind.matrix}}{returns a vector that
 is the contents of a matrix "unwound" column by column}\cr
 \item{\code{rndown5}}{rounds by taking the floor of dbh / 5}\cr
 \item{\code{xydist}}{distance between x,y coordinates}\cr
 \item{\code{xydistvect}}{distance between vector of x,y coordinates, x
 as first element, y as second.  Used with \code{apply}}\cr
 \item{\code{setsize}}{sets point sizes for plotting location of trees in \code{map}}\cr
 \item{\code{quadslope}}{calculates the slope of a quadrate and called from 
 \code{meantopo.quad}}\cr
 \item{\code{calcslope}}{calculates the slope of a quadrate and called from 
 \code{quadslope}}\cr
}
 
 
% \item{\code{select.dbhrange}}{returns rows of a dataframe that have a dbh 
% value greater than or equal to the dbhmin and less than dbhmax}
 % \item{\code{showci}}{draws a line between two confidence intervals and 
% assigns an appropriate colour to that line}
}
\note{
For developers:  The CTFS internal utility functions are found in a file approximately named:
\code{utilties.versionnumber.R}. Some may be located in other files 
because they are used with those functions only.
}

\author{Rick Condit, Suzanne Lao, Pamela Hall}
\seealso{\code{\link{CTFS.quadfunctions}},\code{\link{CTFS.utilities}},
\code{\link{CTFS.readwrite}}, \code{\link{CTFS.groupcats}}
}

\keyword{internal}

