\name{mergecensus}
\alias{mergecensus}
\concept{creating datasets}
\concept{merging censuses}
\concept{status}
\concept{create status}

\title{Merging Two Censuses of CTFS Plot Data}
\description{ 

	THIS FUNCTION IS OBSOLETE - USE \code{create.status} TO ASSIGN THE
	\code{status} VARIABLE TO A CTFS DATAFILE.
	
	Combines the census data of individuals trees from two
	censuses into a single record so that measures of dynamics can be
	computed.  
	}
\usage{
mergecensus(census0, census1, pasoh = FALSE)
	}

\arguments{
  \item{census0}{name of census datafile for first census, must be a
  \strong{dataframe}, must be of same length as \code{census2}}
  \item{census1}{name of census datafile for second census, must be a
  \strong{dataframe}, must be of same length as \code{census1}}
	\item{pasoh}{logical value for indication of site, some \code{status} 
	values at pasoh are treated differently than all other sites.  Use pasoh 
	= TRUE if census is from the Pasoh site.}
	}
\details{
	The variable \code{status} is created in this function.  \code{status} 
	provides information on the condition of a tree in a given census.  
	The current values used are:
	"A" alive, dbh >= 10, pom and date are valid values
	"AS" alive, 0 < dbh < 10 which is too small for growth
	measurements, pom and date are valid values
	"AB" alive, dbh = 0, broken stem with no pom to measure
	"D" dead, dbh and pom are not valid values, both set to NA
	"P" recruit, not yet large enough for inclusion in census, dbh and 
	pom are not valid values, both set to NA
	NA other conditions, missing data, inconsistent data, dbh and pom 
	are not valid values, both set to NA\cr
	
	At Pasoh, there were some differences in the way dbh were coded for 
	trees that were alive but 0 < dbh < 10.  The status of these trees 
	is:
	"A" alive, 0 < dbh < 10, but too small for growth measurements\cr
	Both of the input datasets must already be structured to include the 
	variables and col names referred to in this function.\cr
	The censuses do not have to be sequential, any 2 censuses can be 
	used.  Just make sure to name the generated file appropriately so you 
	know which census data is included. EG: \code{bci9095.full} or 
	\code{bci8590.full}\cr
	
	The \code{site} parameter is only needed, at this time, to
	distinguish "Pasoh" from the other sites as the measurements for
	trees less than 50 mm dbh during the initial census were handled
	differently than at other sites and this affects the computations
	of growth.\cr
	}

\value{
	returns a data frame of the same number of rows as the input files 
	and with new variables to distinguish the censuses from each other
	}
\note{ 
	At this time, \code{mergecensus} only handles the true first 
	census as a value for \code{census0}.  For all other sites, any 2 
	censuses can be used.
	}
\author{Suzanne Lao, Rick Condit and Pamela Hall}
\seealso{\code{\link{CTFS.readwrite}}}

\examples{\dontrun{

mergecensus(site="bci",census0=tst.bci90.full,census1=tst.bci95.full)->tst.bci9095.full
dim(tst.bci9095.full)
[1] 12421    14

names(tst.bci9095.full)
 [1] "tag"     "sp"      "gx"      "gy"      "dbh0"    "dbh1"    "pom0"    "pom1"    "date0"   "date1"   "codes0"  "codes1" 
[13] "status0" "status1"

tst.bci9095.full[1:5,]
     tag     sp    gx    gy dbh0 dbh1 pom0 pom1 date0 date1 codes0 codes1 status0 status1
1 -27784 alsebl  -9.0  -9.0   NA   NA    0    0     0     0      *      *    <NA>    <NA>
2     47 alsebl 984.3 341.6  437  426    2    3  3702  5382      B      B       A       A
3     49 alsebl 985.3 328.9  228  228    1    1  3632  5396      *      *       A       A
4     68 alsebl 985.7 275.8  278  277    1    1  3627  5390      *      *       A       A
5     71 alsebl 999.9 277.6  269  318    1    1  3627  5390      M      *       A       A
}
}
\keyword{manip}
