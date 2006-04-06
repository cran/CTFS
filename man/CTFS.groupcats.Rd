\name{CTFS.groupcats}
\alias{CTFS.groupcats}
\concept{groups}
\concept{categories}
\concept{classes}
\concept{split}
\concept{separate}
\concept{restructure data}

\title{Overview of Functions that Create Categories for Grouping Trees}

\description{
 Provides an overview of a variety of functions that create vectors which
 contain values for each tree that indicate which "group" it belongs to.
 For instance, dbh class, species name, habitat type, quadrate, growth
 form.}

\details{

	FUNCTIONS FOR CREATING CATEGORIES
	
	There are a variety of ways of creating a vector of categories for
	individual trees which can be used for analysis of dynamics by
	those categories.  The functions described here create vectors of
	category values for each tree based on characteristics of the tree
	(species, dbh size, growth form, etc.)  or a location in the plot
	(quadrate, habitat).  For instance, one can compute the mortality
	of a species by habitat by making 2 vectors, one for species and
	one for habitat.
	
	\item{\code{\link{sep.quadinfo}}}{Creates a Vector of Quadrate
	Characteristics for each Tree in a Dataset.}
	\item{\code{\link{sep.sppinfo}}}{Creates a Vector of Species
	Characteristics for each Tree in a Dataset.}
	\item{\code{\link{sep.dbh}}}{Creates a Vector of DBH Classes for
	each Tree in a Dataset.}
	
	MAKING VECTORS FOR GROUPING TREES FOR ANALYSIS
	
	These examples use a test dataset \code{tst.bci90.full}
	\preformatted{
	dim(tst.bci90.full)
	[1] 12421    14
	
	unique(tst.bci90.full$sp)
	[1] "alsebl" "psycde" "socrex"
	}	
	1.  Vector of species names same length as datafile
	\preformatted{
	sppname.vct <- tst.bci90.full$sp
	
	length(sppname.vct)
	[1] 12421
	
	sppname.vct[1:5]
	[1] "alsebl" "alsebl" "alsebl" "alsebl" "alsebl"
	}
	
	NOTE: computing growth for each species by dbh classes can be done 
	using \code{\link{growth.eachspp}}
	
	2.  Vector of dbh classes same length as datafile. \code{cut} returns a vector of factors and \code{Levels} 
		are provided
	\preformatted{
	dbh.vct <- sep.dbh(tst.bci90.full)
	
	dbh.vct[1:5]
	[1]  <NA>      300.10000 100.300   100.300   100.300  
	Levels: 10.100 100.300 300.10000

	dbh.vct <- sep.dbh(tst.bci90.full,dbhclass=c(10,20,30,40,50))
	
	dbh.vct[1:5]
	[1] <NA>     50.10000 50.10000 50.10000 50.10000
	Levels: 10.20 20.30 30.40 40.50 50.10000
	
	table(dbh.vct)
	dbh.vct
   	10.20    20.30    30.40    40.50 50.10000 
    4115     1512      753      436     2337 

	}

	NOTE: computing growth for each dbh class with \code{growth.dbh}.
	
	3. Vector of habitat location of trees
	\preformatted{
	habitat.vct <- sep.quadinfo(tst.bci90.full,bciquad.info,by.col="hab")
	
	length(habitat.vct)
	[1] 12421
	
	habitat.vct[1:10]
	[1] NA  5  6  7  7  7  5  5  6  6
	 }

	4. Vector of growth form of species
	\preformatted{
	grform.vct <- sep.grform(tst.bci90.full,bcispp.info,by.col="grform")

	length(grform.vct)
	[1] 12421

	table(grform.vct)
	grform.vct
    	M     S     T 
 	1133   160 11128 
	
	length(grform.vct)
	[1] 12421
	}
	
}
\note{
For developers:  the functions described here are found in
\dQuote{readandwrite.R}
}
\author{Pamela Hall}
\seealso{
	\code{\link{CTFS.abundance}}, \code{\link{CTFS.basalarea}},
	\code{\link{CTFS.diversity}}, \code{\link{CTFS.growth}},
	\code{\link{CTFS.mortality}}, \code{\link{CTFS.recruitment}},
	\code{\link{CTFS.quadfunctions}},
	\code{\link{CTFS.readwrite}}, \code{\link{CTFS.utilities}}
}
\keyword{univar}
