\name{CTFS.datafiles}
\alias{CTFS.datafiles}
\concept{data structure}
\concept{column names}
\concept{variables}

\title{Overview of Data File Structure for CTFS Package}

\description{
	Provides an overview of the structure of datafiles that is
	required for use with the CTFS Package.  This includes the
	standard format for all CTFS datafiles and discussion of what
	aspects of the structure are flexible for use with non-CTFS
	datasets.
	}

\section{List of Data Files}{	
	
	The general format of census data files names is:
	
	\strong{siteyear.type}

	\code{\link{tst.bci90.full}} Tree based dataset (dataframe) containing
	information from a single census (1990).
	
	\code{\link{tst.bci90.spp}} Tree based dataset (list) containing informatin
	from a single census (1990) composed of dataframes of all trees within a 
	single species.
	
	\code{\link{tst.bci90.quad}} Tree based dataset (list) containing
	information from a single census (1990) composed of dataframes of all
	trees within a single quadrate.
	
	\code{\link{tst.bci90.mult}} Stem based dataset (dataframe) containing the
	dbh measurement of all secondary stems from trees with multiple
	stems.  Multiple stems of each main stem are provided sequentially which 
	means a tree tag number can occur many times in the data file.
	
	Datasets that are based on information that is valid for all
	censuses do not have a \dQuote{year} designation.  The general
	format for these data file names is:
	
	\strong{sitetype.info}
	
	\code{\link{bcispp.info}} Species based dataset (dataframe) containing
	information for each species ever recorded in plot, eg. growth
	form, breeding system, full latin name, etc.	
	
	\code{\link{bciquad.info}} Quadrate based dataset (dataframe) containing
	information on each quadrate (20 x 20 m), eg. elevation, habitat
	designation, etc.
	
	\code{\link{bcielev.info}} Dataset (dataframe) containing
	elevation measures for x, y coordinates.  At BCI this is on a 5 x 
	5 m scale.
}

\section{Information about Data Files}{
	
	Each census data contains exactly the same number of records, one 
	for each tree ever recorded in the plot.  As new enumerations are 
	completed and trees are added to the census due to recruitment,
	these records are added to previous census files as well as the
	file for the year they were recruited into the plot.  Trees that
	die, break or have been lost track of (by error) are retained in
	all census files.  Variable values are given as appropriate for
	such trees (see details on structure of each type of data file).
	
	\code{bcispp.info} contains a record for every species ever
	recorded in the plot.  As species names are changed due to
	taxonomic updates, the names are changed in this datafile and if
	this requires a change in species codes in the census data, then
	the changes are done for all census files.  New information about 
	species can be added to this file, even if not available for every
	species.  New variables can be used for new analyses.
	
	\code{bciquad.info} contains a record for each 20 x 20 m quadrate 
	in a plot.  The total number varies according to the plot shape.  
	The basic information is elevation, slope and convexity but other 
	information can be added.  Retaining the value for \code{habitat} 
	in this file makes analyses using this variable very easy. 
	Multiple versions of a \code{habitat} variable can all be stored
	here.
	
	\code{bcielev.info } contains the x, y coordinates of each point
	that was surveyed for creating a topographic map of the plot. 
	This was done to different degrees of accuracy in each plot.  This
	information is used for drawing topographic maps and for providing
	values in \code{bciquad.info }.
	
	
	To see the details on the contents of each file and its variable
	values, lengths, etc., see the help pages for the named file or
	click on the links above.
}
\references{\url{http://www.ctfs.si.edu}}
\seealso{\code{\link{CTFS.main}}}
\keyword{univar}
