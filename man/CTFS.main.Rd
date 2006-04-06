\name{CTFS.main}
\alias{CTFS.main}
\concept{CTFS}
\concept{Center for Tropical Forest Science}
\concept{overview}
\concept{front page}
\concept{first page}
\concept{list of functions}
\concept{manual}

\title{Center for Tropical Forest Science (CTFS) : Analysis of Long Term Forest 
Plot of Large Size}
\description{
	The CTFS Long Term Forest Plots provide ecologists and foresters a unique
	opportunity to describe and interpret the demographics of individual tree
	species, the spatial patterns of species population dynamics and to
	compare these processes across tropical forests throughout the world.
	The CTFS package provides user friendly access to R based tools for data
	management and analysis.  These tools represent a continuing world-wide
	collaboration of scientists interested in unlocking the mysteries of
	these diverse and dynamic forests.
}

\details{	
	The CTFS project web site is: \url{http://www.ctfs.si.edu/}.
	
	The CTFS package provides functions for analysing tropical forest
	structure, dynamics and species composition including the effects 
	of habitat differences, neighborhood effects.  New analyses will
	be continually added as collaborators produce and publish their
	results.

	The CTFS package is explictily designed to be applied to the data
	collected from the CTFS Long Term Forest Plots. It is also intended to be
	sufficiently general in data requirements and broad in types of
	analyses to be of use to other researchers in the field of forest
	dynamics.
	
	The CTFS plots are a complete enumeration and mapping of all stems
	greater than 1 mm in diameter of all woody species in a many hectare
	plot.  The plot size varies from 16 to 50 ha across sites with a 50 ha
	plot established at most sites.  These plots are reenumerated on average
	every 5 years.  CTFS and collaborators have established over 20 plots
	through out the tropics which includes ever-wet rain forests and
	seasonally dry forests.  To date, over 3 million trees have been
	enumerated, identified, mapped and, at many sites, recensused many times.
	Recensuses and the establishment of new sites continues.
	\url{http://www.ctfs.si.edu/sites/programs/programs.htm}
	
	For a complete description of the method of plot establishment and
	enumeration see	\url{http://ctfs.si.edu/bibliography}.  Look for Condit, R. 1998.
	Tropical Forest Census Plots.  Springer-Verlag, Berlin, and R. G.
	Landes Company, Georgetown, Texas.
	
	Topics of Analysis (Examples of what is and will be available)
	\enumerate{
		\item Forest Structure: size distribution, abundance, basal area, 
	biomass, spatial distribution of trees. Available version 1.0.
		\item Population Dynamics: population abundance change, mortality, recruitment, 
	diameter growth rates.  Available version 1.0.
		\item Density Dependence: analysis of effect of immediate biotic 
	neighborhood on structure and dynamics of populations.  Available version 1.0.
		\item Above Ground Biomass estimations: can be analysed as a characteristic
	of an individual tree, changes over time and space.  Available version
	1.0.
			\item Habitat Specificity: based on abundance and dynamics. 
	Planned for next update.
	}
	
	Manual chapters are also available that focus on the use of R and the 
	CTFS package for the CTFS plots and data structure.  	
	
	R is a language and environment for statistical computing and
	graphics.  R is distributed as free software and compiles and runs
	on a wide variety of UNIX, Windows and Mac operating systems.
	CTFS uses R to analyze datafiles of its large Forest Dynamics Plot
	Network.  More information about R can be found at the R Project
	for Statistical Computing Website:\url{http://www.r-project.org}.
}

\author{package maintainer: Pamela Hall, \email{phall@alum.mit.edu}}
\references{
\url{http://www.ctfs.si.edu}
}
\seealso{
	\code{\link{CTFS.abundance}}, 
	\code{\link{CTFS.basalarea}},
	\code{\link{CTFS.biomass}},
	\code{\link{CTFS.datafiles}},
	\code{\link{CTFS.diversity}},
	\code{\link{CTFS.groupcats}},
	\code{\link{CTFS.growth}},
	\code{\link{CTFS.maps}},
	\code{\link{CTFS.manual}},
	\code{\link{CTFS.mortality}},
	\code{\link{CTFS.quadfunctions}}, 
	\code{\link{CTFS.recruitment}}
%	\code{\link{CTFS.quick.reference}}	
}

\keyword{univar}

