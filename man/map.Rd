\name{map}
\alias{map}
\alias{map1species}
\concept{map}
\concept{topography}
\concept{tree location}
\concept{spatial pattern}

\title{Plot of the Spatial Location of Trees of a Species}
\description{
  This function provides nicely proportioned maps of species distributions 
	within a plot and topographic contour lines if requested.  More than one 
	species can be plotted and will have different coloured symbols for ease 
	of interpretation.
}
\usage{
map(census.spp, spname, elevmat = NULL, plotdim = c(1000, 500), spp.info = NULL, 
	maintitle = NULL, titlepos = "default", addlegend = TRUE, legpos = "default", 
	clrs = "default", symb = "default", size = rep(-1, 3), dbhclass = c(10, 100, 300, 3000), 
	ht = 6, wd = 9, plotside = 4.5, export = "no", bw = FALSE, topoint = 0)
	
map1species(sppdata,plotdim=c(1000,500),sppclrs="black",symbol=16,size=rep(-1,3),
    dbhclass=c(10,100,300,3000),plotside=3,add=FALSE)

}
\arguments{
  \item{census.spp}{data file as a \strong{list} with species as names}
	\item{sppdata}{data file as a \strong{dataframe} of a single species}
  \item{spname}{a species name code in quotes or a vector of names}
  \item{elevmat}{\strong{matrix} of elevation data, preferably at 5 x 5 m scale, 
	see	\code{\link{bcitopo.info}}}
  \item{plotdim}{vector of 2 numbers, maximum x, y dimensions of plot in m}
  \item{spp.info}{name of \strong{dataframe} that contains each
  species and information about it, used to print Latin species name on 
	plot for legend, see \code{\link{bcispp.info}}}
  \item{maintitle}{character string, main title for plot}
  \item{titlepos}{a vector of 2 numbers, position of title on plot in relation to \code{plotdim}}
  \item{addlegend}{logical, TRUE = add legend with color for each species}
  \item{legpos}{position of legend on plot in relation to \code{plotdim}}
  \item{clrs}{colors for symbols when more than one species is plotted.  
	the defaults (8) are: 
	\dQuote{\code{blue}}, \dQuote{\code{red}}, \dQuote{\code{black}}, 
	\dQuote{\code{green}}, \dQuote{\code{yellow}}, 
	\dQuote{\code{lightblue}}, \dQuote{\code{orange}}, \dQuote{\code{pink}}.
	If \code{bw} is TRUE, then the colors 
	range from black to light grey}
	\item{sppclrs}{a single color for mapping 1 species, use for 
	\code{map1species} only}
  \item{symb}{symbol shapes when more than one species is plotted, default 
	is a circle, 16}
	\item{symbol}{default for mapping 1 species, a circle}
  \item{size}{symbol size based on tree dbh and based on number of 
	individuals of most abundant species when more than one species is 
	plotted, see \code{dbhclass} for size classes}
  \item{dbhclass}{range of dbh classes for symbol sizes}
  \item{ht}{height of graphics window, in inches}
  \item{wd}{width of graphics window, in inches}
  \item{plotside}{size of plot within the graphics window }
  \item{export}{character, indicating whether to plot to window or save as 
	file.  Possible values are: \dQuote{\code{no}}, \dQuote{\code{bmp}}, 
	\dQuote{\code{jpg}}, \dQuote{\code{emf}}, \dQuote{\code{pdf}}}.
  \item{bw}{logical, TRUE = use range of grey for symbols, FALSE = use 
	colours}
  \item{topoint}{interval for contour lines}
	\item{add}{logical, must be set to FALSE for plotting 1 species when using 
	\code{map1species}. This variable is set to TRUE by \code{map} which 
	calls \code{map1species} for each species.}
}
\details{
	This function provides nicely proportioned maps of species distributions
	within a plot and topographic contour lines if requested.  The dataset
	input can be the entire site tree data (structured as a list of species)
	or a list containing only the data frames of the species that are to be
	plotted.  The full Latin name can be easily added by specifying the file
	containing them.  The topographic data file MUST be in matrix format.  In
	most cases this would be the matrix component of the list
	\code{\link{bcitopo.info}}.  Files of larger scale can also be used.
	
	The user can alter most of the options that affect the size of the
	window, the symbols, their colours, their size (related to dbh of the
	tree), titles and position of titles and legends.  Take care in
	manipulated these options as some combinations provide unexpected
	results.  The same parameter settings may produce quite different looking
	plots on different machine OSs.  The defaults work well on MSWindows (of
	more recent flavor) and MacOSX.
	
	Any plot that can be created can also be exported as a graphics file.
	Check the available options.  For Windows, the preferred export format is
	winmetafile (\code{emf}) and for MaxOSX its \code{pdf}.  The files are
	written to the current working directory.  File names are created from
	the species name.  If many species are plotted then the file name is
	\code{multisppmap}.  The file extension is appropriate for the type of
	graphics file the user chooses.
}
\value{
  If \code{export=}\dQuote{\code{no}} then a map and all 
	information requested is plotted in the window.  If the user choses to 
	export the file, then the same plot on the screen will be saved in a 
	file in the working directory in the format of the users choice.  See 
	above for file name creation.
}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{maptopo}}, \code{\link{bcitopo.info}}, 
\code{\link{tst.bci90.spp}}, \code{\link{bcispp.info}}, \code{\link{contour}}
}
\examples{\dontrun{
data(tst.bci90.spp)
data(bcispp.info)
data(bcitopo.info)

# minimum needed
map(tst.bci90.spp,"socrex")

# much fancier plot
map(tst.bci90.spp,names(tst.bci90.spp),elevmat=bcitopo.data$mat,spp.info=bcispp.info, maintitle="BCI test data set")

# for different colour scheme, see colors in R help for more colors
map(tst.bci90.spp,spname=names(tst.bci90.spp),clr=rainbow(length(names(tst.bci90.spp))))

# for export just add:
map(tst.bci90.spp,names(tst.bci90.spp),elevmat=bcitopo.data$mat,spp.info=bcispp.info, maintitle="BCI test data set",export="pdf")
}}
\keyword{hplot}
