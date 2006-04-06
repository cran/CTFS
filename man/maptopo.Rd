\name{maptopo}
\alias{maptopo}
\concept{maps}
\concept{topography}
\concept{contour}
\concept{spatial pattern}
\concept{tree location}

\title{Plot a Contour Map of Topography of a CTFS Plot}
\description{
  Draws a contour plot (default to 5 m intervals) for a CTFS Plot.  This 
	function can be run independently or in order to add a contour plot to 
	a plot of the tree location of a species.  The defaults provide a 
	\dQuote{nice} plot.
}
\usage{
maptopo(elevmat, plotdim = c(1000, 500), add = FALSE, export = "no", 
	interval = 0, maintitle = "Topographic Contours", 
	ht = 6, wd = 9, plotside = 4.5, filename = "maptopo")
}

\arguments{
  \item{elevmat}{\strong{matrix} of elevation data, preferably at 5 x 5 m 
	scale, usually the matrix component of the list \code{\link{bcitopo.info}}}
  \item{plotdim}{vector of 2 numbers, maximum x, y dimensions of plot in m}
  \item{add}{logical, TRUE = add contour plot to tree plot, only for use 
	when NOT run independently, FALSE = for making contour plot only}
  \item{export}{character, indicating whether to plot to window or save as 
	file.  Possible values are: \dQuote{\code{no}}, \dQuote{\code{bmp}}, 
	\dQuote{\code{jpg}}, \dQuote{\code{emf}}, \dQuote{\code{pdf}}}.
  \item{interval}{contour line intervals}
  \item{maintitle}{character string, main title for plot}
  \item{ht}{height of graphics window, in inches}
  \item{wd}{width of graphics window, in inches}
  \item{plotside}{size of plot within the graphics window }
  \item{filename}{default filename for export}}
}
\details{
  This function can be run independently to just produce a contour plot 
	or can be called from \code{\link{map}} to add a contour plot to a plot 
	of tree location for a species.
	
	Make sure that \code{elevmat} is a matrix.  In most cases this would be 
	the matrix component of the list \code{\link{bcitopo.info}}.
}
\value{
  If \code{export=}\dQuote{\code{no}} then a map is plotted in a window and all 
	information requested in plotted in the window.  If the user choses to 
	export the file, then the same plot on the screen will be saved in a 
	file in the working directory in the format of the users choice.  See 
	above for file name creation.
}
\author{Rick Condit and Pamela Hall}
\seealso{\code{\link{map}}, \code{\link{bcitopo.info}}, \code{\link{contour}}
}
\examples{\dontrun{
data(bcitopo.info)

# minimum needed
maptopo(bcitopo.info$mat)

# for export just add:
maptopo(bcitopo.info$mat,export="pdf")
}}
\keyword{hplot}
