\name{CTFS.biomass}
\alias{CTFS.biomass}
\concept{biomass}
\concept{growth}
\concept{dynamics}
\concept{wood specific gravity}
\concept{wood density}

\title{Overview of Analysis of the Biomass Dynamics of Tropical Tree Populations}

\description{
 Provides an overview for the analysis of biomass and change in biomass for
 trees.  This overview includes functions that are available, supporting
 functions and ways for using the options of the main functions.
 
 Above ground biomass (AGB) is computed on the tagged tree from a census
 file and from the sum of any multiple stems of the tagged tree that exist.
 Biomass is computed using the specific wood gravity, where available, for
 each species.  The census and multi-stems datasets are used for these
 analyses.  If the multi-stem datasets do not exist, biomass is computed
 only from the main dbh of the tagged tree.
}

\details{

	FUNCTIONS USED TO COMPUTE BIOMASS
	
	\item{\code{\link{biomass}}}{Main function that controls the
	options for computation and calls all the remaining described here.}
	
	\item{\code{\link{pick.dbh}}}{Determines which dbh to use when gorowth is
	not good.}
	
	\item{\code{\link{getwsg.vct}}}{Create a vector of wood
	specific gravities for given species.}
	
	\item{\code{\link{agb.ind}}}{Computes above ground biomass (AGB) for each
	individual which is the sum of the AGB for all recorded stems of the individual.}
	
	\item{\code{\link{agb.mult}}}{Computes above ground biomass for 
	multiple stems which is then added to \code{\link{agb.ind}}.}
	
	\item{\code{\link{gr.raw}}}{Computes growth of individual main stem for
	the purpose of determining if a recorded dbh is reliable.}
	
	\item{\code{\link{badgrowth.tree}}}{Creates a file of trees for which a dbh
	was estimated because the growth rate was outside provided limits, tree 
	was too large to be a recruit or DBH values were missing when the 
	\code{status} equals \dQuote{A} indicating the tree should have been 
	measured.}
	}

\section{Computation of Biomass}{
	The biomass for individual trees is computed by \code{\link{biomass}}
	using the DBH and wood specific gravity (wsg) in the form of models from 
	Table 4	as published in \cite{Chave, J. et.al.(2005). Tree allometry 
	and improved estimation of carbon stocks and balance in tropical 
	forests. Oecologia. 145. pg 87}.
	
	\preformatted{biomass(census1, census2, multfile1 = NULL, multfile2 = NULL, wsgfile = NULL, 
	forest = NULL, badgrowth = "meangrowth", dbhuse = "second", 
	trim = c(-5, 75), recsize = 75, rounddown = FALSE)}

	The quality of the estimated biomass is strongly dependent upon the
	reliability of the measured DBH. Evaluating the reliability of the
	measured DBH is done in two ways.  For trees that are listed as alive in
	both censuses without broken stems or changes in point of measurement,
	the DBH increment (growth rate) between two censuses for each tree is
	computed.  If the growth rate for an individual tree is outside the range
	provided by \option{trim}, then the unreliable DBH is estimated for that
	census using the method determined by the values of \option{badgrowth,
	dbhuse, trim} (see below).  A second form of unreliable DBH is found in
	trees that are very large when first measured as recruits which suggests
	that they may in fact have been missed in a previous census.  The
	parameter \option{recsize} can be specfied by the user.  Only trees that
	are recruits to the second census in the computation are checked against
	this parameter value.  If they are larger or equal to the value of
	\option{recsize} then their DBH in the first census is estimated using
	the mean growth rate of their DBH class and their status is changed to
	\dQuote{A} for the first census.  These trees are therefore no longer
	considered recruits.  \option{recsize} is arbitrarily set equal to 110 mm
	which means a tree would have had to grow at least from 10 to 110 mm in
	the typical 5 year intercensus interval or 20 mm per year.  At BCI, 89\%
	of the species had average annual growth rate less than this between the
	1990 and 1995 census.
	
	Missing DBHs for trees with \code{status} values of \dQuote{A} as also 
	estimated by either adding or substracting the mean DBH increment for 
	their DBH class depending on which census the DBH is missing.	
	
	After any estimates of DBH are completed, the AGB is computed for the
	main stem and any multiple stems as provided by \option{multfile1,
	multfile2} using the wood specific gravity as provided by
	\option{wsgfile}.  The above ground biomass (AGB) of the main stem and
	the multiple stems for a each tree are summed resulting in the total AGB
	for that tree.  The function \code{biomass} returns of \strong{dataframe}
	of the same length as the input census files (\code{census1, census2})
	with the tree \code{tag} and the values of AGB for the main stem, the sum
	of the multiple stems, the total of for the tree for each census.\cr
	
	The usual mode is to provide 2 censuses (see below for computing biomass
	when only one census is available).  \cr
	
	\describe{
	\item{\bold{Diagnosing Unreliable AGB values}}{The user should first run 
	\code{\link{badgrowth.tree}} to obtain a dataframe of trees with 
	unreliable DBHs.  Details of how these are computes is given below.  
	The DBH measurements of these trees should be inspected carefully to 
	determine if there is an error in the census data file or if the DBH is 
	simply missing or unreliable.  Actual errors in the census data should 
	be corrected IN THE CENSUS DATA FILE.  Then the biomass can be computed 
	using this corrected data set and trees with unreliable, not simply 
	improperly recorded, DBHs can be estimated for the purpose of AGB 
	computations.  The user can rerun \code{badgrowth.tree} to retain a list 
	of trees for which DBH was estimated.	 The \code{errtype} provides the 
	justification for DBH estimation.  Remember to record the other 
	parameter values so that the form of DBH estimation is known.}\cr
	%
	\item{\bold{Computing Growth}}{An unreliable DBH can be estimated only if 
	a growth rate can be computed.  If the point of measurement (POM) 
	of a tree has changed between censuses, then the growth rate cannot 
	be computed and no DBH estimation is done.  If one of the DBH is 
	not available (recruit, mortality, missing information), the growth 
	rate cannot be computed and no DBH estimation is done.}\cr
	%
	\item{\bold{Growth Outside Acceptable Range}}{For high growth rates, any
	growth rate that exceeds the maximum provided by \option{trim} is
	considered unreliable.  The function \code{\link{trim.growth}} is used to
	identify unreliable low growth rates (usually large negative values).
	The determination of an unreliable growth is based on a linear model
	estimating the standard deviation of DBH measurements (due to error based
	on repeated measures of the same trees at BCI); the parameters slope and
	intercept define the linear relationship between this error deviation and
	DBH. Any case where the second DBH measure is more than 4 standard
	deviations below the first the growth rate is considered too low and one
	of the DBHs is unreliable.  The default values of slope and intercept are
	based on DBH remeasure tests done in both 1995 and 2000 at BCI. A line
	was fitted through the absolute DBH errors as a function of DBH in both
	years; the average slope and intercept is used here.  The overall affect
	is that reliable growth rates are allowed to be more negative for larger
	trees than for smaller trees.  The minimum value set in \code{trim} is
	used for estimating a DBH that has been determined to be unreliable using
	one of the options for \code{badgrowth} (see below).}\cr
	%
	\item{\bold{Determining the Reliable DBH}}{The value of \option{dbhuse} sets the
	reliable DBH from which the unreliable DBH will be computed.  Values can
	be \code{first} or \code{second} which refer to the census
	files being used.}\cr
	%
	\item{\bold{Form of Estimation of Unreliable DBH}}{There are four forms
	of estimation which are controlled by \option{badgrowth}: }
	
			\describe{
				\item{\option{keep}}{Do not estimate DBH, retain DBH exactly as
					provided in the census files and proceed with AGB based on
					this DBH.}
							
					\item{\option{nogrowth}}{Estimate unreliable DBH by setting it equal
					to the reliable DBH.  This results in the tree having no
					growth.}
							
					\item{\option{meangrowth}}{Estimate the unreliable DBH by
					substracting (or adding) the mean growth of trees of similar
					DBH from/to the reliable DBH.  The mean growth of DBH classes
					is determined internally in the \code{biomass} function 
					using the trees provided by the census files.
					The DBH classes used are finely divided for small trees and
					more grossly so for larger trees.  The DBH classes are:
					10, 20, 50, 100, 200, 300, 400, 500, 600, 700+).}
							
					\item{\option{trim}}{Estimate the unreliable DBH by substracting (or
					adding) the maximum or minimum value of \code{trim}, as
					appropriate, to the reliable DBH.  The maximum value is in mm per 
					year.  The minimum value is the number of standard deviations   
					See \code{\link{trim.growth}} for 
					complete details.}
	}}}}

	\section{Examples of Parameter Settings}{
	
	Example 1: \preformatted{trim(-5,75), dbhuse=second, badgrowth=meangrowth}
	
	Trees with a growth rate of > 75 mm per year and for which the second dbh
	measure is more than 4 standard deviations below the first will have the
	first dbh estimated based on the mean growth rate of the trees in the dbh
	class determined by the second dbh.  The mean growth rate for that dbh
	class will be subtracted from the second dbh to compute the first dbh.
	This estimated dbh will be used to compute all values of AGB for the
	first census.  \bold{This is usually the preferred method of estimation with 
	some variation among sites in the acceptable maximum growth rate set by} 
	\code{trim}.\cr\cr
	
	Example 2: \preformatted{trim(-5,75), dbhuse=first, badgrowth=trim}
	
	Trees with a growth rate of > 75 mm per year and for which the second dbh
	measure is more than 4 standard deviations below the first will have the
	second dbh estimated based on the minimum and maximum growth values in 
	\code{trim} using the first dbh as the reliable dbh.  The trees with the 
	second dbh less than the first will have the second dbh computed as the 
	first dbh plus (-5).  The trees with their second dbh greater than the 
	first will have the second dbh computed as the first dbh pluse (75).  
	This way growth rates will not exceed the values set by \code{trim}.\cr
}

	\section{Diagnostics for DBH and Biomass Estimation}{
		\describe{
			\item{\bold{Listing estimated DBHs}}{Trees for dbhs have been
			estimated can be idenfitied by using \code{\link{badgrowth.tree}}.
			Set the parameters of this function exactly as used for
			\code{\link{biomass}}.  The trees that had a dbh estimated in the
			manner determined by the parameter settings will be listed in a
			\strong{dataframe}.  The initial dbhs and growth rate and the
			estimated dbhs and growth rate will be provided, as well as all the
			other values from the \code{census1, census2} data files.}\cr
			
			\item{\bold{Visualizing computed AGB}}{The AGB for quadrates or
			hectares can be computed and visualized using \code{image}.  See 
			example below for visualization the AGB (in tons) for quadrates.}
			}}

\author{Pamela Hall and Jerome Chave}
\references{\cite{Chave, J. et.al.(2005). Tree allometry 
	and improved estimation of carbon stocks and balance in tropical 
	forests. Oecologia. 145. pg 87}}
\seealso{\code{\link{CTFS.growth}}}
\examples{\dontrun{
plotdim=c(1000,500)
gridsize=20
bci.agb <- biomass(bci90.full,bci95.full,mult1=bci90.mult,mult2=bci95.mult,wsgfile=wsg.ctfs,forest="moist")
bciquad.vct <- gxgy.to.index(bci90.full$gx,bci90.full$gy)
bciquad.agbt1 <- tapply(bci.agb$AGBt1,bciquad.vct,sum,na.rm=TRUE)/1000
x=seq(gridsize/2,plotdim[1]-gridsize/2,gridsize)
y=seq(gridsize/2,plotdim[2]-gridsize/2,gridsize)
agbt1.mat=matrix(bciquad.agbt1,nrow=length(x),ncol=length(y),byrow=TRUE)
image(x,y,z=agbt1.mat,col=topo.colors(6),xlab='AGB 1990 (tons)',xlim=c(0,plotdim[1]),ylim=c(0,plotdim[2]))
}}


\keyword{univar}



