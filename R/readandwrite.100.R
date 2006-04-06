# readandwrite.092.R

# functions to read/write, re-structure and create useful vectors from 
# CTFS data files and supporting files

# sep.data() convert data frame to list by a variable in the data frame
#     new version of splitdata()
# sep.vct.data() convert data frame to list by a vector that defines list components
#     like splitdata() but can use a variable not in the initial data frame
# splitdata() original version is split.data() - without the get() function.  To be phased
# 		out because function name does not meet R specs,
# sep.quadinfo()  creates a vector of variable values derived from sitequad.info
# sep.sppinfo()  creates a vector of variable values derived from sitespp.info
# sep.dbh()  creates a vector of variable values derived from the dbh of a tree
# create.status() new version by Suzanne and Pam
# mergecensus() new version by Suzanne and Pam
# texttordata()  read text file into an RData format
# sep.species() returns a subset of the split species list, either a dataframe or a list of dataframes
# load.species() returns a single dataframe from the split species list
# readelevdata() read text or RData dataframe into a matrix for topographic maps
# elev.to.list() read dataframe of x,y,elevation and convert of a list with the original data frame and a matrix configuration of the elevation data

# assemble.demography() - assembles array returned from demographic functions into a dataframe
# assemble.table.each() - ERROR

###########################
# VERSION INFORMATION
# created by PHall on 4/14/2005
# similar function by Rick Condit

# sep.data()
#
# Separates a data.frame into a list of data.frames by
# the value of variable sepcol which is the column name of the variable by which
# the census1 will be separated.

# The output is a list of data.frames.  Each data.frame is composed of the records of the trees 
# that match a single value of the sepcol variable. For example, if separating on
# species, then there is one data.frame for each species. Other possibilities:
# split by quadrat number (the program create.quadratlistdata in densdepend.r does
# this), growth form, dbh category, habitat.

# Datafile name must be in quotes.

# the subset function is faster than subsetting the rows via [ ]

sep.data=function(census1, sepcol = "sp", handle.na = NA) 
{
	sepvar = which(colnames(census1) == sepcol)
	sepclass = sort(unique(census1[, sepvar]),na.last=handle.na)
	nclass = length(sepclass)
	output = list()
	for (i in 1:nclass) {
		cat("working on ", sepclass[i], "\n")
		if (!is.na(sepclass[i]))
			output[[i]] = subset(census1,census1[, sepvar]==sepclass[i])
		else
			output[[i]]=subset(census1,is.na(census1[, sepvar]))
	}
	names(output) = sepclass
	return(output)
}

##################################
# VERSION INFORMATION
# created by PHall on 4/14/2005
# similar function by Rick Condit

# sep.vct.data()
#
# Separates a data.frame into a list of data.frames by
# the using the values in a vector of equal length as the data.frame
# The vector is made outside this function in any way desired.

# The output is a list of data.frames.  Each data.frame is composed of the
# records of the trees that match a single value of the sepvct variable.
# The value of sepvct is not provided for the data.frame unless it is a
# variable from the data.frame.  For example, if separating on quadrate,
# then there is one data.frame for each quadrate.

# Datafile and sepvect names must be in quotes.

# the subset function is faster than subsetting the rows via [ ]

sep.vct.data=function(census1, sepvct, handle.na=NA) 
{
 	sepclass = sort(unique(sepvct),na.last=handle.na)
	nclass = length(sepclass)
	output = list()
	for (i in 1:nclass) {
	cat("working on ", sepclass[i], "\n")
	if (!is.na(sepclass[i]))
		output[[i]] = subset(census1,sepvct==sepclass[i])
	else
		output[[i]]=subset(census1,is.na(sepvct))
	}
	names(output) = sepclass
	return(output)
}

###########################
# VERSION INFORMATION
# created and modified in CTFS workshops 2002-2004
# as provided on CTFS/workshopweb in April, 2005

# splitdata() renamed version of split.data()
# the use of "split" in a function name separated by a "."  is not allowed
# in R since "split" is a reserved name.  Therefore, this function will be
# removed.  In the meantime it remains as "splitdata" without the
# separating ".".
#
# Splits a single large data.frame into a list of data.frames by
# the split variable, which could be, for example, species or growth form.
# The split variable must be exactly as long as the full dataframe. If desired
# a vector of unique categories can also be submitted, in order to force the
# result to include empty categories (eg, a missing species value).

# The output is a list of data.frames, each composed of the records of the trees 
# that match a single value of the split variable. For example, if splitting on
# species, then there is one data.frame for each species. Other possibilities:
# split by quadrat number (the program create.quadratlistdata in densdepend.r does
# this), growth form, dbh category, habitat.

# It accepts the data as a file, that is, inside quote marks.

splitdata=function(census1,split=census1$sp,allsplit=unique(split))
{
 output=list()

 for(i in 1:length(allsplit))
  {
   cat("working on number ", i, " -- ", allsplit[i], "\n")
   output[[i]]=subset(census1,split==allsplit[i])
  }

 names(output)=allsplit

 return(output)
}
##############################################################
# sep.quadinfo.R

# this functions creates a vector the same length at the full dataset of a
# given site which contains the habitat that each tree is located in based
# on quadrate location.

# gxgy.to.index was changed to return NA instead of -1, so the loop is
# unnecessary
# A "split" or "separate" variable is a categorical variable by
# which trees will be grouped for a function.  This vector is
# useed to compute values of the dynamics functions: growth.r, mortality.r, recruitment.r,
# abundance.r, abundance.change.r etc.

sep.quadinfo=function(census1,table.in,by.col="hab",gridsize=20,plotdim=c(1000,500))
{
	quad=gxgy.to.index(census1$gx,census1$gy,gridsize=gridsize,plotdim=plotdim)
#	colno=(grep(by.col,names(table.in)))
	return(table.in[quad,by.col])
}
###########################
# sep.sppinfo.R

# this functions creates a vector the same length at the full dataset of a
# given site which contains the growth form of each tree as determined by species

# A "split" or "separate" variable is a categorical variable by
# A "split" or "separate" variable is a categorical variable by
# which trees will be grouped for a function.  This vector is
# useed to compute values of the dynamics functions: growth.r, mortality.r, recruitment.r,
# abundance.r, abundance.change.r etc.

# match returns NA for trees with a species name not found on the table.in
# file.  

sep.sppinfo=function(census1,table.in,by.col="grform")
{
	if (!is.null(dim(census1))) sp=census1$sp
	else sp=census1
	spp=match(sp,table.in$sp,nomatch=NA)
	return(table.in[spp,by.col])
}
###########################
# sep.dbh.R

# this functions creates a vector the same length at the full dataset of a
# given site which contains the growth form of each tree as determined by species

# A "split" or "separate" variable is a categorical variable by
# which trees will be grouped for a function.  This vector is
# useed to compute values of the dynamics functions: growth.r, mortality.r, recruitment.r,
# abundance.r, abundance.change.r etc.

# match returns NA for trees with dbh == NA or dbh outside of the
# dbhclasses <NA> is the value in the returned factor vector and this is
# treated diffently than NA since NA is NOT a level in the factor, unless 
# the it is explicitly defined so such as 
# dbhclass <- factor(dbhclass,exclude=NULL)

sep.dbh=function(census1,dbhclass=c(10,100,300))
{
	if (!is.null(dim(census1))) dbh=census1$dbh
	else dbh=census1
  allbreaks=c(dbhclass,10000)
  dbhclass=cut(dbh,breaks=allbreaks,right=FALSE,
              labels=paste(dbhclass,allbreaks[-1],sep="."))
		  
   return(dbhclass)
}

###########################
# create.status()

# Function for creating a "status" variable that provides a character code to
# indicate the condition of the tree for that census.dbh in the census can be
# a valid measurement or a numeric code, change the dbh as needed to NA so
# that a numeric code is not treated as a number for computation.

#  Function accepts a dataframe with data for a single census.

#  Rules for assigning status and changing dbh value
#  dbh >=10  status=A  tree is alive for this census, keep the dbh
#  dbh >0 & dbh<10  status=AS tree is alive, did have dbh=10 during some census but has
#	since shrunk.  Keep the dbh
#  dbh = 0   status=AB  tree is alive but is <10 mm at point of measurement, dbh=NA
#  dbh = -1  status=D  tree is dead, dbh=NA
#  dbh = -2  status=P  tree is not yet large enough for measurment, a future recruit, dbh=NA
#  dbh = -9  status=NA  missing data, unknown as to cause
#  the dbh of all trees that are D, P or NA are assigned NA
#  EXCEPT
#  at Pasoh, dbh = 5 (or any dbh>0 & dbh<10) status=A and dbh=NA 
#  the tree is known to have dbh>=10  but the dbh was not measured (for unknown reasons).  
#  Therefore, the tree should be counted for stem dynamics, but NOT for growth

#  a label for the column names for each variable that reflect a census number can be appended to the column name if desired by the user.

create.status=function(census,pasoh=FALSE)
{
#  locate and assign variable for dbh, ignores any suffix or prefix on variable name
	
	dbh=census[,grep("dbh",colnames(census))]
	
# all values of status initially set to "A"
	status=rep("A",dim(census)[1])

#  compute status value
	
	status[is.na(dbh)]=NA
	status[dbh==(-9) & !is.na(dbh)]=NA
	status[dbh==(-2) & !is.na(dbh)]="P"
	status[dbh==(-1) & !is.na(dbh)]="D"
	status[dbh==0 & !is.na(dbh)]="AB"
	status[dbh>0 & dbh<10 & !is.na(dbh)]="AS"

	if(pasoh) status[dbh==5 & !is.na(dbh)]="A"
	
#  set dbh = NA when status is "P", "D", "NA", "AB"
#  keep dbh when status is "A", "AS"
	
	dbh[dbh<=0 & !is.na(dbh)]=NA
 	if(pasoh) dbh[dbh==5 & !is.na(dbh)]=NA 

# create a new dataframe and label columns, keep variables as characters NOT factors
	
	newcensus<-data.frame(census$tag,census$sp,census$gx,census$gy,dbh,census$pom,census$date,I(census$codes),I(status))
	names(newcensus)<-c("tag","sp","gx","gy","dbh","pom","date","codes","status")				

	return(newcensus)
}
###########################
# VERSION INFORMATION
# authors: Rick Condit, Suzanne Lao
# create date: workshop 2003
# modifications:  Pam Hall
# last modified data: 5/2005
# changes since last version:
#	modified to conform to CTFS R package standards

# mergecensus()
#
# Function to merge two tree datasets that have the EXACT same tag number
#
#  Function for creating a "status" variable that provides a code to 
#  indicate the condition of the tree for that census.
#  Rules for assigning status and changing dbh value
#  dbh >=10  status=A  tree is alive for this census, keep the dbh
#  dbh >0 & dbh<10  status=AS tree is alive, did have dbh=10 during some census but has
#	since shrunk.  Keep the dbh
#  dbh = 0   status=AB  tree is alive but is <10 mm at point of measurement, dbh=NA
#  dbh = -1  status=D  tree is dead, dbh=NA
#  dbh = -2  status=P  tree is not yet large enough for measurment, a future recruit, dbh=NA
#  dbh = -9  status=NA  missing data, unknown as to cause
#  the dbh of all trees that are D, P or NA are assigned NA
#  EXCEPT
#  at Pasoh, dbh = 5 (or any dbh>0 & dbh<10) status=A and dbh=NA 
#  the tree is known to have dbh>=10  but the dbh was not measured (for unknown reasons).  
#  Therefore, the tree should be counted for stem dynamics, but NOT for growth


mergecensus=function(census0,census1,pasoh=FALSE)
{
      status0=status1=character()
      census0$codes=as.character(census0$codes)
      census1$codes=as.character(census1$codes)

#set status for first census file
	status0[census0$dbh>=10]="A"

	status0[census0$dbh>0 & census0$dbh<10]="AS"
	
	status0[census0$dbh==0]="AB"
	census0$dbh[census0$dbh==0]=NA

	status0[census0$dbh==-1]="D"
	census0$dbh[census0$dbh==-1]=NA
	census0$pom[census0$dbh==-1]= NA  

	status0[census0$dbh==-2]="P"
	census0$dbh[census0$dbh==-2]=NA
	census0$pom[census0$dbh==-2]=NA

	status0[census0$dbh==-9]=NA
	census0$dbh[census0$dbh==-9]=NA
	census0$pom[census0$dbh==-9]=NA

	census0$date[census0$date == 0] = NA
	
# set status for second census file
	status1[census1$dbh>=10]="A"

	status1[census1$dbh>0 & census1$dbh<10]="AS"

	status1[census1$dbh==0]="AB"
	census1$dbh[census1$dbh==0]=NA

	status1[census1$dbh==-1]="D"
	census1$dbh[census1$dbh==-1]=NA
	census1$pom[census1$dbh ==-1]=NA

	status1[census1$dbh==-2]="P"
	census1$dbh[census1$dbh==-2]=NA
	census1$pom[census1$dbh==-2]=NA 

	status1[census1$dbh==-9]=NA
	census1$dbh[census1$dbh==-9]=NA
	census1$pom[census1$dbh==-9]=NA 

	census1$date[census1$date == 0] = NA
	
# exception for the above codes
#  At Pasoh, dbh=5 is for trees that were alive AND were >=10 but were not measured
#  so these trees are considered alive but without a valid dbh for growth computation

	if (pasoh) {
		status0[census0$dbh>0 & census0$dbh<10]="A"
		census0$dbh[census0$dbh>0 & census0$dbh<10]=NA
		census0$pom[census0$dbh > 0 & census0$dbh < 10]=NA
		status1[census1$dbh>0 & census1$dbh<10]="A"
		census1$dbh[census1$dbh>0 & census1$dbh<10]=NA
		census1$pom[census1$dbh > 0 & census1$dbh < 10]=NA
	}
		
# merge variables and rearrange column order

	output=data.frame(tag=census0$tag, sp=I(census0$sp), gx=census0$gx, gy=census0$gy,
	dbh0=census0$dbh, dbh1=census1$dbh,
	pom0=census0$pom, pom1=census1$pom,
	date0=census0$date, date1=census1$date,
	codes0=I(census0$codes), codes1=I(census1$codes),
	status0=I(status0), status1=I(status1))	
      
	return(output)
}

##################################################################
# texttordata()
#
# creates an rdata file from a text file
# makes the first column the rownames, effectively
# duplicating the information in the first column

# this can also be used to read in any text file
# as.is=T   keeps the variable type from the text file so that
#           numbers and characters remain as such
# header=T  reads in the first row as the column names
# sep="\t"  is the separater of the values, in this case tabs
# missvalue allows user to define the missing value, can be a string

# this is useful for the spplist file which contains
# information for each species such as genus, family,
# growth form, reproductive system, etc.  Other functions that use
# the resulting rdata file require that species names be one of the columns
# not just rownames.

texttordata=function(text.in,missvalue=NA)
{
 data.out=read.table(text.in,as.is=TRUE,header=TRUE,sep="\t",quote="",
 	comment.char="",na.strings=missvalue) 
 return(data.out)
}

##########################################
# A function for extracting a single species' or a group of species
# dataframes from the large spp dataset (list of dataframes, one per
# species).  The names of the list data file must include the names
# provides in spname.  A dataframe is returned when spname is only one
# name.  A list of dataframes is returned when spnames is > 1 name.

# load.species() only returns a data frame of single species

load.species=function(spname,census.spp)
 {
 index=which(names(census.spp)==spname)
 gc()
 return(census.spp[[index]])
 }

sep.species=function(spname,census.spp)
{
		index=which(names(census.spp)%in%spnames)
		if (length(spname) == 1)
				return(census.spp[[index]])
		else
				return(census.spp[index])
		}			
				

##########################################
# A function to read a table of elevation data into a matrix appropriate for
# mapping. It can read a text table (if text==TRUE) or a data.frame. 
# In either case, there must be x then y coordinates, followed by elevation, 
# at every corner across the plot.


readelevdata=function(elevfile,gridsize=5,text=FALSE)
{
 if(text) elevdata=read.table(elevfile,header=T,sep="\t")
 else elevdata=elevfile
 
 xdim=max(elevdata$x)
 ydim=max(elevdata$y)

 elevmat=matrix(elevdata$elev,nrow=1+ydim/gridsize,ncol=1+xdim/gridsize,byrow=F)

 return(elevmat)
}
######################################
# A function that converts a data frame of x, y coordinates and elevation 
# to a list containing 2 objects: the original data frame and a matrix 
# form of the same data.  The matrix format is used for making topographic 
# maps and computing mean topographic values for quadrates.

# note the names of the columns in the dataframe must be:
# x y elev

# object 1 in list is: the input dataframe
# object 2 in list is: matrix of elevation values in the orientation of the plot

elev.to.list=function(elevfile,gridsize=5)
{ 
 xdim=max(elevfile$x)
 ydim=max(elevfile$y)

 elevmat=matrix(elevfile$elev,nrow=1+ydim/gridsize,ncol=1+xdim/gridsize,byrow=F)

 return(list(col=elevfile,mat=elevmat))
}

#######################
# Takes output of a demographic analysis and converts into one dataframe.
# Only indicated dbh categories are included; be sure that whichdbhcat does 
# not exceed the number of columns in the output. Type is g for growth, m
# for mortality, and a for abundance
 ##

assemble.demography=function(output,type="g",whichdbhcat=1,datefirst="1980-1-1")
{
 noclass=length(output$rate)
 len=1:noclass
 julian0=tojulian(datefirst,"%Y-%m-%d")

 if(type=="g") 
  {
   result=data.frame(output$rate[,whichdbhcat],output$clim[,whichdbhcat],output$N0[,whichdbhcat],
                     output$dbhmean[,whichdbhcat],output$meanyrs[,whichdbhcat],
                     I(fromjulian(output$date0[,whichdbhcat]+julian0,"%d%b%Y")),
                     I(fromjulian(output$date1[,whichdbhcat]+julian0,"%d%b%Y")))
#   colnames(result)=c("rate","clim","N","time","dbhmean","start","end")
	colnames(result)=c("rate","clim","N0","dbhmean","meanyrs","date0","date1")
  }
  
 else if(type=="m")
  {
   result=data.frame(output$rate[,whichdbhcat],output$upper[,whichdbhcat],output$lower[,whichdbhcat],
                     output$N0[,whichdbhcat],output$N0[,whichdbhcat]-output$D[,whichdbhcat],
                     output$dbhmean[,whichdbhcat],output$meanyrs[,whichdbhcat],
                     I(fromjulian(output$date0[,whichdbhcat]+julian0,"%d%b%Y")),
                     I(fromjulian(output$date1[,whichdbhcat]+julian0,"%d%b%Y")))
   colnames(result)=c("rate","upper","lower","N0","S","dbhmean","meanyrs","date0","date1")
  }
                    
 else if(type=="r")
  {
   result=data.frame(output$rate[,whichdbhcat],output$upper[,whichdbhcat],output$lower[,whichdbhcat],
                     output$N1[,whichdbhcat],output$R[,whichdbhcat],
										 output$meanyrs[,whichdbhcat],
                     I(fromjulian(output$date0+julian0,"%d%b%Y")),
                     I(fromjulian(output$date1+julian0,"%d%b%Y")))
#   colnames(result)=c("N2","R","rate","clim.up","clim.low","start","end")
	 colnames(result)=c("rate","upper","lower","N1","R","meanyrs","date0","date1")
    
  }

	else if(type=="a")
  {
   result=data.frame(output$rate[,whichdbhcat],output$prob[,whichdbhcat],
                     output$N0[,whichdbhcat],output$N1[,whichdbhcat],
										 output$N0[,whichdbhcat] - output$S[,whichdbhcat],output$R[,whichdbhcat],
                     output$meanyrs[,whichdbhcat],
                     I(fromjulian(output$date0[,whichdbhcat]+julian0,"%d%b%Y")),
                     I(fromjulian(output$date1[,whichdbhcat]+julian0,"%d%b%Y")))
#   colnames(result)=c("N.1","N.2","interval","little.r","start","end")
	 colnames(result)=c("rate","prob","N0","N1","D","R","meanyrs","date0","date1")
  }
                    

 return(result)
}

###############################

# Takes output of a demographic analysis, converts
# rate, CLs, and N into one dataframe. Only indicated dbh categories
# are included.

assemble.table.each=function(output,type="g",whichdbhcat=1:2)
{
 noclass=length(output$rate)
 len=1:noclass

 if(type=="g")
  return(data.frame(rate=output$rate[,whichdbhcat],clim=output$clim[,whichdbhcat],
                    N0=output$N0[,whichdbhcat]))

 else if(type=="m")
  return(data.frame(rate=output$rate[,whichdbhcat],lower=output$lower[,whichdbhcat],
                    upper=output$upper[,whichdbhcat],N0=output$N0[,whichdbhcat],
                    S=output$S[,whichdbhcat],meanyrs=output$meanyrs[,whichdbhcat]))

}



