# Function to draw map of one or more species from one plot with the 
# option to add contour lines.  Includes maptopo.r that can be used 
# independently to just draw contour lines or to add contour line toa 
# species plot.

# other functions
# map1species.r = maps the trees, to be used only as a call from map.r
# setsize.r = set size of plotting point for trees

# The full Latin binomial can be used for the map legend by including the 
# spp.info file.

# The size of the graphics window and the size of the plot's graph itself
# can be adjusted with arguments ht, wd, and plotside, but there are
# defaults that produce maps of appropriate aspect.

# Symbol sizes and colors can be entered, or there are defaults.
# default colours = c("blue","red","black","green","yellow","lightblue","orange","pink")

# The title is entered as an argument, or left blank. A legend showing 
# species codes or Latin names is optional. The position for title and
# legend can also be entered via legpos, or there are defaults.

# The many arguments allow flexibility. Default values may fail for
# certain plot shapes, species abundances, so users may have to explore.
# Especially difficult are the graph size, through ht, wd, and plotside,
# since these values have to be adjusted depending on the size of the
# computer's monitor.

# Exporting graphic files may be platform dependent.  This version of 
# map() has recommendations emp, bmp, jpg and pdf that have been well 
# tested in some versions of Windows and MacOSX.  Saved graphics files are written to 
# the current working directory

# the following heights and widths are recommended:
# for emf files (in inches): ht=6, wd=9
# for bmp and jpg files (in pixels):
#    for 1000x500 m plots:   ht=500, wd=850
#    for 500x500 m plots:    ht=480, wd=480
#    for pdf ht=6, wd=8.5, plotside=4.
# for pdf (MacOSX):
#			for 1000x500m plots: ht=7,wd=10

map=function(census.spp,spname,elevmat=NULL,plotdim=c(1000,500),
		spp.info=NULL,maintitle=NULL,titlepos="default",addlegend=TRUE,legpos="default",
		clrs="default",symb="default",size=rep(-1,3),dbhclass=c(10,100,300,3000),
		ht=6,wd=9,plotside=4.5,export="no",bw=FALSE,topoint=0)
{
 clrs1=latin=character()
 symb1=numeric()
 spp.size=numeric() 
 nospp=length(spname)
 
 # get full latin name
  if(!is.null(spp.info))
    for(i in 1:nospp) 
      latin[i]=paste(spp.info[spname[i],]$genus,spp.info[spname[i],]$species) 
	else latin=spname

 # export file or cat to console
	 # this needs to be made completely generic
 if(export!="no") on.exit(graphics.off())
 
 if(export=="no") {	
 		if (getOption("device") == "quartz") 
					quartz(width=wd+1,height=ht+1)
			if (getOption("device") == "windows") {
					dev.set(1)
					win.graph(width=wd,height=ht)
			}
	} 
	else {
	  if(nospp==1) filename=paste(spname,export,sep=".")
		else filename=paste("multisppmap",export,sep=".")

		if(export=="bmp") bmp(filename,width=wd,height=ht)
		else if(export=="jpg") jpeg(filename,width=wd,height=ht,quality=100)
		else if(export=="emf") win.metafile(filename,width=wd,height=ht)
		else if(export=="pdf") pdf(filename,width=wd+1,height=ht+1)
	}

 # number and scale of dbh size class symb1
 # scale to species with largest number of individuals or the first 
 # species is only 1 is provided
	
 for (i in 1:nospp) 
 		spp.size[i]=dim(census.spp[[which(names(census.spp) == spname[i])]])[1]	
		
 nogrp=length(dbhclass)-1
 if(size[1]<0) size=setsize(nogrp,max(spp.size))

 # colours used for tree symb1, for multiple species colours are used in the order provided.

 if(clrs[1]=="default")  {
  if (bw==FALSE) 
	    clrs1=c("blue","red","black","green","yellow","lightblue","orange","pink")
   else 
	    clrs1=c(gray(seq(0,.8,len=nospp)))
   }
 else clrs1=clrs

# define symbols for tree for trees
 if(symb[1]=="default")   {
   if(bw==FALSE)  symb1=rep(16,nospp)
   else symb1=rep(c(16,1),nospp/2+1)
  }
 else symb1=symb

# producing map, 1 species at a time, either as separate maps or combined onto one
 for(i in 1:nospp)
  {
	 sppdata = census.spp[[which(names(census.spp) == spname[i])]]
   sppdata=subset(sppdata,gx>=0 & gy>=0 & gx<plotdim[1] & gy<plotdim[2]) 

   if(i==1) 
	 		map1species(sppdata,sppclrs=clrs1[i],symbol=symb1[i],size=size,
                 dbhclass=dbhclass,plotdim=plotdim,plotside=plotside)
   else 
	 		map1species(sppdata,sppclrs=clrs1[i],symbol=symb1[i],dbhclass=dbhclass,
                    size=size,plotdim=plotdim,plotside=plotside,add=TRUE)
  }

# putting elevation topo lines on map
	if (!is.null(elevmat))
   	maptopo(elev=elevmat,plotdim=plotdim,add=TRUE,interval=topoint)

 oldpar=par(xpd=TRUE,font=3)

 # these 2 have to go after the call to plot inside map1species and maptopo
 # title
 if(titlepos[1]=="default")
 		titlepos=c(plotdim[1]/2,1.16*plotdim[2])
 text(x=titlepos[1],y=titlepos[2],labels=maintitle,col="darkblue",font=2,cex=1.25)

# legend
 if(addlegend)
 	if(legpos[1]=="default")
			legpos=c(plotdim[1]/2,-75)
 legend(x=legpos[1],y=legpos[2],legend=latin,xjust=0.5,col=clrs1,pch=16,horiz=TRUE,y.intersp=1.2)
 
 par(oldpar)
}

# Contours for any plotsite.  Can be run independently of map.r. Elevmat
# has elevation data in matrix form.  Plot dimensions are only needed if
# maptopo.r is called independently of map.r

# The ht, wd and plotside parameters are only used if ONLY the topo 
# contour map is being produced.  Otherwise, these values are set in the 
# main map function.  If add =TRUE then a new window is opened and ht, wd, 
# plotsize are needed to control window size.

# running function independently:
#  export = define type of file to be save to current working dir
#  maintitle, ht, wd, plotside, filename are needed for independent 
#  operation
#  add = add topo contour to an existing plot - irrelevant for independent 
#  operation

# running function from call in map.r
#  add = flag for adding contour map to existing plot

maptopo=function(elevmat,plotdim=c(1000,500),add=FALSE,export="no",interval=0,
                 maintitle="Topographic Contours",ht=6,wd=9,plotside=4.5,filename="maptopo")
 {    

 # export file or display in default graphics device
  if (export!="no") on.exit(graphics.off())
	
	if (export=="no" & !add) {	
 		if (getOption("device") == "quartz") 
					quartz(width=wd+1,height=ht+1)
			if (getOption("device") == "windows") {
					dev.set(1)
					win.graph(width=wd,height=ht)
			}
	}
 else {
		  filename=paste(filename,export,sep=".")
			if(export=="bmp") bmp(filename,width=wd,height=ht)
			else if(export=="jpg") jpeg(filename,width=wd,height=ht,quality=100)
			else if(export=="emf") win.metafile(filename,width=wd,height=ht)
			else if(export=="pdf") pdf(filename,width=wd+1,height=ht+1)
	}
	
	# compute scales for elevation and dimensions of plot
	
	lowelev=min(elevmat)
  hielev=max(elevmat)
  r.elev=hielev-lowelev

  xdim=plotdim[1]
  ydim=plotdim[2]
  aspect=xdim/ydim
  gridsize=round(xdim/dim(elevmat)[2],0)
    
  if(!add) oldpar=par(cex=1,col="black",pin=c(aspect*plotside,plotside))

  if(interval==0)
   {
    if(r.elev>100) interval=5
    else interval=2
   }

  contour(x=seq(0,xdim,by=gridsize),y=seq(0,ydim,by=gridsize),
          t(elevmat),levels=seq(lowelev,hielev,by=interval),add=add,
          drawlabels=FALSE,col="black",axes=FALSE,main=maintitle)

# put axis labels on 2 sides			
  if(!add)   {
    axis(side=1,pos=0,at=c(seq(0,xdim,by=100),xdim),cex=1.15)
    axis(side=2,pos=0,at=c(seq(0,ydim,by=100),ydim),cex=1.15)
    axis(side=3,pos=ydim,at=c(seq(0,xdim,by=100),xdim),labels=FALSE)
    axis(side=4,pos=xdim,at=c(seq(0,ydim,by=100),ydim),labels=FALSE)

		#axis(side=3,pos=ydim,at=c(seq(0,xdim,by=100),xdim),cex=1.15)
    #axis(side=4,pos=xdim,at=c(seq(0,ydim,by=100),ydim),cex=1.15)
   } 
   
  if(!add) par(oldpar)
 }

# This does the mapping for a single species, called from map(). It isn't designed
# to be used on its own.

map1species=function(sppdata,plotdim=c(1000,500),sppclrs="black",symbol=16,size=rep(-1,3),
                     dbhclass=c(10,100,300,3000),plotside=3,add=FALSE)
 {
  xdim=plotdim[1]
  ydim=plotdim[2]
  aspect=xdim/ydim

  nogrp=length(dbhclass)-1
  grp=matrix(FALSE,nrow=nogrp,ncol=dim(sppdata)[1])

  if(size[1]<0) size=setsize(nogrp,length(sppdata$gx))

  for(i in 1:nogrp) grp[i,]=(sppdata$dbh>=dbhclass[i] & sppdata$dbh<dbhclass[i+1])

  oldpar=par(pch=symbol,col=sppclrs,col.main="black",
             pin=c(aspect*plotside,plotside),col.axis="black")

  if(add) 
	points(sppdata$gx[grp[1,]],sppdata$gy[grp[1,]],cex=size[1],col=sppclrs)
  else plot(sppdata$gx[grp[1,]],sppdata$gy[grp[1,]],xlim=c(0,xdim),ylim=c(0,ydim),
            axes=FALSE,xlab="",ylab="",cex=size[1],col=sppclrs)

# put axis labels on 2 sides	
  if(!add)
   {
    axis(side=1,pos=0,at=c(seq(0,xdim,by=100),xdim),cex.axis=1.15)
    axis(side=2,pos=0,at=c(seq(0,ydim,by=100),ydim),cex.axis=1.15)
    axis(side=3,at=c(seq(0,xdim,by=100),xdim),labels=FALSE,pos=ydim)
    axis(side=4,at=c(seq(0,ydim,by=100),ydim),labels=FALSE,pos=xdim)
		
#		axis(side=3,pos=ydim,at=c(seq(0,xdim,by=100),xdim),cex=1.15)
#   axis(side=4,pos=xdim,at=c(seq(0,ydim,by=100),ydim),cex=1.15)
   }
  if(nogrp>1)
    for(i in 2:nogrp) 
		points(sppdata$gx[grp[i,]],sppdata$gy[grp[i,]],cex=size[i],col=sppclrs)

  par(oldpar)
 }




# Chooses size of plotting points, depending on a species' abundance. Note that the user
# can override these defaults by submitting sizes.

setsize=function(n,s)
 {
  if(s>1000)     size=0.4+.3*(0:(n-1))
  else if(s>400) size=.55+.3*(0:(n-1))
  else           size=0.7+.3*(0:(n-1))

  if(n==3)
   if(s>1000)
    {
     size[1]=.4
     size[2]=.7
     size[3]=1.1
    }
   else if(s>400)
    {
     size[1]=.55
     size[2]=.8
     size[3]=1.2
    }
   else
    {
     size[1]=.7
     size[2]=1
     size[3]=1.2
    }

  return(size)
 }
