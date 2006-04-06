# VERSION HISTORY
# authors: Rick Condit
# create date: 2004 CTFS workshop
# last modified date: Aug 2005
# 2005 - revised to use single census data files
# variable names changed to be consistent with CTFS package
# census.quad = full dataset restructured as list of quads, each quad a dataframe of the trees in that quad
# census.spp = full dataset restructured as a list of species, each species a dataframe of the trees of that species
# spname = a character string of a species name code as given in census.spp and sitespp.info

# functions included:
# neighbordens.1sp()
# neighbordens.1quad()
# countdens.1quad()
# countdens.1tree()
# fill.distclasses() - moved to utilities
# borderdist() - in here but is an undocumented internal function

# countdens() does computes the number of individuals and basal area of 
# neighbors for conspecifics and heterospecifics using the computed 
# distance from each (potential) neighbor to a target individual.  The 
# function cut() is used to determine the correct neighborhood distance 
# class to place individuals.


neighbordens.1sp=function(spname,census.quad,census.spp,distbreak=c(0,5,10,20,30),
                              gridsize=20,plotdim=c(1000,500))
{
 i=0
 on.exit(print(i))
 # spdata=load.species(spname,census.spp)
# spp.index=which(names(census.spp)==spname)
# spdata=census.spp[[spp.index]]
  spdata=census.spp[[which(names(census.spp)==spname)]]

#maxquad=plotdim[1]*plotdim[2]/gridsize^2
 #maxquad=length(names(census))
 
 for(i in as.numeric(names(census.quad)))
  {
   if(i==1) result=neighbordens.1quad(i,spname,census.quad,distbreak,gridsize,plotdim)
   
   else result=rbind(result,neighbordens.1quad(i,spname,census.quad,distbreak,
                                              gridsize,plotdim))
                                            
   if(i%%50==0) cat("Finished quadrat ", i, " of species ", spname, "\n")                                              
  }
  
 m=match(spdata$tag,rownames(result))
 
 return(data.frame(tag=spdata$tag,result[m,]))
}
   
# Calculate density of conspecific and heterospecific neighbors within distance
# annuli for all individuals in one quadrat of a given species. Distance intervals
# must be entered as well.

# quad is the quad index, a number from 1 to 1250

neighbordens.1quad=function(quad,spname,census.quad,distbreak=distbreak,
                           gridsize=gridsize,plotdim=plotdim)
{	
 neighborquads=c(quad,findborderquads(quad,dist=max(distbreak),
                                      gridsize=gridsize,plotdim=plotdim))

 neighborquads=as.character(neighborquads)                                    
 noquads=length(neighborquads)
 quad=as.character(quad)
 
 focaldata=subset(census.quad[[quad]],sp==spname & !is.na(gx) & !is.na(gy))
 nofocal=dim(focaldata)[1]
 
 if(nofocal>0)
  {
   resultlist=lapply(census.quad[neighborquads],countdens.fullquad,focaldata=focaldata,
                     quad=quad,distbreak=distbreak,plotdim=plotdim)
 
   total=matrix(0,nrow=nofocal,ncol=4*(length(distbreak)-1))

   for(i in 1:length(resultlist)) total=total+resultlist[[i]]
   rownames(total)=focaldata$tag
 
   return(total)
  }
  
}


# This function takes a quadrat number, species name, a data frame of 
# neighbors, and a data frame of focal plants, and calculates the number and 
# BA of neighbors in the given distance categories, for each focal plant,
# separating hetero- and con-specifics neighbors.

countdens.1quad=function(neighdata,focaldata,quad,distbreak,plotdim)
{	
 if (is.null(dim(neighdata[1])))
	 return(array(0,dim=c(dim(focaldata)[1],4*(length(distbreak)-1))))  
 else if (dim(neighdata)[1]==0) 
   	return(array(0,dim=c(dim(focaldata)[1],4*(length(distbreak)-1))))  
 else {  
   
 	neighdata=subset(neighdata,!is.na(gx)&!is.na(gy))
 
 	noneighbor=dim(neighdata)[1]
 
 #dbhcol=which(colnames(neighdata)==paste("dbh",cens,sep=""))
 #statuscol=which(colnames(neighdata)==paste("status",cens,sep=""))
  
 	dbh=neighdata$dbh
 	status=neighdata$status
 
	# all neighbors that are alive or alive and broken as counted as 
	# neighbors.  No other selection is done.  If a user wishes a different 
	# selection criteria, this is probably the logical place to put it.  The 
	# logical vector "inc" will define which neighbor distances and dbhs are 
	# used to compute neighborhood attributes.
	
 	inc=(status=="A" | status=="AB") & !is.na(status)
 	consp=neighdata$sp==focaldata$sp[1]
 
 	neighdist=t(apply(data.frame(focaldata$gx,focaldata$gy),1,xydistvect,
                   pt2=data.frame(neighdata$gx,neighdata$gy)))

 	neighdist=array(neighdist[,inc],dim=c(dim(focaldata)[1],length(inc[inc])))
 
 	neighdens=t(apply(neighdist,1,countdens.1tree,dbh=dbh[inc],consp=consp[inc],
                   distbreak=distbreak))

 	border=borderdist(subset(focaldata,select=c(gx,gy)),plotdim)
 	nodistclass=length(distbreak)-1
 	for(i in length(distbreak):2)
   	{
    	close=border<distbreak[i]
    	neighdens[close,seq((i-1),dim(neighdens)[2],by=nodistclass)]=NA
   	}
   
 	return(neighdens)
	}
}
 
# This works with a single focal individual. It accepts a vector of distances 
# to each neighbor, with a vector of 
# dbhs for each neighbor, and returns the number N and basal area BA in each
# distance category, split into conspecific and heterospecific. The latter
# is determined by a logical vector, T for each dbh of a conspecific. 

countdens.1tree=function(distradius,dbh,consp,distbreak)
{
 distlabels=
   paste(distbreak[1:(length(distbreak)-1)],distbreak[2:length(distbreak)],sep=".")

	 # determine which distance class the distances belong, distances are 
	 # classified as:  low_value <= neighbor_dist < high_value
	 # eg:  for a neighbor that is 7.5 m from target individual, with 
	 # distbreaks = c(0,5,10,20,30), distclass = 5.10 because 5 <= 7.5 < 10
	 # eg: for a neighbor that is 10 m from target individual, with 
	 # distbreaks = c(0,5,10,20,30), distclass = 10.20 because 10 <= 10 < 20
	 
	 distclass=cut(distradius,breaks=distbreak,right=F,labels=distlabels)
 
 inc=distradius<max(distbreak) & distradius>0
 
 Ncon=BAcon=Nhet=BAhet=rep(0,length(distlabels))
 
 if(length(dbh[inc&consp])>0)
  {
   Ncon=drop(tapply(dbh[inc&consp],distclass[inc&consp],length))
   BAcon=drop(tapply(dbh[inc&consp],distclass[inc&consp],basum))
  }
 if(length(dbh[inc&!consp])>0)
  {
   Nhet=drop(tapply(dbh[inc&!consp],distclass[inc&!consp],length))
   BAhet=drop(tapply(dbh[inc&!consp],distclass[inc&!consp],basum))
  }

 Ncon=fill.distclasses(Ncon,distlabels)  
 BAcon=fill.distclasses(BAcon,distlabels)  
 Nhet=fill.distclasses(Nhet,distlabels)  
 BAhet=fill.distclasses(BAhet,distlabels)  

 names(Ncon)=paste("Ncon",distlabels,sep=".")
 names(Nhet)=paste("Nhet",distlabels,sep=".")
 names(BAcon)=paste("BAcon",distlabels,sep=".")
 names(BAhet)=paste("BAhet",distlabels,sep=".")
 
 return(c(Ncon,BAcon,Nhet,BAhet))
} 
 

# Calculates the closest plot border, given a plot size and a vector of
# coordinates. Assumes minimum plot coorindates are zero.

borderdist=function(coord,plotdim)
{
 distradius=distx=disty=numeric()
 
 E=coord$gx>plotdim[1]/2
 W=!E
 N=coord$gy>plotdim[2]/2
 S=!N
 
 distx[E]=plotdim[1]-coord$gx[E]
 distx[W]=coord$gx[W]

 disty[N]=plotdim[2]-coord$gy[N]
 disty[S]=coord$gy[S]
 
 distradius[distx<disty]=distx[distx<disty]
 distradius[distx>=disty]=distx[distx>=disty]

 return(distradius)
}
