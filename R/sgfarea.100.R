# VERSION INFORMATION
# authors: Suzanne Lao, Rick Condit, Pam Hall
# create date: 2003 CTFS workshop
# modifications: Pam Hall
# last modified date: 5/2005
# changes since last version:
#	removed use of get() to access datafile
#	input 2 separate census files instead of merged census files
#	compute start and end dates, add more values to output list
#	correction: removed trees from NA in computation
#	modified to conform to CTFS R package standards

# functions included:
# sqf.area() number of ind, spp, gen, fam in a single quadrate, x,y as parameters
# sgf.randarea() number of ind, spp, gen, fam in a randomly chosen quadrate
# sgf.meanarea() mean number of ind, spp, gen, fam, in N randomly chosen quadrates
# sgf.meanarea.range() mean number of ind, spp, gen, fam, in N randomly chosen quadrates of a range of sizes

# spp.randarea.range() number and plot of species in N randomly chosen quadrates of a range of sizes
# select.randquad() select coordinates of a random quadrate of a given size

# Function for calculating the number of individuals, species, genera, and families 
# in a quadrat with the coordinates specified.
# default is the entire plot

 sgf.area=function(census1,spp.info,corners=c(0,1000,0,500),
                     alivecode=c("A"),mindbh=NULL,
					 unidennames=c("**","UNID","uniden","UNIDEN"))
 {
 
  incdbh=rep(TRUE,dim(census1)[1])
  if(!is.null(mindbh)) incdbh[census1$dbh<mindbh | is.na(census1$dbh)]=FALSE	  
  
  alive=rep(FALSE,dim(census1)[1])
  for(i in 1:length(alivecode)) alive[census1$status==alivecode[i]]=TRUE
  
  incarea=rep(FALSE,dim(census1)[1])
  if(length(corners)!=0) incarea[census1$gx>=corners[1] & census1$gx<corners[2] 
  	& census1$gy>=corners[3] & census1$gy<corners[4]]=TRUE

  no.ind=length(subset(census1$sp,alive & incdbh & incarea 
	  & !unidentified.species(census1$sp,unidenstr=unidennames)))
  
  rownames(spp.info)=spp.info$sp
  
  spp=unique(census1$sp[alive & incdbh & incarea])
  spp=spp[!unidentified.species(spp,unidenstr=unidennames) & !is.na(spp)]

  gen=unique(spp.info[spp,]$genus)
   gen=gen[!is.na(gen)]

  fam=unique(spp.info[spp,]$family)
  fam=fam[!is.na(fam)]

  no.spp=length(spp)
  no.gen=length(gen)
  no.fam=length(fam)

  sgr=no.spp/no.gen
  sfr=no.spp/no.fam
  gfr=no.gen/no.fam
  
  result=(c(no.ind, no.spp, no.gen, no.fam, sgr, sfr, gfr))
  names(result)=c("n.ind", "n.spp", "n.gen", "n.fam", "r.sppgen", "r.sppfam", "r.genfam")
  return(result)
 }

# This computes the numbers of individuals, species, genera, and families in a 
# randomly chosen quadrat of a specified size
# Uses: sgf.area() 
# default is 1 ha quadrate

sgf.randarea=function(census1,spp.info,quaddim=c(100,100),
     alivecode=c("A"),mindbh=NULL,plotdim=c(1000,500),
	 unidennames=c("**","UNID","uniden","UNIDEN"))
 {
   if(quaddim[1]>plotdim[1] | quaddim[2]>plotdim[2]) 
     return("You can't use a quadrat larger than the plot")

   xlo=(plotdim[1]-quaddim[1])*runif(1)
   ylo=(plotdim[2]-quaddim[2])*runif(1)
   xhi=xlo+quaddim[1]
   yhi=ylo+quaddim[2]

   result=sgf.area(census1,spp.info,corners=c(xlo,xhi,ylo,yhi),
		   alivecode=alivecode,mindbh=mindbh,unidennames=c("**","UNID","uniden","UNIDEN"))
   return(result)	   
 }

# This calculates the mean numbers of individuals, species, genera, and families 
# for a given number of randomly chosen quadrats of a specified size
# defaults are 1 ha quadrates and 50 replicates

sgf.meanarea=function(census1,spp.info,quaddim=c(100,100),
	alivecode=c("A"),mindbh=NULL,plotdim=c(1000,500),replicates=50,
    unidennames=c("**","UNID","uniden","UNIDEN"))
 {
   sgf=matrix(ncol=7, nrow=replicates)
   
   for (i in 1:replicates)
    sgf[i,]=sgf.randarea(census1,spp.info,quaddim,alivecode,mindbh,plotdim,unidennames) 
    
  result=c(mean(sgf[,1]),sd(sgf[,1]),
            mean(sgf[,2]),sd(sgf[,2]),
            mean(sgf[,3]),sd(sgf[,3]),
            mean(sgf[,4]),sd(sgf[,4]),
            mean(sgf[,5]),sd(sgf[,5]),
            mean(sgf[,6]),sd(sgf[,6]),
            mean(sgf[,7]),sd(sgf[,7]))
		   
  names(result)=c("mn.nind", "sd.nind", "mn.nspp", "sd.nspp", "mn.ngen", "sd.ngen",
	   "mn.nfam", "sd.nfam", "mn.rsppgen", "sd.rsppgen", "mn.rsppfam", "sd.rsppgen",
	   "mn.rgenfam", "sd.rgenfam")
  return(result)
 }

# This function computes what "sgf.meanarea" does for a range of specified quadrat sizes.

sgf.meanarea.range=function(census1,spp.info,xlengths=c(50,100,200,500),
                ylengths=c(50,100,200,500),alivecode=c("A"),mindbh=NULL,
                plotdim=c(1000,500),replicates=50,unidennames=c("**","UNID","uniden","UNIDEN"))
 {
  no.quadrats=length(xlengths)  
  sgf.per.quad=matrix(ncol=15, nrow=no.quadrats+1)
  
  for(i in 1:no.quadrats)
   {
    cat("working on quadrat ", xlengths[i], "x", ylengths[i], "\n")

    hectares=xlengths[i]*ylengths[i]/10000
       
    sgf.per.quad[i,1]=hectares
    sgf.per.quad[i,2:15]=sgf.meanarea(census1,spp.info,quaddim=c(xlengths[i],ylengths[i]),
                                         alivecode,mindbh,plotdim,replicates,unidennames) 
    }

    cat("working on full plot\n")

    lastquad=no.quadrats+1 
    sgf.per.quad[lastquad,1]=plotdim[1]*plotdim[2]/10000
    sgf.per.quad[lastquad,c(2,4,6,8,10,12,14)]=
            sgf.area(census1,spp.info,corners=c(0,plotdim[1],0,plotdim[2]),alivecode,mindbh,unidennames)

    colnames(sgf.per.quad)=c("area", "mn.nind", "sd.nind", "mn.nspp", "sd.nspp", "mn.ngen", "sd.ngen",
                             "mn.nfam", "sd.nfam", "mn.rsppgen", "sd.rsppgen", 
                             "mn.rsppfam", "sd.rsppgen","mn.rgenfam", "sd.rgenfam")
   return(data.frame(sgf.per.quad))
}

# Function for calculating the number of species in replicated, randomly-placed
# quadrats of various areas. The variable quadside refers to the dimension (side)
# of the square, and can be a vector of as many desired; replicates is the
# number of random draws per dimension. Full plot data are submitted as censdata.
# The program removes unidentified records, based on a default unidennames; see
# the function unidentified.species in utilities.r.

# spp.randarea.range()
# select.randquad() called by spp.randarea.range()


spp.randarea.range=function(census1,quadside=100,shape=1,alivecode=c("A"),
	mindbh=NULL,plotdim=c(1000,500),
    replicates=10,unidennames=c("**","UNID","uniden","UNIDEN"))
{
 toolong=quadside[quadside>plotdim[1]]
 quadside=quadside[quadside*shape<=plotdim[2]]

 if(length(toolong)>0) 
   {
    cat("quadside ", min(toolong))
    if(length(toolong)>1) cat(", ", max(toolong))
    cat (" too big for plot\n")
   }

  alive=rep(FALSE,dim(census1)[1])
  for(i in 1:length(alivecode)) alive[census1$status==alivecode[i]]=TRUE
  
  incdbh=rep(TRUE,dim(census1)[1])
  if(!is.null(mindbh)) incdbh[census1$dbh<mindbh | is.na(census1$dbh)]=FALSE	  
  
  incarea=rep(FALSE,dim(census1)[1])
  if(length(plotdim)!=0) incarea[census1$gx>=0 & census1$gx<plotdim[1] 
  	& census1$gy>=0 & census1$gy<plotdim[2]]=TRUE

 census1=subset(census1,alive & incdbh & incarea &
      !unidentified.species(sp,unidenstr=unidennames))
 
 allquadside=select.randquad(quadside,shape,plotdim,replicates)
 noquad=dim(allquadside)[1]
 
 spp=ind=numeric()
 for(i in 1:noquad)
  {
   data=subset(census1,gx>=allquadside$xlo[i] & gx<allquadside$xhi[i] & 
                        gy>=allquadside$ylo[i] & gy<allquadside$yhi[i])
 
   spp[i]=length(unique(data$sp))
   ind[i]=length(data$sp)

   if(i==1) plot(allquadside$area[i],spp[i],pch=16,ylim=c(1,1200),xlim=c(.01,60),log="xy")
   else points(allquadside$area[i],spp[i],pch=16)
  }

 full=data.frame(area=allquadside$area,n.spp=spp,n.ind=ind)
 spp=tapply(full$n.spp,full$area,mean,na.rm=TRUE)
 SDspp=tapply(full$n.spp,full$area,sd,na.rm=TRUE)
 indiv=tapply(full$n.ind,full$area,mean,na.rm=TRUE)
 SDindiv=tapply(full$n.ind,full$area,sd,na.rm=TRUE)
 no.area=length(spp)
 area=as.numeric(row.names(spp))
 
 spp[no.area+1]=length(unique(census1$sp))
 indiv[no.area+1]=length(census1$sp)
 SDspp[no.area+1]=SDindiv[no.area+1]=NA
 area[no.area+1]=plotdim[1]*plotdim[2]/1e4
 xdim=c(quadside,plotdim[1])
 ydim=c(quadside*shape,plotdim[2])
 
 gc()
 return(list(spparea=data.frame(xdim,ydim,area,mn.nind=indiv,sd.nind=SDindiv,mn.nspp=spp,sd.nspp=SDspp),full=full))
}


# Returns a dataframe of coordinates defining the corners of randomly
# chosen square quadrats.  The variable quadside is a vector of dimensions;
# rep is the number of replicated random quadrats per dimension submitted.
# shape is the value to multiple x dim by to get y dim.  when shape = 1, it
# is a square.

select.randquad=function(quadside,shape = 1,plotdim=c(1000,500),replicates = 100)
{
 noquad=length(quadside)
 allquadside.x=rep(quadside,each=replicates)
 allquadside.y=rep(quadside*shape,each=replicates)
 total=length(allquadside.x)
 
 xlo=runif(total,min=0,max=plotdim[1]-allquadside.x)
 ylo=runif(total,min=0,max=plotdim[2]-allquadside.y)
 xhi=xlo+allquadside.x
 yhi=ylo+allquadside.y

 return(data.frame(xlo,ylo,xhi,yhi,area=(allquadside.x*allquadside.y)/1e4))
}

