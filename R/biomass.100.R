# biomass functions

# the main biomass() functions calls:
# 	pick.dbh()
#   getwsg.vct()
#   agb.ind(), agb.mult(), gr.raw()
#   growth.dbh(), trim.growth()

#  other biomass functions
#   badgrowth.tree() - to list trees with estimated DBHs
#   biomass.change() - compute change in AGB between censuses
#   biomass.dbh.change() - compute change in AGB between census for dbh
#   classes (assessed for each census)

# what biomass() does:
# reads 2 census files, 2 multstem files, wsg(wood density) file and sets other values for computing AGB
# compute growth to determine if the dbh is bad
# 	return a dbh for AGB computation  for census (dbhagb)
# compute AGB for the main stem
#	use wood density (lookup species name)
#   use forest type: wet, moist, dry
#	return AGB for main stem (AGBM)
# compute AGB for multi stems
#	get species from census (match $tag)
# 	use wood density (lookup species name)
#	return AGB for multi stems (AGBS)
# compute total AGB for individual (AGBT)
# return tag, sp, gx, gy, dbh1, dbh2, growth, AGBm0, AGBs0, AGBt0, AGBm1, AGBs1, AGBt1

# take care with setting the values of the arguments.  make sure all
# arguments that have an "=" have some value you set or accept the
# defaults.  Note that the defaults for file names (mult1, mult2,
# wsgfile) are NULL. this means that NO data on multiple stems for either
# census or wood density is used by default.  YOU SHOULD PROVIDE THE FILE
# NAMES FOR EACH OF THESE ARGUMENTS!

# example - set values for estimating dbh that are NOT default but leave
# the default values for determining what is a "bad" growth which will
# require estimating a dbh.  the values set by the user for estimating dbh
# will be used for trees with "bad" growth

# badgrowth "keep", "nogrowth", "meangrowth", "trim"
# # dbhuse "first" "second"
# forest "wet" "moist" "dry"  NULL (default)

biomass=function(census1,census2,mult1=NULL,mult2=NULL,wsgfile=NULL,
	forest=NULL,badgrowth="meangrowth",dbhuse="second",trim=c(-5,75),recsize=110,rounddown=FALSE)
{
	
# check dbh for valid values and estimate invalid dbh when possible
	estdbh=pick.dbh(census1,census2,badgrowth=badgrowth,dbhuse=dbhuse,trim=trim,recsize=recsize,rounddown=rounddown)
	dbh1=estdbh$dbh1
	dbh2=estdbh$dbh2
	status1=estdbh$status1
	pom1=estdbh$pom1
	GRm=estdbh$GRm
	GRraw=estdbh$GRraw
		
# get wood density for main stem based on species 
# need to compute mean for species found at a given site, not for entire wsg file
	  
  if (is.null(wsgfile))
    return("need a wood density file name or a single value for all trees\cr")
  else if (length(wsgfile)==1)	
	wsg.vct=rep(wsgfile,length(dbh1))
  else 
  	wsg.vct=getwsg.vct(census1,wsgfile)
   
# compute AGB for main, mult stems and total for each census
	
  AGBm1=AGBs1=AGBt1=numeric()
  AGBm2=AGBs2=AGBt2=numeric()

  AGBm1=agb.ind(dbh1,wsg.vct,forest=forest)
  AGBm2=agb.ind(dbh2,wsg.vct,forest=forest)
  
  if (is.null(mult1))
  	AGBs1=rep(0,length(census1$tag))
  else
  	AGBs1=agb.mult(census1,mult1,wsgfile,forest=forest)
	
  if (is.null(mult2))
	  AGBs2=rep(0,length(census2$tag))
  else
    AGBs2=agb.mult(census2,mult2,wsgfile,forest=forest)

  AGBs1<-as.vector(AGBs1)
  AGBs2<-as.vector(AGBs2)
  AGBt1=AGBm1+AGBs1
  AGBt2=AGBm2+AGBs2
  
#  merge AGBm and AGBs to outdata and compute AGBt
#  convert all AGB to Mg (metric tons)	
  
  out.data=data.frame(census1$tag,estdbh$errtype,AGBm1/1000,AGBm2/1000,AGBs1/1000,AGBs2/1000,AGBt1/1000,AGBt2/1000)
  names(out.data)=c("tag","errtype","AGBm1","AGBm2","AGBs1","AGBs2","AGBt1","AGBt2")
  return(out.data)
}

##################################################################
# pick.dbh() depends on: gr.raw(), growth(), growth.dbh()
# gr.raw() 

# provides 6 options for determining the dbh of trees that have "bod" growths
# defaults set for use outside of biomass.R

pick.dbh=function(census1,census2,badgrowth="keep",dbhuse="second",trim=c(-5,75),recsize=110,rounddown=FALSE)
{
 dbh1=census1$dbh
 dbh2=census2$dbh
 pom1=census1$pom
 pom2=census2$pom
 date1=census1$date
 date2=census2$date
 status1=census1$status
 status2=census2$status
 
# compute mean annual growth in mm for dbhclasses using trimmed growth as
# conditioned by BCI values and identify the mean growth rate for the dbh
# class of each tree in each census
 

 dbhclass=c(10,20,50,100,200,300,400,500,600,700)
 gr.out=growth.dbh(census1,census2,dbhclass=dbhclass)
 mngr.dbhnames=dimnames(gr.out$rate)[[1]]
 dbh1cl.ind=sep.dbh(dbh1,dbhclass=dbhclass)
 dbh2cl.ind=sep.dbh(dbh2,dbhclass=dbhclass) 
 mv1=match(dbh1cl.ind,mngr.dbhnames,nomatch=NA)
 mv2=match(dbh2cl.ind,mngr.dbhnames,nomatch=NA)
 mngr.dbh1=gr.out$rate[mv1]
 mngr.dbh2=gr.out$rate[mv2]
 time=(date2-date1)/365.25
 
# date1[is.na(date1)] = date2 - (time*365.25)
# date2[is.na(date2)] = date1 + (time*365.25)
# date1[date1==0] = date2 - (time*365.25)
# date2[date2==0] = date1 + (time*365.25)
 
# old code
# dbhclass=c(10,20,50,100,200,300,400,500,600,700,10000)
# ndbh=length(dbhclass)
# gr.out=growth.dbh(census1,census2,dbhclass=dbhclass[-ndbh])
# dbh1cl.ind=cut(dbh1,breaks=dbhclass,right=F,labels=paste(dbhclass[-ndbh],dbhclass[2:ndbh],sep="."))
# dbh2cl.ind=cut(dbh2,breaks=dbhclass,right=F,labels=paste(dbhclass[-ndbh],dbhclass[2:ndbh],sep="."))

# compute growth to determine if the dbh for AGB is "good"
# trees with good growth will retain both census dbh
# this includes trees that have a change in pom because while the growth
# itself is not valid both census dbhs are
# set errtype 
#  "NE" no error
#  "FG" growth too fast, > trim[2]
#  "SG" growth too slow, < trim[1]
#  "M1" missing dbh1 even though status1 == A
#  "M2" missing dbh2 even though status2 == A
#  "BR" too big for recruit
#  "PC" POM changed - not technically an error, but tree cannot be used to compute growth of AGB

 GRm=gr.raw(census1,census2,rounddown=rounddown)
 GRraw=GRm
 good.growth=rep(TRUE,length(GRm))
 errtype=rep("NE",length(GRm))
 
 good.growth[(pom1 == pom2) & (GRm > trim[2])]=FALSE
 errtype[(pom1 == pom2) & (GRm > trim[2])]="FG"
 good.growth[(pom1==pom2)	& (GRm < trim[1])]=FALSE
 errtype[(pom1==pom2)	& (GRm < trim[1])]="SG"

# to generalize to different census, use the census number for the sample function 
# if badgrowth=="keep", do nothing, dbhagb is already assigned to cens dbh 

 if (badgrowth=="nogrowth") {
	 GRm[!good.growth]=0
	 if (dbhuse=="first")  
	 	dbh2[!good.growth]<-dbh1[!good.growth]
	 else if (dbhuse=="second") 
	 	dbh1[!good.growth]<-dbh2[!good.growth]
	 else 
	 	return("invalid dbhuse value provided\cr")
 }
 else if (badgrowth=="meangrowth") {
	 if (dbhuse=="first") {
	 	 GRm[!good.growth]<-mngr.dbh1[!good.growth]
		 dbh2[!good.growth]<-dbh1[!good.growth]+GRm[!good.growth]*time[!good.growth]
	 }
	 else if (dbhuse=="second") {
		 GRm[!good.growth]<-mngr.dbh2[!good.growth]
		 dbh1[!good.growth]<-dbh2[!good.growth]-GRm[!good.growth]*time[!good.growth]
	 }
	 else
	 	return("invalid dbhuse value provided\cr")
 }
 else if (badgrowth=="trim") {
	 GRm[!good.growth & GRm > trim[2]] = trim[2]
	 GRm[!good.growth & GRm < trim[1]] = trim[1]

	 if (dbhuse=="first") 
		dbh2[!good.growth]=dbh1[!good.growth]+GRm[!good.growth]*time[!good.growth]
	 else if (dbhuse=="second") 
		dbh1[!good.growth]=dbh2[!good.growth]-GRm[!good.growth]*time[!good.growth]
     else
		return("invalid dbhuse value provided\cr")
	}
	
# identify trees that are alive, but have a missing dbh in one census and estimate missing dbh.
# these trees are distinct from those with a bad.growth since with a missing dbh, growth cannot be calculated to assess it
# at Pasoh, the status AS is treated as a living tree (status A) and 
# sufficiently large for the census, but the dbh was not measured in the 
# year that status was == AS, so it should be estimated.  
	miss.dbh1=rep(FALSE,length(dbh1))
	miss.dbh1[((status1=="A" | status1=="AS") & is.na(dbh1)) & (status2=="A" & !is.na(dbh2))]=TRUE
	errtype[((status1=="A" | status1=="AS") & is.na(dbh1)) & (status2=="A" & !is.na(dbh2))]="M1"
	dbh1[miss.dbh1]=dbh2[miss.dbh1]-mngr.dbh2[miss.dbh1]*time[miss.dbh1]
	
	miss.dbh2=rep(FALSE,length(dbh2))
	miss.dbh2[(status1=="A" & !is.na(dbh1)) & ((status2=="A" | status2=="AS") & is.na(dbh2))]=TRUE
	errtype[(status1=="A" & !is.na(dbh1)) & ((status2=="A" | status2=="AS") & is.na(dbh2))]="M2"
	dbh2[miss.dbh2]=dbh1[miss.dbh2]+mngr.dbh1[miss.dbh2]*time[miss.dbh2]
	
# identify trees that are recruits to the second census but have a dbh2 that is too large to be a recruit
# estimate dbh1 for bad recruits based on the mean growth rate of the dbhclass for dbh2
	
	bad.rec=rep(FALSE,length(dbh1))
	bad.rec[(status1=="P" & status2=="A") & (dbh2 >= recsize)]=TRUE
	errtype[(status1=="P" & status2=="A") & (dbh2 >= recsize)]="LR"
	dbh1[bad.rec]=dbh2[bad.rec]-mngr.dbh2[bad.rec]*time[bad.rec]
	status1[bad.rec]="A"
	pom1[bad.rec]=pom2[bad.rec]
	
# identify trees with a change in POM, but do not estimate dbh, tree has 2
# valid dbhs, but since POM changed, cannot use them to compute growth
# for all calculation except growth of persistent trees, these trees 
# should be considered alive and present
  errtype[status1 == "A" & status2 == "A" & pom1 != pom2] = "PC"  
 
 result=data.frame(census1$tag,errtype,dbh1,dbh2,status1,census2$status,pom1,census2$pom,GRm,GRraw)
 colnames(result)=c("tag","errtype","dbh1","dbh2","status1","status2","pom1","pom2","GRm","GRraw")
# result=data.frame(tag=census1$tag,errtype=errtype,dbh1=dbh1,dbh2=dbh2,
#				status1=status1,status2=census2$status,pom1=pom1,pom2=census2$pom,
#				GRm=GRm,GRraw=GRraw)
 return(result)		
}

		 
##################################################################
# compute growth for each tree in order to determine a valid dbh for
# biomass computation.  called by pick.dbh().  errlimit and maxgrow set to
# very high values to allow a computation of all growth rates for trees
# with 2 dbhs.

gr.raw=function(census1,census2,rounddown=FALSE)
{
 good=rep(NA,length(census1$dbh))
 time=(census2$date-census1$date)/365.25
 
 if(rounddown)
   {
	sm=rep(FALSE,length(census1$dbh)) 
    sm[census1$dbh<55 | census2$dbh<55]=TRUE
    census1$dbh[sm]=rndown5(census1$dbh[sm])
    census2$dbh[sm]=rndown5(census2$dbh[sm])
	}
 
 gr.ind=(census2$dbh-census1$dbh)/time
 good=trim.growth(census1$dbh,census2$dbh,census1$pom,census2$pom,census1$status,census2$status,
                  time,err.limit=1000,maxgrow=1000)

 gr.ind[!good]=NA
 return(gr.ind)
}

##################################################################
# make a vector of wood density for each tree based on species name (sp code used to match)
# table.in = bci.wood.density
# rows = species, col 2 = wood density
# if the species is not found in the wood density file the density is set
# to the mean of all densities in the table file

getwsg.vct=function(census,wsgfile,by.col="wsg")
{
	
 	sp=census$sp
	spp=match(sp,wsgfile$sp,nomatch=NA)
	wsg.found=wsgfile[spp,by.col]
	
# compute mean wsg for species, not individuals
	
	unq.spp=match(unique(sp),wsgfile$sp,nomatch=NA)
	spp.found=wsgfile[unq.spp,by.col]
	mnwsg=mean(spp.found,na.rm=TRUE)
	wsg.found[is.na(wsg.found)]=mnwsg
	
	return(wsg.found)
	
}

##################################################################
# compute AGB using different models for forest types

agb.ind=function(dbh,wsg.vct,forest=forest)
{
# AGB=p*exp(aprime+b*lndbh+c*lndbh^2+d*lndbh^3)
# aprime=a*RSE^2/2

  RSE=numeric()	
  dbh=dbh/10
  ldbh = log(dbh)
  
  RSE[2]=0.356^2/2
  RSE[1]=0.377^2/2
 
  if (is.null(forest))
  # generic computation for any forest type
  	AGB=wsg.vct*exp(RSE[1]+(-1.589)+(2.284*ldbh)+(0.129*ldbh^2)+(-0.0197*ldbh^3))
  else if (forest=="wet")
   # computation for wet forest type
	AGB=wsg.vct*exp(RSE[2]+(-1.302)+(1.98*ldbh)+(0.207*ldbh^2)+(-0.0281*ldbh^3))
  else if (forest=="moist")
   # computation for moist forest type
	AGB=wsg.vct*exp(RSE[2]+(-1.562)+(2.148*ldbh)+(0.207*ldbh^2)+(-0.0281*ldbh^3))
  else if (forest=="dry")
   # computation for dry forest type
	AGB=wsg.vct*exp(RSE[2]+(-0.73)+(1.784*ldbh)+(0.207*ldbh^2)+(-0.0281*ldbh^3))
  return(AGB)	
}

##################################################################
# agb.mult.r
# merge multfile and census by tag to add the species name from census to multfile
# match to census to get species for wsg.vct
# set wsg for species missing from wsgfile to mean of stems with found species
  
# assign("mult.sp.data",merge(get(multfile),get(datafile)[,1:2],
#		  by="tag",all.x=TRUE),pos=1)
  
agb.mult=function(census,multfile,wsgfile,forest=forest)
{
 
  multfile.sp<-merge(multfile,census[,1:2],by="tag",all.x=TRUE)
  
  if (is.null(wsgfile))
	  return("need a wood density file name or a single value for all trees\cr")
  else if (length(wsgfile)==1)
	wsg.mult.vct<-rep(wsgfile,dim(multfile.sp)[1])
  else 
  	wsg.mult.vct<-getwsg.vct(multfile.sp,wsgfile)

# compute AGB for mult stems

  AGBs.ind=agb.ind(multfile.sp$dbh,wsg.mult.vct,forest=forest)  

# AGBs is the sum of the AGB of all multiple stems for a main stem
# compute the sum of AGB of individual stems for a main stem (tag)
# match the tag number of multstem data to full data
# where there is no match (main stem did not have multstems), set AGB=0
  
  AGBs.tmp=tapply(AGBs.ind, multfile.sp$tag, sum)
  mult.mtch=match(census$tag, names(AGBs.tmp), nomatch=NA)
  AGBs=AGBs.tmp[mult.mtch]
  AGBs[is.na(AGBs)]=0
  
  return(AGBs)
  
}

##################################################################
#  This function identifies trees with "badgrowth" rate, which is any errtype != NE returned by pick.dbh
#  makes an estimation of the dbh that contributed to
#  the bad growth and returns sufficient information about the tree to
#  evaluate this estimation process.
#   pick.dbh=function(census1,census2,badgrowth="keep",dbhuse="first",trim=c(-5,75),recsize=20,rounddown=FALSE)

badgrowth.tree=function(census1,census2,badgrowth="meangrowth",dbhuse="second",trim=c(-5,75),recsize=110,rounddown=FALSE)
{
	estdbh <- pick.dbh(census1,census2,badgrowth=badgrowth,dbhuse=dbhuse,trim=trim,recsize=recsize,rounddown=rounddown)

#	bad.growth=rep(FALSE,length(census1$dbh))
# 	bad.growth[(census1$pom == census2$pom) & (rawGRm > trim[2] | rawGRm < trim[1])]=TRUE

	result=data.frame(census1[c(1:4)],estdbh$errtype,census1$dbh,census2$dbh,estdbh$dbh1,estdbh$dbh2,
		estdbh$GRraw,estdbh$GRm,census1$status,census2$status,estdbh$status1,estdbh$status2,
		census1$pom,census2$pom,estdbh$pom1,estdbh$pom2)
	colnames(result)=c("tag","sp","gx","gy","errtype","dbh1","dbh2","estdbh1","estdbh2","rawGR",
		"estGR","status1","status2","eststatus1","eststatus2","pom1","pom2","estpom1","estpom2")
	
	return(result[estdbh$errtype != "NE",])
}

##################################################################

###############################################################
# biomass.change

# computing mortality, recruitment, ingrowth and persistent biomass for 2 
# censuses
# works like growth(), mortality() etc.  2 split vectors are allowed
#  DO NOT DBHCLASS AS A SPLIT VECTOR, for DBHCLASS use biomass.change.dbh()

# this function relies upon the output from biomass() after any corrections 
# are made to the census1 and census2 datasets and/or corrections based on 
# a choice of dbh estimations.  biomass() provides above 
# ground biomass (AGB) values for any valid dbh.  Status is only changed 
# for "big recruits" (biomass errtype=="LR"), where status1 is changed 
# from = "P" to "A".

# various combinations of status values have to be dealt with here to get 
# a valid measure of AGB dynamics  to whit:

# biomass() computes AGB for every valid dbh (after any needed 
# estimations).  Therefore, this function only needs to identify the 
# presence of a valid AGB for a given census to determine if the tree is 
# alive, dead, a recruit, or persistent.  Ingrowth of AGB can only be 
# calculated from trees that have valid AGB for each census AND have no 
# change in POM.  Trees that are alive in both censuses but have a chang 
# in POM are present, but cannot be used to contribute to growth.  This 
# yields inconsistent results when trying to compute AGB for a given year 
# from previous year AGB, mortality, recruitment and ingrowth.  But the 
# counts of the trees in these categories do add up.

biomass.change=function(biomass.out,split1=NULL,split2=NULL)
{
	AGBt1 = biomass.out$AGBt1
	AGBt2 = biomass.out$AGBt2
	errtype = biomass.out$errtype
	
  if(is.null(split1)) split1=rep("all",dim(biomass.out)[1])
  if(is.null(split2)) split2=rep("all",dim(biomass.out)[1])
	
  class1=sort(unique(split1))
  class2=sort(unique(split2))
  splitboth=list(split1,split2)
 
  splitdead = list(split1[!is.na(AGBt1) & is.na(AGBt2)],split2[!is.na(AGBt1) 
			& is.na(AGBt2)])
	dead = tapply(AGBt1[!is.na(AGBt1) & is.na(AGBt2)],splitdead,sum)
	ndead = tapply(AGBt1[!is.na(AGBt1) & is.na(AGBt2)],splitdead,length)
#	dead = sum(subset(biomass.out,!is.na(AGBt1) & is.na(AGBt2))$AGBt1)

	splitrec = list(split1[is.na(AGBt1) & !is.na(AGBt2)], split2[is.na(AGBt1) 
			& !is.na(AGBt2)])
	rec = tapply(AGBt2[is.na(AGBt1) & !is.na(AGBt2)],splitrec,sum)
	nrec = tapply(AGBt2[is.na(AGBt1) & !is.na(AGBt2)],splitrec,length)
#  rec = sum(subset(biomass.out,is.na(AGBt1) & !is.na(AGBt2))$AGBt2)

  splitaliveboth = list(split1[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],
			split2[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"])
	aliveboth1 = tapply(AGBt1[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth,sum)
	naliveboth1 = tapply(AGBt1[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth,length)
	aliveboth2 = tapply(AGBt2[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth,sum)
	naliveboth2 = tapply(AGBt2[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth,length)
#	aliveboth1 = sum(subset(biomass.out,!is.na(AGBt1) & !is.na(AGBt2) & errtype != "PC" )$AGBt1)
# aliveboth2 = sum(subset(biomass.out,!is.na(AGBt1) & !is.na(AGBt2) & errtype != "PC" )$AGBt2)

  splitalive1 = list(split1[!is.na(AGBt1)],split2[!is.na(AGBt1)])
	alive1 = tapply(AGBt1[!is.na(AGBt1)],splitalive1,sum)
	nalive1 = tapply(AGBt1[!is.na(AGBt1)],splitalive1,length)
#	alive1 = sum(subset(biomass.out, !is.na(AGBt1))$AGBt1)

  splitalive2 = list(split1[!is.na(AGBt2)],split2[!is.na(AGBt2)])
	alive2 = tapply(AGBt2[!is.na(AGBt2)],splitalive2,sum)
	nalive2 = tapply(AGBt2[!is.na(AGBt2)],splitalive2,length)
#	alive2 = sum(subset(biomass.out, !is.na(AGBt2))$AGBt2)

	splitalivePC = list(split1[errtype == "PC"],split2[errtype == "PC"])
	alivePC1 = tapply(AGBt1[errtype == "PC"],splitalivePC,sum)
	nalivePC1 = tapply(AGBt1[errtype == "PC"],splitalivePC,length)
	alivePC2 = tapply(AGBt2[errtype == "PC"],splitalivePC,sum)
	nalivePC2 = tapply(AGBt2[errtype == "PC"],splitalivePC,length)
#	alivePC1 = sum(subset(biomass.out, errtype == "PC" )$AGBt1)
# alivePC2 = sum(subset(biomass.out, errtype == "PC" )$AGBt2)

	alivegr = aliveboth2 - aliveboth1
	alivegrPC = alivePC2 - alivePC1

	dead=fill.dimension(dead,class1,class2,fill=0)
	ndead=fill.dimension(ndead,class1,class2,fill=0)
	rec=fill.dimension(rec,class1,class2,fill=0)
	nrec=fill.dimension(nrec,class1,class2,fill=0)
	aliveboth1=fill.dimension(aliveboth1,class1,class2,fill=0)
	naliveboth1=fill.dimension(naliveboth1,class1,class2,fill=0)
	aliveboth2=fill.dimension(aliveboth2,class1,class2,fill=0)
	naliveboth2=fill.dimension(naliveboth2,class1,class2,fill=0)
	alive1=fill.dimension(alive1,class1,class2,fill=0)
	nalive1=fill.dimension(nalive1,class1,class2,fill=0)
	alive2=fill.dimension(alive2,class1,class2,fill=0)
	nalive2=fill.dimension(nalive2,class1,class2,fill=0)
	alivePC1=fill.dimension(alivePC1,class1,class2,fill=0)
	nalivePC1=fill.dimension(nalivePC1,class1,class2,fill=0)
	alivePC2=fill.dimension(alivePC2,class1,class2,fill=0)
	nalivePC2=fill.dimension(nalivePC2,class1,class2,fill=0)
	
	alivegr=fill.dimension(alivegr,class1,class2,fill=0)
	alivegrPC=fill.dimension(alivegrPC,class1,class2,fill=0)

	result=list(alive1=alive1,alive2=alive2,dead=dead,rec=rec,
			aliveboth1=aliveboth1,aliveboth2=aliveboth2,alivePC1=alivePC1,alivePC2=alivePC2,
			nalive1=nalive1,nalive2=nalive2,ndead=ndead,nrec=nrec,
			naliveboth1=naliveboth1,naliveboth2=naliveboth2,nalivePC1=nalivePC1,nalivePC2=nalivePC2,
			alivegr=alivegr,alivegrPC=alivegrPC)
  gc()
	return(result)
}


###################################################################
# biomass.change.dbh

# computing biomass.change but for dbh classes assessed independently for 
# each census, in the same manner as for abundance

# use biomass() to compute AGB values but run twice each time with a 
# different dbh vct based on different censuses.  

# dbh from census 1 is used for AGB alive1, dead, aliveboth1 and alivePC1
# dbh from census 2 is used for AGB alive2, rec, aliveboth2 and alivePC2
# dbh from census 1 is used for AGB growth: alivegr and alivegrPC
#  this function repeats all the code from biomass.change because 
#  different split variables have to be used since not all trees have a 
#  valid dbh for each census.

biomass.change.dbh=function(census1, census2, biomass.out,dbhclass=c(10,100,300),split=NULL)
{
	AGBt1 = biomass.out$AGBt1
	AGBt2 = biomass.out$AGBt2
	errtype = biomass.out$errtype

	if(is.null(split)) split=rep("all",dim(biomass.out)[1])
	dbh1.vct = sep.dbh(census1,dbhclass=dbhclass)		
	dbh2.vct = sep.dbh(census2,dbhclass=dbhclass)

  class1=sort(unique(dbh1.vct))
  class2=sort(unique(dbh2.vct))
	class3=sort(unique(split))

  splitdead = list(split[!is.na(AGBt1) & is.na(AGBt2)],dbh1.vct[!is.na(AGBt1) 
			& is.na(AGBt2)])
	dead = tapply(AGBt1[!is.na(AGBt1) & is.na(AGBt2)],splitdead,sum)
	ndead = tapply(AGBt1[!is.na(AGBt1) & is.na(AGBt2)],splitdead,length)

	splitrec = list(split[is.na(AGBt1) & !is.na(AGBt2)], dbh2.vct[is.na(AGBt1) 
			& !is.na(AGBt2)])
	rec = tapply(AGBt2[is.na(AGBt1) & !is.na(AGBt2)],splitrec,sum)
	nrec = tapply(AGBt2[is.na(AGBt1) & !is.na(AGBt2)],splitrec,length)

  splitaliveboth1 = list(split[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],
			dbh1.vct[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"])
	splitaliveboth2 = list(split[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],
			dbh2.vct[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"])

	aliveboth1 = tapply(AGBt1[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth1,sum)
	naliveboth1 = tapply(AGBt1[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth1,length)
	aliveboth2 = tapply(AGBt2[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth2,sum)
	naliveboth2 = tapply(AGBt2[!is.na(AGBt1) & !is.na(AGBt2) & errtype != 
			"PC"],splitaliveboth2,length)

  splitalive1 = list(split[!is.na(AGBt1)],dbh1.vct[!is.na(AGBt1)])
	alive1 = tapply(AGBt1[!is.na(AGBt1)],splitalive1,sum)
	nalive1 = tapply(AGBt1[!is.na(AGBt1)],splitalive1,length)

  splitalive2 = list(split[!is.na(AGBt2)],dbh2.vct[!is.na(AGBt2)])
	alive2 = tapply(AGBt2[!is.na(AGBt2)],splitalive2,sum)
	nalive2 = tapply(AGBt2[!is.na(AGBt2)],splitalive2,length)

	splitalivePC1 = list(split[errtype == "PC"],dbh1.vct[errtype == "PC"])
	splitalivePC2 = list(split[errtype == "PC"],dbh2.vct[errtype == "PC"])

	alivePC1 = tapply(AGBt1[errtype == "PC"],splitalivePC1,sum)
	nalivePC1 = tapply(AGBt1[errtype == "PC"],splitalivePC1,length)
	alivePC2 = tapply(AGBt2[errtype == "PC"],splitalivePC2,sum)
	nalivePC2 = tapply(AGBt2[errtype == "PC"],splitalivePC2,length)

	alivegr = aliveboth2 - aliveboth1
	alivegrPC = alivePC2 - alivePC1

	dead=fill.dimension(dead,class3,class1,fill=0)
	ndead=fill.dimension(ndead,class3,class1,fill=0)
	rec=fill.dimension(rec,class3,class2,fill=0)
	nrec=fill.dimension(nrec,class3,class2,fill=0)
	aliveboth1=fill.dimension(aliveboth1,class3,class1,fill=0)
	naliveboth1=fill.dimension(naliveboth1,class3,class1,fill=0)
	aliveboth2=fill.dimension(aliveboth2,class3,class2,fill=0)
	naliveboth2=fill.dimension(naliveboth2,class3,class2,fill=0)
	alive1=fill.dimension(alive1,class3,class1,fill=0)
	nalive1=fill.dimension(nalive1,class3,class1,fill=0)
	alive2=fill.dimension(alive2,class3,class2,fill=0)
	nalive2=fill.dimension(nalive2,class3,class2,fill=0)
	alivePC1=fill.dimension(alivePC1,class3,class1,fill=0)
	nalivePC1=fill.dimension(nalivePC1,class3,class1,fill=0)
	alivePC2=fill.dimension(alivePC2,class3,class2,fill=0)
	nalivePC2=fill.dimension(nalivePC2,class3,class2,fill=0)
	
	alivegr=fill.dimension(alivegr,class3,class1,fill=0)
	alivegrPC=fill.dimension(alivegrPC,class3,class1,fill=0)

	result=list(alive1=alive1,alive2=alive2,dead=dead,rec=rec,
			aliveboth1=aliveboth1,aliveboth2=aliveboth2,alivePC1=alivePC1,alivePC2=alivePC2,
			nalive1=nalive1,nalive2=nalive2,ndead=ndead,nrec=nrec,
			naliveboth1=naliveboth1,naliveboth2=naliveboth2,nalivePC1=nalivePC1,nalivePC2=nalivePC2,
			alivegr=alivegr,alivegrPC=alivegrPC)
  gc()
	return(result)
}
















