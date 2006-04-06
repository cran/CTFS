# VERSION INFORMATION
# authors: Rick Condit
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
# growth()
# growth.eachspp()
# growth.dbh()
# growth.indiv()
# trim.growth()

# depends: utilities2005.R
# fill.dimensions()
# rounddown()

# can be used with
# assemble.demography()

# Overview for computing growth rates

# Growth is the main function, and is constructed like recruitment and
# mortality.  It requires two complete datasets, one per census, with dbh,
# pom, and date for every individual of all species in at least 2 censuses.
# It can then calculate mortality based on up to user-submitted factors.
# The two datasets have exactly the same individuals, in exactly the same
# order, one individual per row.

# Growth relies on status and pom. The function trim.growth handles all
# checks for trees to include -- excludes code R, excludes extreme values
# based on err.limit and maxgrow, excludes trees below a minimum dbh in the
# first census.

# It can calculate growth for any subsection of plot data, as determined by
# the variables split1 and split2.

growth=function(census1,census2,rounddown=FALSE,method='I',stdev=FALSE,
                mindbh=10,err.limit=4,maxgrow=75,split1=NULL,split2=NULL)
{
 if(is.null(split1)) split1=rep("all",dim(census1)[1])
 if(is.null(split2)) split2=rep("all",dim(census2)[1])

# if(is.null(census2$codes)) census2$codes=rep(".",length(census2$dbh))

 timeint=(census2$date-census1$date)/365.25

 if(rounddown)
  {
   sm=(census1$dbh<55 | census2$dbh<55)
   census1$dbh[sm]=rndown5(census1$dbh[sm])
   census2$dbh[sm]=rndown5(census2$dbh[sm])
  }

 if(method=='I') growth=(census2$dbh-census1$dbh)/timeint
 else if(method=='E') growth=(log(census2$dbh)-log(census1$dbh))/timeint
 good=trim.growth(census1$dbh,census2$dbh,census1$pom,census2$pom,census1$status,census2$status,
                  timeint,err.limit=err.limit,maxgrow=maxgrow)
 good[census1$dbh<mindbh]=FALSE
 growth[!good]=NA

 class1=sort(unique(split1))
 class2=sort(unique(split2))
 splitgood=list(split1[good],split2[good])

 mean.grow=tapply(growth[good],splitgood,mean,na.rm=T)
 sd.grow=tapply(growth[good],splitgood,sd,na.rm=T)
 N0=tapply(growth[good],splitgood,length)
 meandbh=tapply(census1$dbh[good],splitgood,mean,na.rm=T) 
 meantimeint=tapply(timeint[good],splitgood,mean,na.rm=T)
 startdate=tapply(census1$date[good],splitgood,mean,na.rm=T)
 enddate=tapply(census2$date[good],splitgood,mean,na.rm=T)

 mean.grow=fill.dimension(mean.grow,class1,class2,fill=NA)
 sd.grow=fill.dimension(sd.grow,class1,class2,fill=NA)
 N0=fill.dimension(N0,class1,class2,fill=0)
 meandbh=fill.dimension(meandbh,class1,class2,fill=NA)
 meantimeint=fill.dimension(meantimeint,class1,class2,fill=NA)
 startdate=fill.dimension(startdate,class1,class2,fill=NA)
 enddate=fill.dimension(enddate,class1,class2,fill=NA)

 ci.grow=sd.grow
 ci.grow[N0==0]=NA
 ci.grow[N0>0]=sd.grow[N0>0]*qt(0.975,N0[N0>0])/sqrt(N0[N0>0])

 if(stdev==F) 
    result=list(N0=N0,rate=mean.grow,clim=ci.grow,dbhmean=meandbh,
                meanyrs=meantimeint,date0=startdate,date1=enddate)
 else if(stdev==T) 
    result=list(N0=N0,rate=mean.grow,sd=sd.grow,dbhmean=meandbh,
                meanyrs=meantimeint,date0=startdate,date1=enddate)
 
 gc()
 return(result)
}


# A cover function to calculate growth for each species in given dbh categories.

growth.eachspp=function(census1,census2,dbhclass=c(10,100,300),
	rounddown=FALSE,method="I",stdev=FALSE,err.limit=4,maxgrow=75)
{
# allbreak=c(classbreak,10000)
# dbhclass=cut(census1$dbh,breaks=allbreak,right=F,
#              labels=paste(classbreak,allbreak[-1],sep="."))

 dbh1.vct=sep.dbh(census1$dbh,dbhclass=dbhclass)
 sp=census1$sp
 
 result=growth(census1,census2,rounddown=rounddown,
               method=method,err.limit=err.limit,maxgrow=maxgrow,
               stdev=stdev,split1=sp,split2=dbh1.vct)
 gc()
 return(result)
}


# A cover function to calculate forest-wide growth in given dbh categories.

growth.dbh=function(census1,census2,dbhclass=c(10,100,300),
	rounddown=FALSE, method="I",stdev=FALSE,err.limit=4,maxgrow=75)
{
# allbreak=c(classbreak,10000)
 dbh1.vct=sep.dbh(census1$dbh,dbhclass=dbhclass)
 
 result=growth(census1,census2,rounddown=rounddown,
               method=method,err.limit=err.limit,maxgrow=maxgrow,
               stdev=stdev,split1=dbh1.vct)

 gc()
 return(result)
}


# Growth.indiv calculates the growth rate of each individual, returning 2
# different growth measures in a data.frame.  The data returned only
# included records where growth could be calculated, with tag number,
# coordinates, initial dbh, and the growth data.

growth.indiv=function(census1,census2,rounddown=FALSE,mindbh=10,err.limit=4,maxgrow=75)
{
# if(is.null(census2$codes)) census2$codes=rep(".",length(census2$dbh))

 timeint=(census2$date-census1$date)/365.25

 if(rounddown)
  {
   sm=(census1$dbh<55 | census2$dbh<55)
   census1$dbh[sm]=rndown5(census1$dbh[sm])
   census2$dbh[sm]=rndown5(census2$dbh[sm])
  }

 incgrowth=(census2$dbh-census1$dbh)/timeint
 expgrowth=(log(census2$dbh)-log(census1$dbh))/timeint
 good=trim.growth(census1$dbh,census2$dbh,census1$pom,census2$pom,census1$status,census2$status,
                  timeint,err.limit=err.limit,maxgrow=maxgrow)
 good[census1$dbh<mindbh]=FALSE
 incgrowth[!good]=expgrowth[!good]=NA

 growthdata=data.frame(tag=census1$tag,sp=census1$sp,gx=census1$gx,gy=census1$gy,
                       dbh1=census1$dbh,dbh2=census2$dbh,years=timeint,
                       incgr=incgrowth,expgr=expgrowth)

 return(growthdata)
}

##############################################################

# Trim.growth is a function to exclude growth rates. It is based on 
# a linear model estimating the standard deviation of dbh measures (due to error, that
# is); the parameters slope and intercept define the linear relationship between
# this error deviation and dbh. Any case where the second dbh measure is more than
# 4 standard deviations below the first is marked false so as to exclude from
# growth measurements. The default values of slope and intercept are based on dbh
# remeasure tests done in both 1995 and 2000. A line was fitted through the absolute 
# dbh errors as a function of dbh in both years; the average slope and intercept is
# used here. It also excludes any growth rate > 75 mm per yr, cases
# where the POM changes, and if either dbh<10.

# trim.growth() relies upon the variable "code" to identify "reprouts" with
# "R".  However, this is not used consistently among sites and furthermore
# is not well defined and "R" may be used in a two letter code to indicate
# some other condition.  The condition of resprouts should be captured in a
# combination of status and pom.  If the dbh measurement of a census is in
# a different location the next census, then dbh growth cannot be calculated for that stem.
# Location is defined as different pom on the same stem, a different
# stem for a multiple stemmed individual or a resprout.  A resprout could
# be an entire new stem from basal roots or a new shoot from a broken stem.

trim.growth=function(dbh1,dbh2,pom1,pom2,status1,status2,timeint,
                     slope=0.006214,intercept=.9036,err.limit=4,maxgrow=75)				 
#trim.growth=function(dbh1,dbh2,pom1,pom2,code2,timeint,
#                     slope=0.006214,intercept=.9036,err.limit=4,maxgrow=75)
{
 stdev.dbh1=slope*dbh1+intercept

 growth=(dbh2-dbh1)/timeint

 bad.neggrow=which(dbh2<=(dbh1-err.limit*stdev.dbh1)) 
 bad.posgrow=which(growth>maxgrow)

 accept=rep(TRUE,length(dbh1))
# accept[grep("R",code2)]=FALSE
 accept[status1 != "A" | status2 != "A"]=FALSE
 accept[is.na(growth)]=FALSE
 accept[pom1!=pom2]=FALSE
 accept[dbh1<10 | dbh2<10]=FALSE

 accept[bad.neggrow]=FALSE
 accept[bad.posgrow]=FALSE

 return(accept)
}


