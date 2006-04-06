# VERSION INFORMATION ver 1.0
# authors: Rick Condit
# create date: 2003 CTFS workshop
# modifications: Pam Hall

# functions included:
# mortality()
# mortality.eachspp()
# mortality.dbh()
# mortality.calculation()

# Overview of functions for calculating mortality rates.

# mortality() is the main function, and is constructed like 
# growth() and recruitment(). It requires two complete datasets, one per census,
# with variables named dbh, pom, status, and date for every individual of all species in at least 2 censuses. 
# It can then calculate mortality based on up to user-submitted factors. The two
# datasets have exactly the same individuals, in exactly the same order, one 
# individual per row.

# Mortality is based on status; any status indicating a live tree can be
# submitted.  Any tree without an alivecode in census 2 is considered dead.
# Individuals whose status is NA in either census are deleted from all
# calculations.

mortality=function(census1,census2,alivecode=c("A","AB","AS"),split1=NULL,split2=NULL)
{
 if(is.null(split1)) split1=rep("all",dim(census1)[1])
 if(is.null(split2)) split2=rep("all",dim(census2)[1])

 inc=!is.na(census1$status) & !is.na(census2$status)
 census1=census1[inc,]
 census2=census2[inc,]
 split1=split1[inc]
 split2=split2[inc]

 timeint=(census2$date-census1$date)/365.25

 alive1=alive2=rep(FALSE,dim(census1)[1])
 alive1[census1$status=="A"]=TRUE
 for(i in 1:length(alivecode)) alive2[census2$status==alivecode[i]]=TRUE

 class1=sort(unique(split1))
 class2=sort(unique(split2))
 splitN0=list(split1[alive1],split2[alive1])
 splitS=list(split1[alive1&alive2],split2[alive1&alive2])

 N0=tapply(census1$dbh[alive1],splitN0,length)
 S=tapply(census1$dbh[alive1&alive2],splitS,length)
 meantimeint=tapply(timeint[alive1],splitN0,mean,na.rm=T)
 meandbh=tapply(census1$dbh[alive1],splitN0,mean,na.rm=T) 
 startdate=tapply(census1$date[alive1],splitN0,mean,na.rm=T)
 enddate=tapply(census2$date[alive1],splitN0,mean,na.rm=T)

 N0=fill.dimension(N0,class1,class2)
 S=fill.dimension(S,class1,class2)
 meantimeint=fill.dimension(meantimeint,class1,class2,fill=NA)
 meandbh=fill.dimension(meandbh,class1,class2,fill=NA)
 startdate=fill.dimension(startdate,class1,class2,fill=NA)
 enddate=fill.dimension(enddate,class1,class2,fill=NA)

 if(sum(N0)==0) {
   result=list(N0=rep(NA,length(class1)),S=rep(NA,length(class1)),D=rep(NA,length(class1)),
               rate=rep(NA,length(class1)),
               lower=rep(NA,length(class1)),upper=rep(NA,length(class1)),
               dbhmean=rep(NA,length(class1)),meanyrs=rep(NA,length(class1)),
               date0=rep(NA,length(class1)),date1=rep(NA,length(class1)))
		gc()
		return(result)
    }

# Rick's changes in 9/05 which deletes find.climits.smallN and creates mortality.calculation()
	
 m=mortality.calculation(N0=N0,S=S,meantimeint=meantimeint)

# lower.ci=upper.ci=N
# noS=(S==0)
# lower.ci[!noS]=find.climits(N[!noS],(N-S)[!noS],kind="lower")
# lower.ci[noS]=find.climits.smallN(N[noS],(N-S)[noS],kind="lower")
 
# noD=(N-S==0)
# upper.ci[!noD]=find.climits(N[!noD],(N-S)[!noD],kind="upper")
# upper.ci[noD]=find.climits.smallN(N[noD],(N-S)[noD],kind="upper")

# mort.rate=(log(N)-log(S))/meantimeint
# upper.rate=(log(N)-log(N-upper.ci))/meantimeint
# lower.rate=(log(N)-log(N-lower.ci))/meantimeint

# mort.rate[S==0]=upper.rate[S==0]=Inf
# upper.rate[upper.ci==N]=Inf
# mort.rate[N==0]=lower.rate[N==0]=upper.rate[N==0]=NA

# result=list(N=N,D=N-S,rate=mort.rate,lower=lower.rate,upper=upper.rate,
#	 dbhmean=meandbh,timeint=meantimeint,date1=startdate,date2=enddate)
# change order from Ricks 9/05 version to keep standard with other dynamic functions

 result=list(N0=N0,S=S,D=N0-S,rate=m$rate,lower=m$lowerCI,upper=m$upperCI,dbhmean=meandbh,
		 meanyrs=meantimeint,date0=startdate,date1=enddate)
 gc()
 return(result)
}



# A cover function to calculate mortality for each species in given dbh categories.

mortality.eachspp=function(census1,census2,dbhclass=c(10,100,300),alivecode=c("A","AB","AS"))
{
# allbreaks=c(dbhbreaks,10000)
# dbhclass=cut(census1$dbh,breaks=allbreaks,right=F,
#              labels=paste(dbhbreaks,allbreaks[-1],sep="."))
 
 sp=census1$sp
 dbh.vct=sep.dbh(census1,dbhclass=dbhclass)				
 result=mortality(census1,census2,alivecode=alivecode,
	 split1=sp,split2=dbh.vct)
 gc()
 return(result)
}


# A cover function to calculate forest-wide mortality in given dbh categories.

mortality.dbh=function(census1,census2,dbhclass=c(10,100,300),alivecode=c("A","AB","AS"))
{
# allbreaks=c(dbhclass,10000)
# dbhbreaks=cut(census1$dbh,breaks=allbreaks,right=F,
#              labels=paste(dbhclass,allbreaks[-1],sep="."))

 dbh.vct=sep.dbh(census1,dbhclass=dbhclass)				
 result=mortality(census1,census2,alivecode=alivecode,
	 split2=dbh.vct)
 gc()
 return(result)
}

# A simple function to calculate mortality and confidence limits from N 
# (number alive at the outset), S (number of survivors), and timeint (time interval).
# They can be vectors or scalars. 

mortality.calculation=function(N0,S,meantimeint)
{
 lower.ci=find.climits(N0,S,kind="lower")
 upper.ci=find.climits(N0,S,kind="upper")

# find.climits returns the upper and lower number DEAD trees.
# to calculate mortality, S has to be computed from the upper and lower
# number of D, hence, S=N0-upperD and S=N0-lowerD

 mort.rate=(log(N0)-log(S))/meantimeint
 upper.rate=(log(N0)-log(lower.ci))/meantimeint
 lower.rate=(log(N0)-log(upper.ci))/meantimeint

 mort.rate[S==0]=upper.rate[S==0]=Inf
 upper.rate[upper.ci==N0]=Inf
 lower.rate[lower.ci==N0]=0
 mort.rate[N0==0]=lower.rate[N0==0]=upper.rate[N0==0]=NA

 if(is.null(dim(N0)))
   return(data.frame(N0=N0,S=S,D=N0-S,rate=mort.rate,lowerCI=lower.rate,upperCI=upper.rate,
			 meantimeint=meantimeint))
 else
   return(list(N0=N0,S=S,D=N0-S,rate=mort.rate,lowerCI=lower.rate,upperCI=upper.rate,
			 meantimeint=meantimeint))
}

# lower.ci=upper.ci=N
# noS=(S==0)
# lower.ci[!noS]=find.climits(N[!noS],(N-S)[!noS],kind="lower")
# lower.ci[noS]=find.climits.smallN(N[noS],(N-S)[noS],kind="lower")
# noD=(N==S)
# upper.ci[!noD]=find.climits(N[!noD],(N-S)[!noD],kind="upper")
# upper.ci[noD]=find.climits.smallN(N[noD],(N-S)[noD],kind="upper")

