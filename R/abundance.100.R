# VERSION INFORMATION ver 1.0
# authors: Rick Condit
# created date: 2003 CTFS workshop
# modifications: Pam Hall

# functions included:
# abundance()
# abundance.quad()
# abundance.spp()
# abundance.change()
# abundance.change.dbh() 
# meaninterval.quad()
# abundance.change.calculation()
# abund.ch.prob()

# Overview of functions for calculating abundance and change in population size

# abundance() and its derivative functions, require a single census
# dataset.  The total number of individuals is computed including trees
# with status == "A" and mindbh = 10, or as specified by the user.  It can
# then calculate abundance based on up to 2 user-submitted factors.
# abundance.quad() is a cover function that calculate abundance per user
# defined quadrate.  abundance.spp() is a cover function that calculates
# abundance for each species for dbh classes.  abundance by dbh classes can
# be computed using abundance() and a vector that defines the dbh class of
# each tree.

# abundance.change() is the calculation of the rate of change in abundance.
# It requires 2 census files and determines whether a tree should be
# included based on its status and mindbh for each census INDEPENDENTLY.
# This is NOT a single cohort, but the change in the structure of the
# forest between censuses.  A computation of the probability that N0 = N1
# is also computed.

abundance=function(census1,alivecode=c("A"),mindbh=10,
	split1=NULL,split2=NULL)
{

 if(is.null(split1)) split1=rep("all",dim(census1)[1])
 if(is.null(split2)) split2=rep("all",dim(census1)[1])

 inc=rep(TRUE,dim(census1)[1])
 if(!is.null(mindbh)) inc[census1$dbh<mindbh | is.na(census1$dbh)]=FALSE

 alive=rep(FALSE,length(census1$dbh))
 for(i in 1:length(alivecode)) alive[census1$status==alivecode[i]]=TRUE

 class1=sort(unique(split1))
 class2=sort(unique(split2))
 splitA=list(split1[alive&inc],split2[alive&inc])

 N=tapply(census1$dbh[alive&inc],splitA,length)
 meandate=tapply(census1$date[alive&inc],splitA,mean,na.rm=TRUE)
 meandbh=tapply(census1$dbh[alive&inc],splitA,mean,na.rm=TRUE) 
 
 N=fill.dimension(N,class1,class2,fill=0)
 meandate=fill.dimension(meandate,class1,class2,fill=NA)
 meandbh=fill.dimension(meandbh,class1,class2,fill=NA)
 
 result=list(N=N,meandbh=meandbh,meandate=meandate)
 gc()
 return(result)
}

# Finds abundance for each species per user defined area, using abundance().  The user
# defined area are conventionally 1 hectare (gridsize=100) or 1 quadrate
# (gridsize=20).

# This returns the list of arrays of abundance().  The first component of
# this list is a matrix of species per user defined area.

abundance.quad=function(census1,mindbh=10,plotdim=c(1000,500),gridsize=100)
{
 sp=census1$sp
 quadno=gxgy.to.index(census1$gx,census1$gy,gridsize=gridsize,plotdim=plotdim)
  
 result=abundance(census1,mindbh=mindbh,split1=sp,split2=quadno)
 gc()
 return(result)
}


# A cover function to calculate total abundance for each species in given
# dbh categories.

abundance.spp=function(census1,alivecode=c("A"),dbhclass=c(10,100,300))
{
# allbreaks=c(dbhbreaks,10000)
# dbhclass=cut(census1$dbh,breaks=allbreaks,right=F,
#              labels=paste(dbhbreaks,allbreaks[-1],sep="."))
 dbh1.vct = sep.dbh(census1,dbhclass=dbhclass)		
 sp=census1$sp

 abund1=abundance(census1,alivecode=alivecode,split1=sp,split2=dbh1.vct)
 gc()
 return(abund1)
}                 


# Finds abundance in two censuses and population change between them.
# Accepts two dataframes, each a complete dataset for a single census.
# Works like abundance with a qualification about dates: mean census date
# is computed from the dates for all individuals alive in either census.

# Plants recruited between the two censuses have a first census date equal to
# the date on which the quadrat they later appear in was censused in the first
# census.

# cannot use dbh as the splitting vector in this function since the dbh of
# the tree has to be examined in both censuses to compute changes in
# abundance by dbh class correctly.  Use abundance.change.dbh to accomplish
# this task properly.

abundance.change=function(census1,census2,alivecode=c("A"),
	mindbh=10,split1=NULL,split2=NULL)
{
 if(is.null(split1)) split1=rep("all",dim(census1)[1])
 if(is.null(split2)) split2=rep("all",dim(census2)[1])

 inc1=inc2=incboth=rep(TRUE,length(census1$dbh))
 
 if(!is.null(mindbh))
  {
   inc1[census1$dbh<mindbh | is.na(census1$dbh)]=FALSE 
   inc2[census2$dbh<mindbh | is.na(census2$dbh)]=FALSE   
   incboth[(census1$dbh<mindbh & census2$dbh<mindbh) | {is.na(census1$dbh) & is.na(census2$dbh)}]=FALSE
  }

  alive1=alive2=rep(FALSE,dim(census1)[1])
	
 for(i in 1:length(alivecode))
  {
   alive1[census1$status==alivecode[i]]=TRUE
   alive2[census2$status==alivecode[i]]=TRUE
  }
 aliveboth=alive1 | alive2
 
 splitA1=list(split1[inc1&alive1],split2[inc1&alive1])
 splitA2=list(split1[inc2&alive2],split2[inc2&alive2])
 splitAboth=list(split1[incboth&aliveboth],split2[incboth&aliveboth])
 splitS=list(split1[alive1&alive2&inc1],split2[alive1&alive2&inc1])

 timeint=(census2$date-census1$date)/365.25
 
 N0=tapply(census1$dbh[inc1&alive1],splitA1,length)
 N1=tapply(census2$dbh[inc2&alive2],splitA2,length)
 S=tapply(census1$dbh[alive1&alive2&inc1],splitS,length)
 date1=tapply(census1$date[incboth&aliveboth],splitAboth,mean,na.rm=T)
 date2=tapply(census2$date[incboth&aliveboth],splitAboth,mean,na.rm=T)

 meantimeint=(date2-date1)/365.25

 class1=sort(unique(split1))
 class2=sort(unique(split2))
 
 N0=fill.dimension(N0,class1,class2)
 N1=fill.dimension(N1,class1,class2)
 S=fill.dimension(S,class1,class2)
 meantimeint=fill.dimension(meantimeint,class1,class2,fill=NA)
 date1=fill.dimension(date1,class1,class2,fill=NA)
 date2=fill.dimension(date2,class1,class2,fill=NA)
 
 a=abundance.change.calculation(N0=N0,N1=N1,S=S,meantimeint=meantimeint)
 
 rate=fill.dimension(a$rate,class1,class2,fill=NA)
 prob=fill.dimension(a$prob,class1,class2,fill=NA)

 result=list(N0=N0,N1=N1,S=S,R=N1-S,rate=rate,prob=prob,
		 meanyrs=meantimeint,date0=date1,date1=date2)
 gc()
 return(result)
}  

# Calculates population change between two censuses in a series of dbh classes.
# Only one additional splitting variable is allowed.
# a rate of change in abundance is computed, but not a probability of this
# value being different from 0 (no change in population size).  This is
# because the loss of a tree from a dbh class can occur in 2 ways, it can
# die or it can grow out of that dbh class and into the next.  The growth
# of a tree is not taken into account by the probability computed for the
# change in abundance when dbh is not one of the split variables.

abundance.change.dbh=function(census1,census2,alivecode=c("A"),dbhclass=c(10,100,300),
                        split=NULL)
{
 if(is.null(split)) split=rep("all",dim(census1)[1])
 dbh1.vct = sep.dbh(census1,dbhclass=dbhclass)		
# dbhclass=cut(census1$dbh,breaks=classbreak,right=FALSE)
 abund1=abundance(census1,alivecode=alivecode,split1=dbh1.vct,split2=split)
 
 dbh2.vct = sep.dbh(census2,dbhclass=dbhclass)
# dbhclass=cut(census2$dbh,breaks=classbreak,right=FALSE)
 abund2=abundance(census2,alivecode=alivecode,split1=dbh2.vct,split2=split)

 meantimeint=(abund2$meandate-abund1$meandate)/365.25
 rate = (log(abund2$N)-log(abund1$N))/meantimeint
 
 class1=sort(unique(dbh1.vct))
 class2=sort(unique(split))
 
 meantimeint=fill.dimension(meantimeint,class1,class2,fill=NA)
 rate=fill.dimension(rate,class1,class2,fill=NA)
 
 return(list(N0=abund1$N,N1=abund2$N,rate=rate,
				 meanyrs=meantimeint,date0=abund1$meandate,date1=abund2$meandate))
}

# functions to calculate change in abundance and test of significance for
# the hypothesis that the number of individuals has not changed over the
# interval measured.  this works on n dimensional arrays

abundance.change.calculation=function(N0=N0,N1=N1,S=S,meantimeint=meantimeint)
{

 N0 <- as.array(N0)
 N1 <- as.array(N1)
 S <- as.array(S)
 meantimeint <-as.array(meantimeint)
 
 rate=(log(N1)-log(N0))/meantimeint
 R = N1 - S
 recrate = R / N1
 survrate = S / N0
 
 change.prob <- mapply(abund.ch.prob,N0=N0,N1=N1,recrate=recrate,survrate=survrate)
 change.prob <- array(data=change.prob,dim=dim(N1),dimnames=dimnames(N1))
 result=list(rate=rate,prob=change.prob)
 return(result)
}

		
abund.ch.prob=function(N0,N1,recrate,survrate)
{

# recrate = R / N1 where R is the number of recruits
		# R = N1 - S
# R = number of expected recruits (rate * initial N)
		#  R is labeled as pi in the descriptive paper, but given that pi has a
		#  fixed meaning in software, I have substituted a different variable
		#  label
# survrate = S / N0 where S is the number of survivors, 
		#  S = N0 - D 
# This function works on a vector of any length.  It does not work on a 2
# dimensional or higher array.  To use abundance.change.calculation or use
# mapply(abund.ch.prob...)  as shown in abundance.change.calculation.
		
 R=recrate*N0
 
 part1=N1[N1<=N0] # population decrease
 part2=N1[N1>N0]  # population increase
 
 prob=numeric()
 i=1
 
 for(x in part1)
  {
   n1=0:x
   
   probvect=dpois(n1,lambda=R[i],log=T)+dbinom(x-n1,N0[i],prob=survrate[i],log=T)
   prob[i]=sum(exp(probvect))

   i=i+1
  }
  
 for(x in part2)
  {
#   k=x-N0
   n1=(x-N0):x

   probvect=dpois(n1,lambda=R[i],log=T)+dbinom(x-n1,N0[i],prob=survrate[i],log=T)
   prob[i]=sum(exp(probvect))   
   
   i=i+1
  }

 return(prob)
}


# Finds mean census interval per species per quadrat, given two fulldatasets.
# It will consider all alive in the first census, the second, or all alive in
# either census, depending on the last argument

meaninterval.quad=function(census1,census2,mindbh=10,plotdim=c(1000,500),
              gridsize=100,xrange=c(0,plotdim[1]),yrange=c(0,plotdim[2]),
              whichcensus="both")
{
 goodloc=census1$gx>=xrange[1] & census1$gy>=yrange[1] &
         census1$gx<xrange[2] & census1$gy<yrange[2] &
         census2$gx>=xrange[1] & census2$gy>=yrange[1] &
         census2$gx<xrange[2] & census2$gy<yrange[2]
 gooddate=census1$date>0 & census2$date>0
                            
 if(whichcensus=="first") inc=(census1$dbh>=mindbh & goodloc & gooddate)
 else if(whichcensus=="second") inc=(census2$dbh>=mindbh & goodloc & gooddate)
 else if(whichcensus=="both") inc=((census1$dbh>=mindbh | census2$dbh>=mindbh) & 
                                    goodloc & gooddate)
 
 census1=subset(census1,inc)
 census2=subset(census2,inc)

 interval=(census2$date-census1$date)/365.25
 
 quadno=gxgy.to.index(census1$gx,census1$gy,gridsize=gridsize,plotdim=plotdim)
 
 meaninterval=tapply(interval,list(census1$sp,quadno),mean)
 
 return(meantime=meaninterval)
}

