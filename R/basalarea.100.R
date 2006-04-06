# VERSION INFORMATION ver 1.0
# authors: Rick Condit
# create date: 2003 CTFS workshop
# modifications: Pam Hall

# functions included:
# basal.area()
# ba.change()
# ba.indiv()
# convert.ba.mult()
# ba()
# basum()

# Overview of functions for calculating basal area 
# functions for calculating basal area, structured in the same fashion as
# abundance() functions.  Basal area is computed assuming the tree is a
# cylinder.  Multiple stems area added to the calculation using the
# multi-stem file.

# Functions for calculating total basal area based on 1 or 2 grouping variables,
# submitted as vectors, each having exactly the same length as the number of rows in the
# full dataset. The first vector should be the one with the most categories
# (as in the case of species).

# If no grouping variable is desired, the vectors should still be submitted as 
# constants.

basal.area=function(census1,mult1=NULL,alivecode=c("A"),mindbh=10,
                    split1=NULL,split2=NULL)
{

 if(is.null(split1)) split1=rep("all",dim(census1)[1])
 if(is.null(split2)) split2=rep("all",dim(census1)[1])

 #inc=rep(FALSE,length(census1$dbh))
 #if(!is.null(mindbh)) inc[census1$dbh>=mindbh]=TRUE
 #else inc=rep(TRUE,length(census1$dbh))

 inc=rep(TRUE,dim(census1)[1])
 if(!is.null(mindbh)) inc[census1$dbh<mindbh | is.na(census1$dbh)]=FALSE

 alive=rep(FALSE,dim(census1)[1])
 for(i in 1:length(alivecode)) alive[census1$status==alivecode[i]]=TRUE
 
 class1=sort(unique(split1))
 class2=sort(unique(split2))
 splitB=list(split1[inc&alive],split2[inc&alive])

 if(!is.null(mult1))
  {
   bam=convert.ba.mult(census1,mult1,mindbh=mindbh)
   m=match(census1$tag,bam$tag)
   bamall=bam$ba[m]
   bamall[is.na(bamall)]=0
   baind=ba(census1$dbh)+bamall
  }
 else baind=ba(census1$dbh)
 
 batot=tapply(baind[alive&inc],splitB,sum,na.rm=TRUE)
 meandate=tapply(census1$date[alive&inc],splitB,mean,na.rm=TRUE)

 batot=fill.dimension(batot,class1,class2,fill=0)
 meandate=fill.dimension(meandate,class1,class2,fill=NA)

 result=list(basalarea=batot,meandate=meandate)
 gc()
 return(result)
 
}

# Calculates total basal area in two different censuses, and it's rate of change.
# Works exactly like abundance, except for multiple stems, which are stored in separate
# files. See function below convert.ba.mult for how the multiple stems are merged.

ba.change=function(census1,census2,mult1=NULL,mult2=NULL,alivecode=c("A"),mindbh=10,
                   split1=NULL,split2=NULL)
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
 
 splitB1=list(split1[inc1&alive1],split2[inc1&alive1])
 splitB2=list(split1[inc2&alive2],split2[inc2&alive2])
 splitBboth=list(split1[incboth&aliveboth],split2[incboth&aliveboth])

 if(!is.null(mult1))
  {
   bam0=convert.ba.mult(census1,mult1,mindbh=mindbh)
   m=match(census1$tag,bam0$tag)
   bam0all=bam0$ba[m]
   bam0all[is.na(bam0all)]=0
   batot0=ba(census1$dbh)+bam0all
  }
 else batot0=ba(census1$dbh)

 if(!is.null(mult2))
  {
   bam1=convert.ba.mult(census2,mult2,mindbh=mindbh)
   m=match(census2$tag,bam1$tag)
   bam1all=bam1$ba[m]
   bam1all[is.na(bam1all)]=0
   batot1=ba(census2$dbh)+bam1all
  }
 else batot1=ba(census2$dbh)
 
 timeint=(census2$date-census1$date)/365.25
 
 ba0=tapply(batot0[inc1&alive1],splitB1,sum,na.rm=TRUE)
 ba1=tapply(batot1[inc2&alive2],splitB2,sum,na.rm=TRUE)
 date1=tapply(census1$date[incboth&aliveboth],splitBboth,mean,na.rm=TRUE)
 date2=tapply(census2$date[incboth&aliveboth],splitBboth,mean,na.rm=TRUE)

 meantimeint=(date2-date1)/365.25

 class1=sort(unique(split1))
 class2=sort(unique(split2))
 
 ba0=fill.dimension(ba0,class1,class2)
 ba1=fill.dimension(ba1,class1,class2)
 date1=fill.dimension(date1,class1,class2,fill=NA)
 date2=fill.dimension(date2,class1,class2,fill=NA)
 meantimeint=fill.dimension(meantimeint,class1,class2,fill=NA)
 
 ba.rate=(log(ba1)-log(ba0))/meantimeint
 ba.rate=fill.dimension(ba.rate,class1,class2,fill=NA)

 result=list(BA0=ba0,BA1=ba1,rate=ba.rate,
		 meanyrs=meantimeint,date0=date1,date1=date2)
 gc()
 return(result)
}  


# computes the basal area for individual, including multstems dataset if
# available and returns a vector of basal areas of the same length at the input dataset

ba.indiv=function(census1,mult1,alivecode=c("A"),mindbh=10)
{
 inc=rep(TRUE,dim(census1)[1])
 if(!is.null(mindbh)) inc[census1$dbh<mindbh | is.na(census1$dbh)]=FALSE

if(!is.null(mult1))
  {
   bam=convert.ba.mult(census1,mult1,mindbh=mindbh)
   m=match(census1$tag,bam$tag)
   bamall=bam$ba[m]
   bamall[is.na(bamall)]=0
   baind=ba(census1$dbh)+bamall
  }
 else baind=ba(census1$dbh)

 baind[!inc]=NA
 result=as.vector(baind)
 gc()
 return(result)
}

# This takes a multiple stem dataframe, sums the BA of all stems for each
# tag, and returns a dataframe that has one tag per row, with its total BA of multiple stems.

convert.ba.mult=function(census,mult,mindbh=10)
{
	
#inc=rep(FALSE,length(census1$dbh))
#if(!is.null(mindbh)) inc[census1$dbh>=mindbh]=TRUE
#else inc=rep(TRUE,length(census1$dbh))

 inc=rep(TRUE,dim(census)[1])
 if(!is.null(mindbh)) inc[census$dbh<mindbh | is.na(census$dbh)]=FALSE

 ba=tapply(mult$dbh[inc],mult$tag[inc],basum)
 
 m=match(as.numeric(names(ba)),mult$tag)
 result=data.frame(tag=as.numeric(names(ba)),ba)
 
 gc()
 return(result)
}  

# converts a circle's diameter to it's area

ba=function(dbh) return(pi*(dbh/2000)^2)


# converts a circle's diameter into the total (summed) area

basum=function(dbh,mindbh=10)
{
 dbh=dbh[dbh>=mindbh]

 if(length(dbh)==0) return(0)
 return(sum(ba(dbh),na.rm=TRUE))
}
