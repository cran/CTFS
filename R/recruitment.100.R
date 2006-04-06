# VERSION INFORMATION ver 1.0
# authors: Rick Condit
# create date: 2003 CTFS workshop
# modifications: Pam Hall

# functions included:
# recruitment()
# recruitment.eachspp()
# recruitment.calculation()

# Overvie of functions for calculating recruitment rates.

# Recruitment() is the main function, and is constructed like growth() and
# mortality().  It requires two complete datasets, one per census, with
# dbh, pom, and date for every individual of all species in at least 2
# censuses.  The two datasets have exactly the same individuals, in exactly
# the same order, one individual per row.  It can then calculate
# recruitment rates based on up to user-submitted factors.

# Recruitment is based on status and dbh.  Any status indicating a live
# tree can be submitted in the variable alivecode.  Survivors are all
# individuals alive in both censuses, with status==A in the first census,
# and larger than the minimum dbh in the first census.  The total
# population in the second census includes all those alive, above the
# minimum dbh, plus any other survivors.

# As in mortality, individuals whose status is NA in either census are deleted
# from all calculations.

recruitment=function(census1,census2,mindbh=10,alivecode=c("A","AB","AS"),
                     split1=NULL,split2=NULL)
{
 if(is.null(split1)) split1=rep("all",dim(census1)[1])
 if(is.null(split2)) split2=rep("all",dim(census2)[1])

 inc=!is.na(census1$status) & !is.na(census2$status)
 census1=census1[inc,]
 census2=census2[inc,]
 split1=split1[inc]
 split2=split2[inc]
 
 timeint=(census2$date-census1$date)/365.25

 survivor=alive1=alive2=rep(FALSE,dim(census1)[1])
 alive1[census1$status=="A"]=TRUE
 for(i in 1:length(alivecode))
  {
   survivor[census1$status=="A" & census2$status==alivecode[i]]=TRUE
   alive2[census2$status==alivecode[i]]=TRUE
  }
 S.inc = survivor & census1$dbh>=mindbh
 N1.inc = (alive2 & census2$dbh>=mindbh) | S.inc

 class1=sort(unique(split1))
 class2=sort(unique(split2))
 splitS=list(split1[S.inc],split2[S.inc])
 splitN0=list(split1[alive1],split2[alive1]) 
 splitN1=list(split1[N1.inc],split2[N1.inc])
   
 S=tapply(census2$dbh[S.inc],splitS,length)
 N1=tapply(census2$dbh[N1.inc],splitN1,length)
 meantimeint=tapply(timeint[N1.inc],splitN1,mean,na.rm=T)
 startdate=tapply(census1$date[alive1],splitN0,mean,na.rm=T)
 enddate=tapply(census2$date[N1.inc],splitN1,mean,na.rm=T)
 
 S=fill.dimension(S,class1,class2)
 N1=fill.dimension(N1,class1,class2)
 meantimeint=fill.dimension(meantimeint,class1,class2,fill=NA)
 startdate=fill.dimension(startdate,class1,class2,fill=NA)
 enddate=fill.dimension(enddate,class1,class2,fill=NA)

 if(sum(N1)==0) {
   result=list(N1=rep(NA,length(class1)),S=rep(NA,length(class1)),R=rep(NA,length(class1)),
               rate=rep(NA,length(class1)),lowerCI=rep(NA,length(class1)),
               upperCI=rep(NA,length(class1)),meanyrs=rep(NA,length(class1)),
               date0=rep(NA,length(class1)),date1=rep(NA,length(class1)))
	 gc()
	 return(result)
   }

   r=recruitment.calculation(N1=N1,S=S,meantimeint=meantimeint)            
# lower.ci=upper.ci=N2CI
# lower.ci=find.climits(N1,S,kind="lower")
# upper.ci=find.climits(N1,S,kind="upper")

# rec.rate=(log(N1)-log(S))/timeint
# upper.rate=(log(N1)-log(lower.ci))/timeint
# lower.rate=(log(N1)-log(upper.ci))/timeint

# rec.rate[S==0]=upper.rate[S==0]=Inf
# upper.rate[lower.ci==0]=Inf
# rec.rate[N1==0]=lower.rate[N1==0]=upper.rate[N1==0]=NA

 result=list(N1=N1,S=S,R=N1-S,rate=r$rate,lower=r$lowerCI,upper=r$upperCI,
             meanyrs=meantimeint,date0=startdate,date1=enddate)
 gc()
 return(result)
}


# A cover function to calculate recruitment for each species in given dbh
# categories.

recruitment.eachspp=function(census1,census2,mindbh=10,alivecode=c("A","AB","AS"))
{
 sp=census1$sp
 
 result=recruitment(census1,census2,mindbh=mindbh,alivecode=alivecode,split1=sp)
 gc()
 return(result)
}

recruitment.calculation=function(N1,S,meantimeint)
{
 lower.ci=upper.ci=N1
 lower.ci=find.climits(N1,S,kind="lower")
 upper.ci=find.climits(N1,S,kind="upper")

 # find.climits returns the upper and lower values of DEAD trees.
 # recruitment = N1/S  S has to be computed from the upper and lower
 # number of D, hence, 


 rec.rate=(log(N1)-log(S))/meantimeint
 upper.rate=(log(N1)-log(lower.ci))/meantimeint
 lower.rate=(log(N1)-log(upper.ci))/meantimeint

 rec.rate[S==0]=upper.rate[S==0]=Inf
 upper.rate[lower.ci==0]=Inf
 rec.rate[N1==0]=lower.rate[N1==0]=upper.rate[N1==0]=NA

 if(is.null(dim(N1)))
   return(data.frame(N1=N1,S=S,R=N1-S,rate=rec.rate,lowerCI=lower.rate,upperCI=upper.rate,
			 meantimeint=meantimeint))
 else
   return(list(N1=N1,S=S,R=N1-S,rate=rec.rate,lowerCI=lower.rate,upperCI=upper.rate,
			 meantimeint=meantimeint))
}


# noR=(N1==S)
# lower.ci[!noR]=find.climits(N1[!noR],S[!noR],kind="lower")
# lower.ci[noR]=find.climits.smallN(N1[noR],S[noR],kind="lower")

# noS=(S==0)
# upper.ci[!noS]=find.climits(N1[!noS],S[!noS],kind="upper")
# upper.ci[noS]=find.climits.smallN(N1[noS],S[noS],kind="upper")
