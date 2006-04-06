# VERSION INFORMATION
# authors: Rick Condit, Pamela Hall
# create date: 2003 CTFS workshop and later
# last modified date: 5/2005

# Revision information
# Revised July 7, 2004 by PH
# functions for finding guadrates, gx,gy location, for neighborhood analyses,
# list of functions included:

# quad functions
# rowcol.to.index=function(rowno,colno,gridsize=20,plotdim=c(1000,500))
# gxgy.to.index=function(gx,gy,gridsize=20,plotdim=c(1000,500))
# index.to.rowcol=function(index,gridsize=20,plotdim=c(1000,500))
# index.to.gxgy=function(index,gridsize=20,plotdim=c(1000,500))
# gxgy.to.rowcol=function(gx,gy,gridsize=20,plotdim=c(1000,500))
# gxgy.to.hectindex=function(gx,gy,plotdim=c(1000,500))
# findborderquads=function(index,dist=20,gridsize=20,plotdim=c(1000,500))

# date functins
# tojulian()
# fromjulian()

# confidence intervals
# find.climits() 

# filling arrays to maintain all categories
# fill.dimension() 
# fill.1dimension()
# fill.distclasses()

# misc
# unidentified.species()
# unwind.matrix()
# rndown5()

# distance between pairs of points
# xydist()
# xydistvect()

#############
rowcol.to.index=function(rowno,colno,gridsize=20,plotdim=c(1000,500))
{
 badrc=(rowno<=0 | colno<=0 | rowno>plotdim[2]/gridsize | colno>plotdim[1]/gridsize)

 rowno=rowno-1
 colno=colno-1
 maxrow=floor(plotdim[2]/gridsize)
 index=colno*maxrow+rowno+1
 if(length(badrc[badrc>0])) index[badrc]=NA
 return(index)
}
#############

gxgy.to.index=function(gx,gy,gridsize=20,plotdim=c(1000,500))
{
 badgxgy=(gx<0 | gy<0 | gx>=plotdim[1] | gy>=plotdim[2] | 
          is.na(gx) | is.na(gy))

 colno=1+floor(gx/gridsize)
 rowno=1+floor(gy/gridsize)
 if(length(badgxgy[badgxgy>0])) colno[badgxgy]=rowno[badgxgy]=NA

 return(rowcol.to.index(rowno,colno,gridsize,plotdim))
}
##############

index.to.rowcol=function(index,gridsize=20,plotdim=c(1000,500))
{
 index=index-1

 badindex=(index<0 | index>=plotdim[1]*plotdim[2]/(gridsize^2))

 maxrow=floor(plotdim[2]/gridsize)
 rowno=index%%maxrow
 colno=floor((index-rowno)/maxrow)
 row=rowno+1
 col=colno+1

 if(length(badindex[badindex>0])) row[badindex]=col[badindex]=NA

 return(data.frame( row=row, col=col))
}
##############

index.to.gxgy=function(index,gridsize=20,plotdim=c(1000,500))
{
 badindex=(index<=0 | index>plotdim[1]*plotdim[2]/(gridsize^2))

 rc=index.to.rowcol(index,gridsize,plotdim)
 gx=gridsize*(rc$col-1)
 gy=gridsize*(rc$row-1)
 
 if(length(badindex[badindex>0])) gx[badindex]=gy[badindex]=NA

 return(data.frame(gx=gx,gy=gy))
}
##############

gxgy.to.rowcol=function(gx,gy,gridsize=20,plotdim=c(1000,500))
{
 index=gxgy.to.index(gx,gy,gridsize,plotdim)
 return(index.to.rowcol(index,gridsize,plotdim))
}
##############

gxgy.to.hectindex=function(gx,gy,plotdim=c(1000,500))
{
	
 badgxgy=(gx<0 | gy<0 | gx>=plotdim[1] | gy>=plotdim[2] | 
          is.na(gx) | is.na(gy))
 ha.rowno=floor(gy/100)
 ha.colno=floor(gx/100)
 max.ha.row=plotdim[2]/100
 if(length(badgxgy[badgxgy>0])) ha.colno[badgxgy]=ha.rowno[badgxgy]=NA
 
 return(ha.colno*max.ha.row+ha.rowno+1)
}

##################################
# Returns indices of neighboring quadrats, for a given quadrat index.
# It includes all neighboring quadrats out to a distance dist. 

findborderquads=function(index,dist=20,gridsize=20,plotdim=c(1000,500))
{
 bound.index=numeric(8)
 no.boundaries=0

 no.row=index.to.rowcol(index,gridsize,plotdim)$row
 no.col=index.to.rowcol(index,gridsize,plotdim)$col
 maxrow=plotdim[2]/gridsize
 maxcol=plotdim[1]/gridsize

 layers=1+floor(dist/gridsize)
 
 for(i in (no.row-layers):(no.row+layers))
  for(j in (no.col-layers):(no.col+layers))
   if(i!=no.row | j!=no.col)
     if(i>=1 & i<=maxrow & j>=1 & j<=maxcol)
        {
         no.boundaries=no.boundaries+1
         bound.index[no.boundaries]=rowcol.to.index(i,j,gridsize,plotdim)
        }

 return( bound.index[bound.index>0] )
} 

##############################################################
# This function starts with a character representation of a date in the form of numbers for month, day year and returns
# a julian. For instance: x="7/3/1988" and returns 3Jul88

tojulian=function(x,dateform="%m/%d/%Y",centurycut=2020)
{
 x=strptime(x,format=dateform)
 
 year=x$year+1900
 year[year>centurycut & !is.na(year)]=year[year>centurycut & !is.na(year)]-100

 month=x$mon+1
 day=x$mday
 return(mdy.date(month=month,day=day,year=year))
}

##############################################################
# This function starts with a julian and returns any character representation of a date.

fromjulian=function(j,dateform="%m/%d/%Y")
{
 s=date.mmddyyyy(j)
 d=strptime(s,format="%m/%d/%Y")
 return(strftime(d,dateform))
}


##############################################################
# This calculates confidence limits around a number of survivors S, out of N individuals.
# The calculation is simple with good Bayesian understanding, using the beta
# distribution.
# Modification of Ricks version so that it returns CI about survivorship.
# Then the computation of the CI can be done as desired for any given
# fucntion call.

find.climits=function(N,S,alpha=.05,kind="upper")
{
 if(kind=="lower")
  {
   result=N*(1-qbeta(1-alpha/2,shape1=N-S+1,shape2=S+1))
   result[S==0]=0
  }
 else if(kind=="upper")
  {
   result=N*(1-qbeta(alpha/2,shape1=N-S+1,shape2=S+1))
   result[S==N]=N[S==N]
  }
  
 return(result)
 
}

##############################################################
# This calculates confidence limits around a number of deaths, D, out of N individuals.
# The calculation is simple with good Bayesian understanding, using the beta
# distribution.

# to use this properly, different S (number of survivors) has to be given
# for CI about recruitment AND the computation of the CI has to be
# calculated differently.  I find this confusing and prefer to have
# identical code for both recruitment and mortality.  I have modified this
# function and recruitment.calculation and mortality.calculation to make
# them consistent.
#find.climits.Rick=function(N,D,alpha=.05,kind="upper")
#{
# if(kind=="lower")
#  {
#   result=N*(1-qbeta(1-alpha/2,shape1=N-D+1,shape2=D+1))
#   result[D==0]=0
#  }
# else if(kind=="upper")
#  {
#   result=N*(1-qbeta(alpha/2,shape1=N-D+1,shape2=D+1))
#   result[D==N]=N[D==N]
#  }
#  
# return(result)
# 
#}

##############################################################
# find.climits.old and find.climits.smallN.old have been supersceded by a new
# version of find.climits.  These older versions are retained in the source
# code for reference.
# This calculates confidence limits around a 
# number of deaths, D, out of N individuals. It uses the binomial distribution 
# around D out of N, finding the value L at which the the binomial probability
# to the left of L is below .025 (and then U at which probability .025 is left above. 
# In R, this works simply with the qbinom function.

#find.climits.old=function(N,D,alpha=.05,kind="upper")
#{
# if(kind=="upper") return(qbinom(1-alpha/2,N,D/N))
# else if(kind=="lower") return(qbinom(alpha/2,N,D/N))
#}

##############################################################
# This is the way I used in the original mortality paper for estimating
# confidence limits on number of deaths, when sample size was N<500 or D<5.
# It searches for the p=D/N at which the probability of the observed data
# is just .025 (or .975). This has the advantage of producing limits even
# if D=0 or D=N (the 'standard' method cannot).

#find.climits.smallN.old=function(N,D,alpha=.05,kind="upper",delta=.00005)
#{
# if(!is.array(N)) N=array(N,dim=length(N))
# if(!is.array(D)) D=array(D,dim=length(D))
# best=best.tmp=p=array(-99,dim=dim(N))
# toolow=toohigh=justright=array(F,dim=dim(N))

# if(kind=="upper") { best=lo=D;  hi=N }
# else if(kind=="lower") { best=lo=array(0,dim=dim(N)); hi=D }
# else return("Enter either upper or lower\n")

# if(kind=="lower") good=D>0
# else if(kind=="upper") good=D<N

# if(kind=="upper") { a1=alpha/2-delta;  a2=alpha/2+delta }
# else if(kind=="lower") { a1=1-alpha/2-delta;  a2=1-alpha/2+delta }

# p[good]=pbinom(D[good],N[good],best[good]/N[good])

# while(length(p[(p>a2 | p<a1) & good])>0)
#  {
#   if(kind=="lower") p[good]=pbinom(D[good]-1,N[good],best[good]/N[good])
#   if(kind=="upper") p[good]=pbinom(D[good],N[good],best[good]/N[good])

#   toolow[good]=p[good]>a2
#   toohigh[good]=p[good]<a1

#   justright[good]=p[good]>=a1 & p[good]<=a2

#   best.tmp[toohigh&good]=(best+lo)[toohigh&good]/2
#   hi[toohigh&good]=best[toohigh&good]

#   best.tmp[toolow&good]=(best+hi)[toolow&good]/2
#   lo[toolow&good]=best[toolow&good]

#   best.tmp[justright&good]=best[justright&good]

#   best=best.tmp
#  }

# if(kind=="lower") best[!good]=0
# else if(kind=="upper") best[!good]=N[!good]

# return(best)
#}

##############################################################
# This checks whether the results of tapply in mortality(), growth() and
# recruit() have all dimensions.  Dataarray is the output of tapply, and
# will have dim=2, with colnames=class2 and rownames=class1.  If a dimension is
# missing, it will be filled with the value given by fill.

fill.dimension=function(dataarray,class1,class2,fill=0)
{
 if(dim(dataarray)[1]<length(class1))
  {
   m=match(class1,dimnames(dataarray)[[1]])
   dataarray=array(dataarray[m,],
                   dim=c(length(class1),length(class2)),
                   dimnames=list(class1,class2))
  }

 if(dim(dataarray)[2]<length(class2))
  {
   m=match(class2,dimnames(dataarray)[[2]])
   dataarray=array(dataarray[,m],
                   dim=c(length(class1),length(class2)),
                   dimnames=list(class1,class2))
   colnames(dataarray)=class2
  }
  
  dataarray[is.na(dataarray)]=fill
  return(dataarray)
 }
  
##############################################################
# Filling a 1D array resulting from tapply.

fill.1dimension=function(dataarray,class1,fill=0)
{
 if(dim(dataarray)[1]<length(class1))
  {
   m=match(class1,dimnames(dataarray)[[1]])
   dataarray=dataarray[m]
   
   names(dataarray)=class1
  }

  dataarray[is.na(dataarray)]=fill
  return(dataarray)
 }
  
##############################################################
# Fills a 1D arrary resulting from tapply for distance classes
# This is used by the density dependence funtions: calc.neighbordens.sp
# Ricks comments:  
# This checks whether the results of tapply in countdens has all distclasses.
# Dataarray is the dropped output of tapply, and will have dim=1, with 
# names=distclasses. If a class is missing, it must be filled with a zero.

fill.distclasses=function(dataarray,distclasses)
{
 if(length(dataarray)[1]<length(distclasses))
  {
   m=match(distclasses,names(dataarray))
   dataarray=dataarray[m,]
  }

 dataarray[is.na(dataarray)]=0
 return(dataarray)
}



##############################################################

# Takes a string of species names or codes and returns a logical vector indicating
# with TRUE those that are unidentified.

unidentified.species=function(spnames,unidenstr=c("**","UNID","unid","UNIDEN","uniden"))
{
 remove=numeric()
 unidname=rep(FALSE,length(spnames))
 
 if(is.null(unidenstr)) return(unidname)
 
 for(i in 1:length(unidenstr))
   remove=c(remove,grep(unidenstr[i],spnames,fixed=TRUE))

 unidname[remove]=TRUE
 unidname[is.na(spnames)]=TRUE
 return(unidname)
} 


##############################################################
# unwind.matrix.R

# unwinds a matrix into a vector
# presumes order of unwinding by column

unwind.matrix=function(matrix.in)

{
output=vector()
for (i in 1:dim(matrix.in)[2])
 {
	for (j in 1:dim(matrix.in)[1])
	{
		k=((i-1)*dim(matrix.in)[1])+j
		output[k]<-matrix.in[j,i]
	}
}

return(as.vector(output))

}

#############################################################
# rndown5
rndown5=function(s) 
{
	return(5*floor(s/5))
}

############################################################
# Distance between x-y coordinates

xydist=function(x1,y1,x2,y2) return( sqrt( (x1-x2)^2 + (y1-y2)^2 ) )

#############################################################
# Likewise, but accepts each set of coordinates as a vector of length 2,
# with the first element the x coordinates, the second y. For use with 
# apply.

xydistvect=function(pt1,pt2) 
  return( sqrt( (pt1[1]-pt2[1])^2 + (pt1[2]-pt2[2])^2 ) )


