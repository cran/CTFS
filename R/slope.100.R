# slope.092.R

# meantopo.quad 
# quadslope 
# calcslope 
# depends upon: findborderquads(), rowcol.to.index() both in utilities.R

# functions for calculating slope of all quadrats in a plot calcslope takes
# 3 elevations and finds the slope of the plane through them quadslope
# divides the 4 corners of a quadrat into 4 different groups of 3 stakes,
# takes the slope of each, then averages these were first written in C++
# see slopeelev.cpp for more on the geometry allquadratslopes goes through
# all 20x20 m quadrats in a plot and finds the slope, mean elevation, and
# convexity of each convexity is the mean elevation of one 20x20 m quadrat
# relative (minus) the mean of its immediate neighbors

# Helene added a section to correct convexity in edge quadrats

meantopo.quad=function(elevfile,gridsize=20,plotdim=c(1000,500))
{
 rw=cl=0
 on.exit(cat(rw," ",cl))
 
 columns=1+max(elevfile$col$x)/gridsize
 rows=1+max(elevfile$col$y)/gridsize
 totalquads=(columns-1)*(rows-1)

 elevdata=elevfile$col[elevfile$col$x%%gridsize==0 & elevfile$col$y%%gridsize==0,]
 elevmat=matrix(elevdata$elev,nrow=rows,ncol=columns,byrow=F)

 meanelev=convex=convex2=slope=numeric()
 corner=sideht=numeric()

 for(c in 1:(columns-1)) 
  for(r in 1:(rows-1))
   {
    quad.index=rowcol.to.index(r,c,gridsize=gridsize,plotdim=plotdim)

    corner[1]=elevmat[r,c]
    corner[2]=elevmat[r+1,c]
    corner[3]=elevmat[r+1,c+1]
    corner[4]=elevmat[r,c+1]

    meanelev[quad.index]=mean(corner)
    slope[quad.index]=quadslope(corner,gridsize=gridsize)[1]
   }


 for(i in 1:totalquads)
   {
    neighbor.quads=findborderquads(i,gridsize=gridsize,plotdim=plotdim)
    meanelev.neighbor=mean(meanelev[neighbor.quads])
    convex[i]=meanelev[i]-meanelev.neighbor
   }

# correcting convexity in edge quadrats, based on center of the 20x20 rather
# than surrounding 20x20s

 for(cl in 1:(columns-1))
   for(rw in 1:(rows-1))
   {
    if((cl==1) | (cl==(columns-1)) | (rw==1) | (rw==(rows-1))) 
    { 
     quad.index=rowcol.to.index(rw,cl,gridsize=gridsize,plotdim=plotdim)
     convex[quad.index]=elevfile$mat[4*rw-1,4*cl-1]-meanelev[quad.index]
    }
   }

 return(data.frame(meanelev=meanelev,convex=convex,slope=slope))
}



quadslope=function(cornerelev,gridsize=20)
{
 slope=numeric(4)
 z=numeric(3)

 for(j in 1:4)
  {
   post=1

   for(k in (j+1):(j+3))
      {
       if(k>4) m=k%%4
       else m=k

       z[post]=cornerelev[m]
       post=post+1
      }

   slope[j]=calcslope(z,gridsize)
  }

 return( c(mean(slope),sqrt(var(slope))) )
}



calcslope=function(z, gridsize=20)
{
 z2=z[3]-z[2]
 z1=z[1]-z[2]

 if(z2==0)
  {
   if(z1==0) return(0)
   else denom=sqrt( 1+(gridsize/z1)^2 )

   theta1 = acos( (-gridsize/z1)/denom )
   theta2 = acos( (gridsize/z1)/denom )
  }
 
 else 
  {
   denom = sqrt( 1+(z1/z2)^2+(gridsize/z2)^2 )

   theta1 = acos( (-gridsize/z2)/denom )
   theta2 = acos( (gridsize/z2)/denom )
  }

 if(theta1<=theta2) return(180*theta1/pi)
 else return(180*theta2/pi)
}
