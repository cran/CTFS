\name{texttordata}
\alias{texttordata}
\concept{tab delimited text}
\concept{reading data into R}

\title{Creates CTFS R dataset from any tab delimited text file}
\description{
Creates an rdata file from a text file and makes the first column the
rownames, effectively duplicating the information in the first column.
	}
\usage{
texttordata(text.in,missvalue=NA)	}

\arguments{
  \item{text.in}{any text file with defined columns, file name must be
  in quotes. Either change to the directory in which the file is
  located or provide the full path name and file name in quotes.  Use
  \code{getwd} to find out the full path name if necessary to use,
  or use \code{setwd} to change to the directory where the text file
  is located.}
  \item{missvalue}{any value that should be considered missing and
  replaced with NA. \code{missvalue} can be a vector of strings which
  are to be interpreted as NA values.  Blank or empty fields are also
  considered to be missing values.}
	}
\details{
  Rownames are numbers.  Column names are taken from the first row in
  the text file.  This function uses the R function \code{read.table}.
  If the default values of \code{texttordata} are not sufficient,
  read the help pages for \code{read.table}.
   
  \preformatted{
  read.table(text.in,as.is=TRUE,header=TRUE,sep="\t",quote="",comment.char="",na.strings=missvalue)
  }

  \code{as.is} indicates that variable formats should be retained
  \code{header} indicates that the first row of the text file provides
  the column (variable) names
  \code{sep} indicates the separator, default is \code{tab}
	}

\value{
 Returns a \strong{dataframe} with rownames as an additional column. 
}

\author{Rick Condit}
\seealso{\code{\link{CTFS.readwrite}}}
	
\examples{\dontrun{
> texttordata("bcispp.info.txt",missvalue=c("."))->info.test
> str(info.test)
`data.frame':	318 obs. of  8 variables:
 $ sp      : chr  "acacme" "acaldi" "acalma" "ade1tr" ...
 $ genus   : chr  "Acacia" "Acalypha" "Acalypha" "Adelia" ...
 $ species : chr  "melanoceras" "diversifolia" "macrostachya" "triloba" ...
 $ family  : chr  "Fabaceae:Mimos." "Euphorbiaceae" "Euphorbiaceae" "Euphorbiaceae" ...
 $ grform  : chr  "U" "S" "U" "U" ...
 $ repsize : int  4 2 2 10 4 20 16 4 4 20 ...
 $ breedsys: chr  "B" "M" "M" "D" ...
 $ maxht   : int  6 6 5 5 15 15 -1 6 10 15 ...
> info.test[1:10,]
       sp      genus       species          family grform repsize breedsys maxht
1  acacme     Acacia   melanoceras Fabaceae:Mimos.      U       4        B     6
2  acaldi   Acalypha  diversifolia   Euphorbiaceae      S       2        M     6
3  acalma   Acalypha  macrostachya   Euphorbiaceae      U       2        M     5
4  ade1tr     Adelia       triloba   Euphorbiaceae      U      10        D     5
5  aegipa  Aegiphila    panamensis     Verbenaceae      M       4        B    15
6  alchco  Alchornea costaricensis   Euphorbiaceae      T      20        D    15
7  alchla  Alchornea     latifolia   Euphorbiaceae      M      16        D    -1
8  alibed  Alibertia        edulis       Rubiaceae      U       4        D     6
9  allops Allophylus  psilospermus     Sapindaceae      M       4        P    10
10 alsebl     Alseis     blackiana       Rubiaceae      T      20        B    15
\ldots
}
}
\keyword{manip}
