###########################################################################################
## Title: Inference for Aggregate Efficiency: Theory and Guidelines for Practitioners
## Authors: Leopold Simar, Valentin Zelenyuk, and Shirong Zhao
## Date: August 7, 2023
## The programming codes used in this paper involve some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
##############
## The following codes will report the results for the aggregate efficiency in
## Section 5.1 of the above mentioned paper.
## Here we use VRS-DEA method in the input orientation, and the dimension is 4
## x is np times n matrix for inputs 
## y is nq times n matrix for outputs
## wx is the cost in this rice example
###########################################################################################

require(readxl)
require(FEAR)
source("./Functions/coverage.agg.vrs.input.illu.R")
if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)

###############################################################
np=3
nq=1
################################################################
rice <- read_excel("./Data/rice.xls", 
                   col_types = c("numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric"))

y=matrix(rice$PROD,nrow=1)
x=t(cbind(rice$AREA,rice$LABOR,rice$NPK))
wx=rice$AREA*rice$AREAP+rice$LABOR*rice$LABORP+rice$NPK*rice$NPKP
year=rice$YEARDUM
#
res=matrix(NA,nrow=9,ncol=12)
res0=matrix(NA,nrow=9,ncol=12)
for (i in 1:9) {
  ii=which(year==i)
  if (length(ii)>0){
    yi=matrix(y[,ii],nrow=1)
    xi=x[,ii]
    wxi=wx[ii]
  } else {
    yi=y
    xi=x
    wxi=wx
  }
  tt=coverage.agg.vrs.input.illu(x=xi,y=yi,wx=wxi,L=20)
  # without data sharpening
  res0[i,1:3]=tt$estimate0
  res0[i,4:6]=tt$sig0
  res0[i,7:8]=tt$bounds10[3,]
  res0[i,9:10]=tt$bounds20[3,]
  res0[i,11:12]=tt$bounds30[3,]
  # with data sharpening
  res[i,1:3]=tt$estimate
  res[i,4:6]=tt$sig
  res[i,7:8]=tt$bounds1[3,]
  res[i,9:10]=tt$bounds2[3,]
  res[i,11:12]=tt$bounds3[3,]
}

round(res0,3)
round(res,3)
year=seq(1:8)+1989

### construct the Table Without Data Sharpening ###
tex=formatC(year,width=6,digits=0,format="f")
tex=append(tex,"Pooled")
for (k in 1:ncol(res0)) {
  tex = paste(tex,"&",formatC(res0[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")
write(tex,file="./Output/coverage-rice.tex")

### construct the Table With Data Sharpening ###
tex=formatC(year,width=6,digits=0,format="f")
tex=append(tex,"Pooled")
for (k in 1:ncol(res)) {
  tex = paste(tex,"&",formatC(res[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")
write(tex,file="./Output/coverage-rice-sharpening.tex")

