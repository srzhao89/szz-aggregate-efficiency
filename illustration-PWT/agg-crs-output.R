##########################################################################################
## Title: Inference for Aggregate Efficiency: Theory and Guidelines for Practitioners
## Authors: Leopold Simar, Valentin Zelenyuk, and Shirong Zhao
## Date: August 7, 2023
## The programming codes used in this paper involve some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
##########################################################################################

require(readxl)
require(FEAR)
source("./Functions/coverage.agg.crs.output.illu.R")

if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)

################################################
np=2
nq=1
################################################################################
pwt100 <- read_excel("Data/pwt100.xlsx", sheet = "Data")
################################################################################
bhz<-c("Albania","Argentina","Armenia","Australia", "Austria", "Azerbaijan","Belarus",
       "Belgium", "Bolivia (Plurinational State of)","Brazil","Bulgaria","Canada","Chile",
       "China","Colombia","Costa Rica","Croatia","Czech Republic","Denmark","Dominican Republic",
       "Ecuador","Estonia","Finland","France","Germany","Greece","Guatemala","Honduras",
       "China, Hong Kong SAR","Hungary","Iceland","India","Indonesia" ,"Ireland","Israel","Italy",
       "Jamaica","Japan","Kazakhstan","Kenya","Republic of Korea","Kyrgyzstan","Latvia", 
       "Lithuania","North Macedonia","Madagascar","Malawi","Malaysia","Mauritius","Mexico",
       "Republic of Moldova","Morocco","Netherlands","New Zealand","Nigeria","Norway","Panama",
       "Paraguay","Peru","Philippines","Poland","Portugal","Romania","Russian Federation","Sierra Leone" ,
       "Singapore","Slovakia","Slovenia","Spain","Sri Lanka","Sweden","Switzerland","Syrian Arab Republic",
       "Taiwan", "Tajikistan","Thailand","Turkey","Ukraine","United Kingdom","Uruguay","United States",
       "Venezuela (Bolivarian Republic of)","Zambia","Zimbabwe")
ii=which(pwt100$country %in% bhz)
df<-pwt100[ii,]
year=seq(1990,2019,1)
df<-df[df$year>=min(year),]
y=matrix(df$cgdpo,nrow=1)
x=t(cbind(df$emp,df$cn))
py=df$cgdpo

res=matrix(NA,nrow=length(year),ncol=12)
res0=matrix(NA,nrow=length(year),ncol=12)
for (i in 1:length(year)) {
  
  ii=which(df$year==i+min(year)-1)
  yi=matrix(y[,ii],nrow=1)
  xi=x[,ii]
  pyi=py[ii]
  
  tt=coverage.agg.crs.illu(x=xi,y=yi,py=pyi,L=100)
  # without data sharpening
  res0[i,1:3]=tt$estimate0
  res0[i,4:6]=tt$sig0
  res0[i,7:8]=tt$bounds10[2,]
  res0[i,9:10]=tt$bounds20[2,]
  res0[i,11:12]=tt$bounds30[2,]
  # with data sharpening
  res[i,1:3]=tt$estimate
  res[i,4:6]=tt$sig
  res[i,7:8]=tt$bounds1[2,]
  res[i,9:10]=tt$bounds2[2,]
  res[i,11:12]=tt$bounds3[2,]
}

round(res0,3)
round(res,3)


### construct the Table Without Data Sharpening ###
tex=formatC(year,width=6,digits=0,format="f")
for (k in 1:ncol(res0)) {
  tex = paste(tex,"&",formatC(res0[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")
write(tex,file="./Output/pwt-coverage-agg-100kappa.tex")

### construct the Table With Data Sharpening ###
tex=formatC(year,width=6,digits=0,format="f")
for (k in 1:ncol(res)) {
  tex = paste(tex,"&",formatC(res[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")
write(tex,file="./Output/pwt-coverage-agg-sharpening-100kappa.tex")

## Plot the estimate of the standard deviation
pdf(file = "./Figures/pwt-sigma-agg-100kappa.pdf",width=10, height=8)

sol1=res0[,4]
sol2=res0[,5]
sol3=res0[,6]

sol4=res[,4]
sol5=res[,5]
sol6=res[,6]

y.min=min(sol1,sol2,sol3,sol4,sol5,sol6)
y.max=max(sol1,sol2,sol3,sol4,sol5,sol6)
plot(year,sol1, type="l",lty=1,lwd=2,col = 1,
     ylim=c(y.min-0.2,y.max),
     ylab=expression(paste("The estimates of standard deviations")),
     xlab="year")
lines(year,sol2, type="l",lty=2,lwd=2,col = 2)
lines(year,sol3, type="l",lty=3,lwd=2,col = 3)
lines(year,sol4, type="l",lty=4,lwd=2,col = 4)
lines(year,sol5, type="l",lty=5,lwd=2,col = 5)
lines(year,sol6, type="l",lty=6,lwd=2,col = 6)
legend('bottom',legend = c("Sol1","Sol2","Sol3","Sol4","Sol5","Sol6"), col = c(1,2,3,4,5,6), 
       lty=c(1,2,3,4,5,6),lwd =c(2,2,2,2,2,2),xpd = TRUE,horiz = TRUE, cex = 0.9, bty = 'n')

dev.off()

## Plot the Figures for the dynamics of the efficiency and its CIs.
require(plotrix)

pdf(file = "./Figures/pwt-CI-agg-100kappa.pdf",width=10, height=8)

eff.agg0=res0[,3]
U0=res0[,12]
L0=res0[,11]

eff.agg=res[,3]
U=res[,12]
L=res[,11]
y.min=min(L0,L)
y.max=max(U0,U)
plotCI(year, eff.agg, ui=U, li=L, 
       col=1,
       lwd=2,
       pch=1,
       ylim=c(y.min-0.1,y.max),
       ylab="The aggregate inefficiency estimates and CIs")
plotCI(year, eff.agg0, ui=U0, li=L0, 
       col=2,
       lwd=2,
       pch=2,
       add = TRUE,
       )
legend('bottom',legend = c("Sol6","Sol3"), col = c(1,2),
       lwd=c(2,2),
       pch=c(1,2),
       xpd = TRUE,horiz = TRUE, cex = 1, 
       bty = 'n')
dev.off()


## Plot the estimate of the aggregate efficiency
pdf(file = "./Figures/pwt-meanbias-agg-100kappa.pdf",width=10, height=8)

case1=res0[,1]
case2=res0[,3]
case3=res[,1]
case4=res[,3]

y.min=min(case1,case2,case3,case4)
y.max=max(case1,case2,case3,case4)
plot(year,case1, type="l",lty=1,lwd=2,col = 1,
     ylim=c(y.min-0.1,y.max),
     ylab=expression(paste("The aggregate inefficiency estimates")),
     xlab="year")
lines(year,case2, type="l",lty=2,lwd=2,col = 2)
lines(year,case3, type="l",lty=3,lwd=2,col = 3)
lines(year,case4, type="l",lty=4,lwd=2,col = 4)
legend('bottom',legend = c("Standard","Bias-Corrected","Standard+Sharpened","Bias-Corrected+Sharpened"), 
       col = c(1,2,3,4), 
       lty=c(1,2,3,4),lwd =c(2,2,2,2),ncol=2,cex = 1, bty = 'n')

dev.off()

mean(case2)
mean(case4)