###########################################################################################
## Title: Inference for Aggregate Efficiency: Theory and Guidelines for Practitioners
## Authors: Leopold Simar, Valentin Zelenyuk, and Shirong Zhao
## Date: August 7, 2023
## The programming codes used in this paper involve some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
##############
## The following codes will report the results for the simple mean efficiency in
## Section 5.2 of the above mentioned paper.
## Here we use CRS-DEA method in the output orientation, and the dimension is 3
## x is np times n matrix for inputs 
## y is nq times n matrix for outputs
###########################################################################################

coverage.simple.crs.illu <- function(x,y,L=20) {
  np=nrow(x)
  nq=nrow(y)
  n=ncol(x)
  na=floor(n/2)
  nb=n-na
  kappa=2/(np+nq)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  tau0=n^(-kappa)
  # evaluate efficiency using CRS-DEA, output Direction
  d0=FEAR::dea(XOBS=x,YOBS=y,XREF=x,YREF=y, 
            METRIC=1,ORIENTATION=2, RTS=3)
  # using data sharpening method in NSZ2022
  ii=which(d0>=1-tau0)
  epison=runif(length(ii),1-tau0,1)
  d=d0
  d[ii]=d0[ii]*epison
  d0=1/d0
  d=1/d
  epison=1/epison
  # compute bias corrections via generalized jackknife:
  tbar=rep(0,n)
  tbar0=rep(0,n)
  ind=c(1:n)
  for (j in 1:L) {
    if (j==1) {
      ind1=c(1:n)
      x.b=x
      y.b=y
    } else {
      ind1=sample(ind,size=n)
      x.b[,1:n]=x[,ind1]
      y.b[,1:n]=y[,ind1]
    }
    #
    da0=FEAR::dea(XOBS=matrix(x.b[,1:na],nrow=np),
                  YOBS=matrix(y.b[,1:na],nrow=nq),
                  XREF=matrix(x.b[,1:na],nrow=np),
                  YREF=matrix(y.b[,1:na],nrow=nq), 
                  METRIC=1,ORIENTATION=2, RTS=3)
    da0=1/da0
    db0=FEAR::dea(XOBS=matrix(x.b[,(na+1):n],nrow=np),
                  YOBS=matrix(y.b[,(na+1):n],nrow=nq),
                  XREF=matrix(x.b[,(na+1):n],nrow=np),
                  YREF=matrix(y.b[,(na+1):n],nrow=nq), 
                  METRIC=1,ORIENTATION=2, RTS=3)
    db0=1/db0
    #
    tbar0[ind1[1:na]]=tbar0[ind1[1:na]] +
      da0 - d0[ind1[1:na]]
    tbar0[ind1[(na+1):n]]=tbar0[ind1[(na+1):n]] +
      db0 - d0[ind1[(na+1):n]]
  }
  #
  # tbar/tbar0 contains the bias for eff of each obs (with and without data sharpening)
  tbar0=(1/L)*bc.fac*tbar0
  tbar=tbar0
  tbar[ii]=tbar0[ii]*epison
  # mean bias
  mu.bias=mean(tbar) 
  mu.bias0=mean(tbar0) 
  # d.bc/d.bc0 is bias-corrected eff for each obs
  d.bc=d-tbar
  d.bc0=d0-tbar0
  # mean eff and bias-corrected mean eff
  mu=mean(d)
  mu.bc=mu-mu.bias
  mu0=mean(d0)
  mu.bc0=mu0-mu.bias0
  # variance
  # KSW2015 original method
  var=var(d) 
  # SZZ2023 method
  var.bc=var(d.bc) 
  # SZ2020 method
  var.bc.sz=var+mu.bias^2 
  #
  # KSW2015 original method
  var0=var(d0) 
  # SZZ2023 method
  var.bc0=var(d.bc0) 
  # SZ2020 method
  var.bc.sz0=var0+mu.bias0^2 
  # KSW2015 original method
  sig1=sqrt(var)
  sig10=sqrt(var0)
  # SZ2020 correction method
  sig2=sqrt(var.bc.sz)
  sig20=sqrt(var.bc.sz0)
  # SZZ2023 method
  sig3=sqrt(var.bc)
  sig30=sqrt(var.bc0)
  # for non-data sharpening
  ts10=sig10/sqrt(n)
  ts20=sig20/sqrt(n)
  ts30=sig30/sqrt(n)
  # for data sharpening
  ts1=sig1/sqrt(n)
  ts2=sig2/sqrt(n)
  ts3=sig3/sqrt(n)
  #
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  #
  bounds1=matrix((mu.bc-ts1*crit),nrow=3,ncol=2)
  bounds2=matrix((mu.bc-ts2*crit),nrow=3,ncol=2)
  bounds3=matrix((mu.bc-ts3*crit),nrow=3,ncol=2)
  #
  bounds10=matrix((mu.bc0-ts10*crit),nrow=3,ncol=2)
  bounds20=matrix((mu.bc0-ts20*crit),nrow=3,ncol=2)
  bounds30=matrix((mu.bc0-ts30*crit),nrow=3,ncol=2)
  # make a list of results to return to calling routine and then quit:
  res=list(bounds1=bounds1,bounds2=bounds2,bounds3=bounds3,
           bounds10=bounds10,bounds20=bounds20,bounds30=bounds30,
           sig=c(sig1,sig2,sig3),
           sig0=c(sig10,sig20,sig30),
           estimate=c(mu,mu.bias,mu-mu.bias),
           estimate0=c(mu0,mu.bias0,mu0-mu.bias0)
  )
  return(res)
}
