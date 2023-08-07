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

coverage.agg.vrs.input.illu <- function(x,y,wx,L=20) {
  
  np=nrow(x)
  nq=nrow(y)
  n=ncol(x)
  na=floor(n/2)
  nb=n-na
  kappa=2/(np+nq+1)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  tau0=n^(-kappa)
  # evaluate efficiency using VRS-DEA, input Direction
  d0= FEAR::dea(XOBS=x,YOBS=y,XREF=x,YREF=y, 
                 METRIC=2,ORIENTATION=1, RTS=1)
  # using data sharpening method in NSZ2022
  ii=which(d0>=1-tau0)
  epison=runif(length(ii),1-tau0,1)
  d=d0
  d[ii]=d0[ii]*epison
  # compute bias corrections via generalized jackknife:
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
    ###
    da0=FEAR::dea(XOBS=matrix(x.b[,1:na],nrow=np),
                  YOBS=matrix(y.b[,1:na],nrow=nq),
                  XREF=matrix(x.b[,1:na],nrow=np),
                  YREF=matrix(y.b[,1:na],nrow=nq), 
                  METRIC=2,ORIENTATION=1, RTS=1)
    db0=FEAR::dea(XOBS=matrix(x.b[,(na+1):n],nrow=np),
                  YOBS=matrix(y.b[,(na+1):n],nrow=nq),
                  XREF=matrix(x.b[,(na+1):n],nrow=np),
                  YREF=matrix(y.b[,(na+1):n],nrow=nq), 
                  METRIC=2,ORIENTATION=1, RTS=1)
    #
    tbar0[ind1[1:na]]=tbar0[ind1[1:na]] +
      da0 - d0[ind1[1:na]]
    tbar0[ind1[(na+1):n]]=tbar0[ind1[(na+1):n]] +
      db0 - d0[ind1[(na+1):n]]
  }
  #
  z2=wx
  z1=d*z2
  z10=d0*z2
  # tbar/tbar0 contains the bias for eff of each obs (with and without data sharpening)
  tbar0=(1/L)*bc.fac*tbar0
  tbar=tbar0
  tbar[ii]=tbar0[ii]*epison
  # mean bias
  mu1.bias=mean(tbar*z2)
  mu.bias=mean(tbar*z2)/mean(z2) 
  mu1.bias0=mean(tbar0*z2)
  mu.bias0=mean(tbar0*z2)/mean(z2)
  # d.bc and d.bc0 is bias-corrected eff for each obs
  d.bc=d-tbar
  d.bc0=d0-tbar0
  # z1.bc/z1.bc0 is bias-corrected z1 for each obs
  z1.bc=d.bc*z2
  z1.bc0=d.bc0*z2
  # mean eff and bias-corrected mean eff
  mu1=mean(z1)
  mu1.bc=mean(z1.bc)
  mu2=mean(z2)
  mu=mu1/mu2
  muk=mean(z1[1:nk])/mean(z2[1:nk])
  mu.bc=mu-mu.bias
  #
  mu10=mean(z10)
  mu1.bc0=mean(z1.bc0)
  mu0=mu10/mu2
  muk0=mean(z10[1:nk])/mean(z2[1:nk])
  mu.bc0=mu0-mu.bias0
  # variance and covariance
  var1=var(z1)
  var2=var(z2)
  cov12=cov(z1,z2)
  var1.bc=var(z1.bc)
  cov12.bc=cov(z1.bc,z2)
  #
  var10=var(z10)
  cov120=cov(z10,z2)
  var1.bc0=var(z1.bc0)
  cov12.bc0=cov(z1.bc0,z2)
  # variance without data sharpening
  # SZ2018 original method
  var=mu^2*(var1/(mu1^2)+var2/(mu2^2)-2*cov12/(mu1*mu2)) 
  # SZ2020 method
  var.sz=mu^2*((var1+mu1.bias^2)/(mu1^2)+var2/(mu2^2)-2*cov12/(mu1*mu2))  
  # our proposed method
  var.szz=mu.bc^2*(var1.bc/(mu1.bc^2)+var2/(mu2^2)-2*cov12.bc/(mu1.bc*mu2)) 
  # variance with data sharpening
  # SZ2018 original method
  var0=mu0^2*(var10/(mu10^2)+var2/(mu2^2)-2*cov120/(mu10*mu2)) 
  # SZ2020 method
  var.sz0=mu0^2*((var10+mu1.bias0^2)/(mu10^2)+var2/(mu2^2)-2*cov120/(mu10*mu2))  
  # our proposed method
  var.szz0=mu.bc0^2*(var1.bc0/(mu1.bc0^2)+var2/(mu2^2)-2*cov12.bc0/(mu1.bc0*mu2)) 
  #
  # SZ2018 original method
  sig1=sqrt(var)
  sig10=sqrt(var0)
  # SZ2020 correction method
  sig2=sqrt(var.sz)
  sig20=sqrt(var.sz0)
  # SZZ2023
  sig3=sqrt(var.szz)
  sig30=sqrt(var.szz0)
  #######################################################
  # for non-data sharpening
  ts10=sig10/sqrt(nk)
  ts20=sig20/sqrt(nk)
  ts30=sig30/sqrt(nk)
  # for data sharpening
  ts1=sig1/sqrt(nk)
  ts2=sig2/sqrt(nk)
  ts3=sig3/sqrt(nk)
  #
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  bounds1=matrix((muk-mu.bias-ts1*crit),nrow=3,ncol=2)
  bounds2=matrix((muk-mu.bias-ts2*crit),nrow=3,ncol=2)
  bounds3=matrix((muk-mu.bias-ts3*crit),nrow=3,ncol=2)
  #
  bounds10=matrix((muk0-mu.bias0-ts10*crit),nrow=3,ncol=2)
  bounds20=matrix((muk0-mu.bias0-ts20*crit),nrow=3,ncol=2)
  bounds30=matrix((muk0-mu.bias0-ts30*crit),nrow=3,ncol=2)
  # make a list of results to return to calling routine and then quit:
  res=list(bounds1=bounds1,bounds2=bounds2,bounds3=bounds3,
           bounds10=bounds10,bounds20=bounds20,bounds30=bounds30,
           sig=c(sig1,sig2,sig3),
           sig0=c(sig10,sig20,sig30),
           estimate=c(muk,mu.bias,muk-mu.bias),
           estimate0=c(muk0,mu.bias0,muk0-mu.bias0)
  )
  return(res)
}
