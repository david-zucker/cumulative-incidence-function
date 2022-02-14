#CIF ESTIMATES - SIMULATION STUDY

#Program for simulation study reported in 
#"Revisiting the Cumulative Incidence Function With Competing Risks Data
#by David Zucker and Malka Gorfine

#14 FEB 2022

#case of normally distributed covariate

#confidence intervals by weighted bootstrap

#bias calculation is cut off at the lowest grid point 
#that is greater than the pcut quantile of the last event time

#supremum statistic for confidence band calculation
#is cut off at the last event time

#PREAMBLE

#setwd("C:/Users/owner/Dropbox/WORK/Zucker/Malka CIF")
setwd("C:/Users/owner/Desktop")
library(survival)
library(audio)
options(max.print = 999999)
options(width=250)
set.seed(7903718)

#BASIC PARAMETER SETTINGS

#ADJUSTABLE PARAMETERS
scenario = 1.5
type = 1
ttyp = 1
n=75
nrep = 1000
nboot = 1000
rr = 6
beta1 = log(rr)
beta2 = log(rr)
znew = 1.68
censoring = F
outflg=T

#ADDITIONAL PARAMETERS
tgmax = 10
tgrid = seq(0,tgmax,by=0.05)
pcut = 0.90
nquad = 10
didl1 = 1e-4
covpr = 0.95

## SURVIVAL CURVE CONFIGURATIONS IN SIMULATIONS ##

#Cause A is the cause with the higher event rate
#Cause B is the cause with the lower event rate

#1 = Cause A hazard increasing,  Cause B hazard increasing
#2 = Cause A hazard increasing,  Cause B hazard decreasing
#3 = Cause A hazard increasing,  Cause B hazard up-and-down
#4 = Cause A hazard decreasing,  Cause B hazard increasing
#5 = Cause A hazard decreasing,  Cause B hazard decreasing
#6 = Cause A hazard decreasing,  Cause B hazard up-and-down
#7 = Cause A hazard up-and-down, Cause B hazard increasing
#8 = Cause A hazard up-and-down, Cause B hazard decreasing
#9 = Cause A hazard up-and-down, Cause B hazard up-and-down

## SUMMARY STATISTICS FOR PAPER ##
# max bias up to quantile of last event time
# sd at  quantile of last event time
# final total cif
# coverage rate for confidence band up to last event time
# mean half-width of confidence band

#TRANSFORMATION FUNCTION
trnsfn = function(u,typ.arg) {

#TYPES
#1 = IDENTITY
#2 = MODIFIED COMPLEMENTARY LOGARITHM  
#3 = LOGIT
#4 = LOG
#5 = ARCSIN ROOT

if (typ.arg==1) {
  ans=u
}   
    
if (typ.arg==2) {
  uu = pmin(u,1)
  ans = -log(1-uu+didl1)  
}  

if (typ.arg==3) {
  uu = pmin(u,1)
  ans = log(1-uu+didl1) - log(uu+didl1)  
}  
  
if (typ.arg==4) {
  ans = -log(u+didl1)  
}
  
if (typ.arg==5) {
  uu = pmax(u,0)
  uu = pmin(uu,1)
  ans = asin(sqrt(uu))
}  
  
return(ans)

}  

#CIF CURVE ESTIMATES
cif = function(t.o,delta.o,z.o,wts.o,beta1,beta2,znew) {
  
  #observation times are assumed to be ordered

  #SETUPS
  n = length(t.o)
  rr.o1 = exp(beta1*z.o)
  rr.o2 = exp(beta2*z.o)
  wrr.o1 = wts.o*rr.o1
  wrr.o2 = wts.o*rr.o2
  bresden1 = cumsum(wrr.o1[n:1])
  bresden1 = bresden1[n:1]
  bresden2 = cumsum(wrr.o2[n:1])
  bresden2 = bresden2[n:1]
  type1 = (delta.o==1)
  type2 = (delta.o==2)
  thz1 = exp(beta1*znew)
  thz2 = exp(beta2*znew)

  #METHOD 3 - NEW ESTIMATOR
  g1 = 1 - (1-wrr.o1*type1/bresden1)^(thz1/wrr.o1)
  g2 = 1 - (1-wrr.o2*type2/bresden2)^(thz2/wrr.o2)
  gtot = g1 + g2
  fac = cumprod(1-gtot)
  fac = c(1,fac[1:(n-1)])
  f1 = cumsum(wts.o*fac*g1)
  f2 = cumsum(wts.o*fac*g2)
  #if (max(f1)>1) {browser()}
  t.o = c(0,t.o)
  f1 = c(0,f1)
  f2 = c(0,f2)
  
  #METHOD 1
  bres.jmp = wts.o*(thz1*type1/bresden1 + thz2*type2/bresden2)
  fac.b = exp(-cumsum(bres.jmp))
  fac.b = c(1,fac.b[1:(n-1)])
  f1b = cumsum(fac.b*thz1*wts.o*type1/bresden1)
  f2b = cumsum(fac.b*thz2*wts.o*type2/bresden2)
  f1b = c(0,f1b)
  f2b = c(0,f2b)
  
  #METHOD 2
  bj2 = pmin(bres.jmp,1)
  fac.c = cumprod(1-bj2)
  fac.c = c(1,fac.c[1:(n-1)])
  f1c = cumsum(fac.c*thz1*wts.o*type1/bresden1)
  f2c = cumsum(fac.c*thz2*wts.o*type2/bresden2)
  f1c = c(0,f1c)
  f2c = c(0,f2c)
  
  #RETURN RESULTS
  ans = list(t.o=t.o, f1=f1, f2=f2, f1b=f1b, f2b=f2b, f1c=f1c, f2c=f2c)
  return(ans)
  
}

#ESTIMATE beta1, beta2, AND CIF FOR GIVEN WEIGHTS
estcmp = function(t.o,delta.o,z.o,znew,wts.o) {
  
  #CREATE SURVIVAL OBJECTS
  delta1.o = as.numeric((delta.o==1))
  srvobj1 = Surv(t.o,delta1.o)
  delta2.o = as.numeric((delta.o==2))
  srvobj2 = Surv(t.o,delta2.o)

  #ESTIMATE BETAS
  beta1.hat = 0
  beta2.hat = 0
  crash1 <<- F
  crash2 <<- F
  warn1 <<- F
  warn2 <<- F
  cox1 = tryCatch(coxph(srvobj1 ~ z.o, weights=wts.o),
    error = function(qq) {
      crash1 <<- T
      print(qq)
      return(qq)
    },  
    warning = function(qq) {
      warn1 <<- T
      print(qq)
      return(qq)
    })
  if (!crash1) {beta1.hat = coef(cox1)}
  cox2 = tryCatch(coxph(srvobj2 ~ z.o, weights=wts.o),
    error = function(qq) {
      crash2 <<- T
      print(qq)
      return(qq)
    }, 
    warning = function(qq) {
      warn2 <<- T
      print(qq)
      return(qq)
    })
  if (!crash2) {beta2.hat = coef(cox2)}
  crash = (crash1 | crash2)
  warn = (warn1 | warn2)
  
  #COMPUTE CIF CURVES OVER OBSERVED FOLLOW-UP ITMES
  cif.raw = cif(t.o,delta.o,z.o,wts.o,beta1.hat,beta2.hat,znew)

  #SET UP VECTORS FOR CIF CURVES OVER TIME GRID
  f1fin = rep(NA,ngrid)
  f2fin = rep(NA,ngrid)
  f1bfin = rep(NA,ngrid)
  f2bfin = rep(NA,ngrid)
  f1cfin = rep(NA,ngrid)
  f2cfin = rep(NA,ngrid)
  
  #COMPUTE CIF CURVES OVER TIME GRID
  for (ig in 1:ngrid) {
    ixt = which(cif.raw$t.o <= tgrid[ig])
    ixt1 = max(ixt)
    f1fin[ig] = cif.raw$f1[ixt1]
    f2fin[ig] = cif.raw$f2[ixt1]
    f1bfin[ig] = cif.raw$f1b[ixt1]
    f2bfin[ig] = cif.raw$f2b[ixt1]
    f1cfin[ig] = cif.raw$f1c[ixt1]
    f2cfin[ig] = cif.raw$f2c[ixt1]
  }  
  
  #END-OF-STUDY TOTAL CIF
  f.sum.end = cif.raw$f1[n+1] + cif.raw$f2[n+1]
  f.sum.end.b = cif.raw$f1b[n+1] + cif.raw$f2b[n+1]
  f.sum.end.c = cif.raw$f1c[n+1] + cif.raw$f2c[n+1]
  
  #RETURN RESULTS
  ans = list(t.o=cif.raw$t.o, f1=cif.raw$f1, f2=cif.raw$f2, f1b=cif.raw$f1b, 
    f2b=cif.raw$f2b, f1c=cif.raw$f1c, f2c=cif.raw$f2c, f1fin=f1fin, f2fin=f2fin,
    f1bfin=f1bfin, f2bfin=f2bfin, f1cfin=f1cfin, f2cfin=f2cfin, f.sum.end=f.sum.end, 
    f.sum.end.b=f.sum.end.b, f.sum.end.c=f.sum.end.c, crash=crash, warn=warn)
  return(ans)

}

#SIMPSON'S RULE INTEGRATION
#THE R FUNCTION integrate WAS GIVING ME PROBLEMS
simp.int = function(f,pars1,pars2,flg,xl,xu,npts) {
  h = (xu-xl)/npts
  x = xl + (0:npts)*h
  wts.simp = rep(2,npts-1) + 2*((1:(npts-1)) %% 2)
  wts.simp = c(1,wts.simp,1)
  y = rep(0,npts+1)
  for (i in 1:(npts+1)) {
    y[i] = f(x[i],pars1,pars2,flg)
  }
  ans = (h/3)*sum(wts.simp*y)
  return(ans)
}

#HAZARD FUNCTION
hazfn = function(t,pars) {
  p = pars[1]
  a = pars[2]
  b = pars[3]
  sig = pars[4]
  tt = t+a
  hazval = sig*p*(tt^(p-1))/(1+b*(tt^p))
  return(hazval)
}  

#CUMULATIVE HAZARD FUNCTION
cmlhaz = function(tval,pars) {
  if (tval==0) {return(0)}
  p = pars[1]
  a = pars[2]
  b = pars[3]
  sig = pars[4]
  weiflg = pars[5]
  tval1 = tval + a
  if (weiflg) {
    chval = sig*((tval1^p) - (a^p))
  }
  if (!weiflg) {
    chval = (sig/b)* (log(1+b*(tval1^p)) - log(1+b*(a^p)))
  }    
  return(chval)
}  

#INTEGRAND FOR TRUE CIF COMPUTATION 
cif.tru.int = function(tt,pars1,pars2,flg) {
  stru = exp(-rr1.tru*cmlhaz(tt,pars1)) * exp(-rr2.tru*cmlhaz(tt,pars2))
  if (flg==1) {haz = rr1.tru * hazfn(tt,pars1)}
  if (flg==2) {haz = rr2.tru * hazfn(tt,pars2)}
  intgrnd = stru*haz   
  return(intgrnd)
}  

###############################################################################
# MAIN PROGRAM
###############################################################################

tim0 = proc.time()

#SURVIVAL CURVE PARAMETERS
#pars = (p,a,b,sig,weiflg)

if (type == 1) {
  typ.lbl = '1 = Cause A hazard increasing,  Cause B hazard increasing'
  cens.rate = 0.35 #0.40
  sigsum =  0.08   #0.12
  sigfrac = 0.65
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(3,0,0,sig1,TRUE)
  pars2 = c(3,0,0,sig2,TRUE)
}  

if (type == 2) {
  typ.lbl = '2 = Cause A hazard increasing,  Cause B hazard decreasing'
  cens.rate = 0.36
  sigsum = 0.66
  sigfrac = 0.08
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(3,0,0,sig1,TRUE)
  pars2 = c(0.3,0.4,0,sig2,TRUE)
}  

if (type == 3) {
  typ.lbl = '3 = Cause A hazard increasing,  Cause B hazard up-and-down'
  cens.rate = 0.55
  sigsum = 0.45
  sigfrac = 0.44
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(3,0,0,sig1,TRUE)
  pars2 = c(3,0,0.75,sig2,FALSE)
}  

if (type == 4) {
  typ.lbl = '4 = Cause A hazard decreasing,  Cause B increasing'
  cens.rate = 1.13
  sigsum = 3.07
  sigfrac = 0.84
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(0.3,0.4,0,sig1,TRUE)
  pars2 = c(3,0,0,sig2,TRUE)
}  

if (type == 5) {
  typ.lbl = '5 = Cause A hazard decreasing, Cause B decreasing'
  cens.rate = 1.9
  sigsum = 3.00 
  sigfrac = 0.65
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(0.5,0.4,0,sig1,TRUE)
  pars2 = c(0.5,0.4,0,sig2,TRUE)
}  

if (type == 6) {
  typ.lbl = '6 = Cause A hazard decreasing, Cause B up-and-down'
  cens.rate = 0.74
  sigsum = 2.13
  sigfrac = 0.86
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(0.3,0.4,0,sig1,TRUE)
  pars2 = c(3,0,0.75,sig2,FALSE)
}  

if (type == 7) {
  typ.lbl = '7 = Cause A hazard up-and-down, Cause B increasing'
  cens.rate = 0.66
  sigsum = 0.76
  sigfrac = 0.80
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(3,0,0.75,sig1,FALSE)
  pars2 = c(3,0,0,sig2,TRUE)
}  

if (type == 8) {
  typ.lbl = '8 = Cause A hazard up-and-down, Cause B decreasing'
  cens.rate = 0.44
  sigsum = 0.79
  sigfrac = 0.37
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(3,0,0.75,sig1,FALSE)
  pars2 = c(0.3,0.4,0,sig2,TRUE)
}  

if (type == 9) {
  typ.lbl = '9 = Cause A hazard up-and-down, Cause B up-and-down'
  cens.rate = 0.83
  sigsum = 1.50 #0.74
  sigfrac = 0.65
  sig1 = sigfrac*sigsum
  sig2 = (1-sigfrac)*sigsum
  pars1 = c(3,0,0.75,sig1,FALSE)
  pars2 = c(3,0,0.75,sig2,FALSE)
}  

### SETUPS FOR SIMULATION RESULTS ###

ngrid = length(tgrid)
f1tru = rep(0,ngrid)
f2tru = rep(0,ngrid)
at.risk.mat = matrix(0,nrep,ngrid)
ev1.mat = matrix(0,nrep,ngrid)
ev2.mat = matrix(0,nrep,ngrid)
cens.mat = matrix(0,nrep,ngrid)
t.end.vec = rep(0,n)

f1mat = NULL
f2mat = NULL
f1bmat = NULL
f2bmat = NULL
f1cmat = NULL
f2cmat = NULL

cvr1 = NULL
cvr2 = NULL
cvr1b = NULL
cvr2b = NULL
cvr1c = NULL
cvr2c = NULL

d1crit.vec = NULL
d2crit.vec = NULL
d1bcrit.vec = NULL
d2bcrit.vec = NULL
d1ccrit.vec = NULL
d2ccrit.vec = NULL

f.sum.end = rep(0,nrep)
f.sum.end.b = rep(0,nrep)
f.sum.end.c = rep(0,nrep)

errtik1 = 0
errtik2 = 0
wrntik1 = 0
wrntik2 = 0
cxcr.vec = NULL
cxwn.vec = NULL

## END OF SETUPS FOR SIMULATION RESULTS ##

#TRUE CIF COMPUTATION OVER GRID
rr1.tru = exp(beta1*znew)
rr2.tru = exp(beta2*znew)
t.old = 0
for (r in 2:ngrid) {
  tcur = tgrid[r]
  f1tru[r] = f1tru[r-1] + simp.int(cif.tru.int,pars1,pars2,1,t.old,tcur,nquad)
  f2tru[r] = f2tru[r-1] + simp.int(cif.tru.int,pars1,pars2,2,t.old,tcur,nquad)
  t.old = tcur
}
ftot.at.end = f1tru[ngrid] + f2tru[ngrid]
f1tru = f1tru / ftot.at.end
f2tru = f2tru / ftot.at.end
f1tru.tr = trnsfn(f1tru,ttyp)
f2tru.tr = trnsfn(f2tru,ttyp)

n1 = round(1.25*n)

#SIMULATION LOOP
for (irep in 1:nrep) {

if (irep %% 10 == 0) {print(irep)}

#GENERATE DATA

#z = runif(n1) - 0.5
z = rnorm(n1,0,2)
rr1 = exp(beta1*z)
rr2 = exp(beta2*z)

u1 = runif(n1)
trg1 = -log(u1)/(rr1*pars1[4])
p = pars1[1]
a = pars1[2]
b = pars1[3]
bfac = 1 + b*(a^p)
p.inv = 1/p
if (pars1[5]==TRUE) {
  t1 = (trg1+(a^p))^p.inv - a
}
if (pars1[5]==FALSE) {
  t1 = ((bfac*exp(b*trg1)-1)/b)^p.inv - a
}

u2 = runif(n1)
trg2 = -log(u2)/(rr2*pars2[4])
p = pars2[1]
a = pars2[2]
b = pars2[3]
p.inv = 1/p
bfac = 1 + b*(a^p)
if (pars2[5]==TRUE) {
  t2 = (trg2+(a^p))^p.inv - a
}
if (pars2[5]==FALSE) {
  t2 =  ((bfac*exp(b*trg2)-1)/b)^p.inv - a
}  

te.ini = pmin(t1,t2)
delta = ifelse(t1<t2,1,2)
ixgud = which(te.ini <= tgmax)
ixgud = ixgud[1:n]
z = z[ixgud]
te = te.ini[ixgud]
delta = delta[ixgud]

tc = Inf
if (censoring) {
  u3 = runif(n)
  tc = -log(u3)/cens.rate
}
t = pmin(te,tc)
delta = ifelse(te<tc,delta,0)

#EVENT COUNTS
for (r in 1:ngrid) {
  at.risk.mat[irep,r] = sum(t >= tgrid[r])
  ev1.mat[irep,r] = sum((t <= tgrid[r]) & (delta==1)) 
  ev2.mat[irep,r] = sum((t <= tgrid[r]) & (delta==2))
  cens.mat[irep,r] = sum((t <= tgrid[r]) & (delta==0))
}

na.flg = sum(is.na(at.risk.mat[irep,]))
if (na.flg > 1) {
  print("NA encountered")
  print(irep)
  browser()
}

#END OF DATA GENERATION

#OBSERVATIONS SORTED BY FOLLOW-UP TIME
ord = order(t)
t.o = t[ord]
delta.o = delta[ord]
z.o = z[ord]

#IDENTIFY THE END
ix.end = max(which(delta.o > 0))
t.end = t.o[ix.end]
t.end.vec[irep] = t.end
ix.grd.end = min(which(tgrid >= t.end))
ix.range = 1:ix.grd.end

#COMPUTE CIF ESTIMATES FOR ORIGINAL DATA
wts.o = rep(1,n)
cif.orig = estcmp(t.o,delta.o,z.o,znew,wts.o)
coxcrash1 = (cif.orig$crash)
errtik1 = errtik1 + coxcrash1
wrntik1 = wrntik1 + cif.orig$warn

if (!coxcrash1) {

  #APPEND RESULTS FOR THIS REPLICATION TO THE RESULTS MATRICES
  f1mat = rbind(f1mat,cif.orig$f1fin)
  f2mat = rbind(f2mat,cif.orig$f2fin)
  f1bmat = rbind(f1bmat,cif.orig$f1bfin)
  f2bmat = rbind(f2bmat,cif.orig$f2bfin)
  f1cmat = rbind(f1cmat,cif.orig$f1cfin)
  f2cmat = rbind(f2cmat,cif.orig$f2cfin)
  
  #END-OF-STUDY TOTAL CIF
  f.sum.end[irep] = cif.orig$f.sum.end
  f.sum.end.b[irep] = cif.orig$f.sum.end.b
  f.sum.end.c[irep] = cif.orig$f.sum.end.c

  #TAKE TRANSFORMATION
  f1.tr = trnsfn(cif.orig$f1fin,ttyp)
  f2.tr = trnsfn(cif.orig$f2fin,ttyp)
  f1b.tr = trnsfn(cif.orig$f1bfin,ttyp)
  f2b.tr = trnsfn(cif.orig$f2bfin,ttyp)
  f1c.tr = trnsfn(cif.orig$f1cfin,ttyp)
  f2c.tr = trnsfn(cif.orig$f2cfin,ttyp)

  #BOOTSTRAP LOOP
  f1mat.boot = NULL
  f2mat.boot = NULL
  f1bmat.boot = NULL
  f2bmat.boot = NULL
  f1cmat.boot = NULL
  f2cmat.boot = NULL

  cxcr = 0
  cxwn = 0
  for (ib in 1:nboot) {
    wts.o = rexp(n,1)
    wts.o = wts.o / mean(wts.o)
    cif.boot = estcmp(t.o,delta.o,z.o,znew,wts.o)
    coxcrash2 = (cif.boot$crash)
    cxcr = cxcr + coxcrash2
    cxwn = cxwn + (cif.boot$warn)
    if (!coxcrash2) {
      f1mat.boot =  rbind(f1mat.boot,cif.boot$f1fin)
      f2mat.boot =  rbind(f2mat.boot,cif.boot$f2fin)
      f1bmat.boot = rbind(f1bmat.boot,cif.boot$f1bfin)
      f2bmat.boot = rbind(f2bmat.boot,cif.boot$f2bfin)
      f1cmat.boot = rbind(f1cmat.boot,cif.boot$f1cfin)
      f2cmat.boot = rbind(f2cmat.boot,cif.boot$f2cfin)
    }  
  }

  errtik2 = errtik2 + (cxcr > 0)
  cxcr.vec = c(cxcr.vec,cxcr)
  wrntik2 = wrntik2 + (cxwn > 0)
  cxwn.vec = c(cxwn.vec,cxwn)
  nboot.adj = nboot - cxcr

  #TAKE TRANSFORMATION
  f1mat.boot.tr = trnsfn(f1mat.boot,ttyp)
  f2mat.boot.tr = trnsfn(f2mat.boot,ttyp)
  f1bmat.boot.tr = trnsfn(f1bmat.boot,ttyp)
  f2bmat.boot.tr = trnsfn(f2bmat.boot,ttyp)
  f1cmat.boot.tr = trnsfn(f1cmat.boot,ttyp)
  f2cmat.boot.tr = trnsfn(f2cmat.boot,ttyp)

  #SUPREMUM STATISTICS IN BOOTSTRAP SAMPLES
  sup1 = NULL
  sup2 = NULL
  sup1b = NULL
  sup2b = NULL
  sup1c = NULL
  sup2c = NULL
  for (ib in 1:nboot.adj) { 
    sup1.new =  max(abs(f1mat.boot.tr[ib,ix.range]-f1.tr[ix.range]))
    sup2.new =  max(abs(f2mat.boot.tr[ib,ix.range]-f2.tr[ix.range]))
    sup1b.new = max(abs(f1bmat.boot.tr[ib,ix.range]-f1b.tr[ix.range]))
    sup2b.new = max(abs(f2bmat.boot.tr[ib,ix.range]-f2b.tr[ix.range]))
    sup1c.new = max(abs(f1cmat.boot.tr[ib,ix.range]-f1c.tr[ix.range]))
    sup2c.new = max(abs(f2cmat.boot.tr[ib,ix.range]-f2c.tr[ix.range]))
    sup1 = rbind(sup1,sup1.new)
    sup2 = rbind(sup2,sup2.new)
    sup1b = rbind(sup1b,sup1b.new)
    sup2b = rbind(sup2b,sup2b.new)
    sup1c = rbind(sup1c,sup1c.new)
    sup2c = rbind(sup2c,sup2c.new)
  }  

  #CRITICAL VALUES FOR BOOTSTRAP CONFIDENCE BANDS
  d1crit = quantile(sup1,covpr,na.rm=T)
  d2crit = quantile(sup2,covpr,na.rm=T)
  d1bcrit = quantile(sup1b,covpr,na.rm=T)
  d2bcrit = quantile(sup2b,covpr,na.rm=T)
  d1ccrit = quantile(sup1c,covpr,na.rm=T)
  d2ccrit = quantile(sup2c,covpr,na.rm=T)
  d1crit.vec = c(d1crit.vec,d1crit)
  d2crit.vec = c(d2crit.vec,d2crit)
  d1bcrit.vec = c(d1bcrit.vec,d1bcrit)
  d2bcrit.vec = c(d2bcrit.vec,d2bcrit)
  d1ccrit.vec = c(d1ccrit.vec,d1ccrit)
  d2ccrit.vec = c(d2ccrit.vec,d2ccrit)

  #COVERAGE STATUS OF BOOTSTRAP CONFIDENCE BANDS
  sup1.orig = max(abs(f1.tr[ix.range]-f1tru.tr[ix.range]))
  sup2.orig = max(abs(f2.tr[ix.range]-f2tru.tr[ix.range]))
  sup1b.orig = max(abs(f1b.tr[ix.range]-f1tru.tr[ix.range]))
  sup2b.orig = max(abs(f2b.tr[ix.range]-f2tru.tr[ix.range]))
  sup1c.orig = max(abs(f1c.tr[ix.range]-f1tru.tr[ix.range]))
  sup2c.orig = max(abs(f2c.tr[ix.range]-f2tru.tr[ix.range]))
  cvr1 = c(cvr1, (sup1.orig <= d1crit))
  cvr2 = c(cvr2, (sup2.orig <= d2crit))
  cvr1b = c(cvr1b, (sup1b.orig <= d1bcrit))
  cvr2b = c(cvr2b, (sup2b.orig <= d2bcrit))
  cvr1c = c(cvr1c, (sup1c.orig <= d1ccrit))
  cvr2c = c(cvr2c, (sup2c.orig <= d2ccrit))

}
  
}
# END OF SIMULATION LOOP #

#COMPUTE MEAN EVENT COUNTS OVER GRID
at.risk = colMeans(at.risk.mat)
ev1 = colMeans(ev1.mat)
ev2 = colMeans(ev2.mat)
cens = colMeans(cens.mat)
at.risk = round(at.risk,digits=1)
ev1 = round(ev1,digits=1)
ev2 = round(ev2,digits=1)
cens = round(cens,digits=1)

#COMPUTE MEAN OF ESTIMATES OVER GRID
f1mean = colMeans(f1mat,na.rm=T)
f2mean = colMeans(f2mat,na.rm=T)
f1bmean = colMeans(f1bmat,na.rm=T)
f2bmean = colMeans(f2bmat,na.rm=T)
f1cmean = colMeans(f1cmat,na.rm=T)
f2cmean = colMeans(f2cmat,na.rm=T)

#COMPUTE SD OF ESTIMATES OVER GRID
f1sd = apply(f1mat,2,sd,na.rm=T)
f2sd = apply(f2mat,2,sd,na.rm=T)
f1bsd = apply(f1bmat,2,sd,na.rm=T)
f2bsd = apply(f2bmat,2,sd,na.rm=T)
f1csd = apply(f1cmat,2,sd,na.rm=T)
f2csd = apply(f2cmat,2,sd,na.rm=T)

#COMPUTE MEAN ERROR OF ESTIMATES OVER GRID
f1dif = f1mean-f1tru
f2dif = f2mean-f2tru
f1bdif = f1bmean-f1tru
f2bdif = f2bmean-f2tru
f1cdif = f1cmean-f1tru
f2cdif = f2cmean-f2tru

#SUM OF CAUSE A CDF AND CAUSE B CDF
f.sum.mean = f1mean + f2mean
fb.sum.mean = f1bmean + f2bmean
fc.sum.mean = f1cmean + f2cmean

### OUTPUT RESULTS OVER GRID ###

ans = cbind(tgrid,at.risk,ev1,ev2,cens,f1tru,f1bdif,f1cdif,f1dif,
  f2tru,f2bdif,f2cdif,f2dif,fb.sum.mean,fc.sum.mean,f.sum.mean)
ans1 = round(ans,digits=4)
print(ans1)
cat('\n')
print('errtik1=')
print(errtik1)
print('errtik2=')
print(errtik2)
print('wrntik1=')
print(wrntik1)
print('wrntik2=')
print(wrntik2)
cat('\n')
print('summary of bootstrap failures')
print(summary(cxcr.vec))
cat('\n')
print('summary of bootstrap warnings')
print(summary(cxwn.vec))
cat('\n')

## SUMMARY STATISTICS FOR PAPER ##

t.cut = quantile(t.end.vec,probs=pcut)
ix.grd.end.fin = min(which(tgrid >= t.cut))
ix.range.fin = 1:ix.grd.end.fin

#MAX BIAS UP TO CUTOFF TIME
f1.max.bias = max(abs(f1dif[ix.range.fin]))
f1b.max.bias = max(abs(f1bdif[ix.range.fin]))
f1c.max.bias = max(abs(f1cdif[ix.range.fin]))
f2.max.bias = max(abs(f2dif[ix.range.fin]))
f2b.max.bias = max(abs(f2bdif[ix.range.fin]))
f2c.max.bias = max(abs(f2cdif[ix.range.fin]))
f1.max.bias1 = round(f1.max.bias,digits=4)
f1b.max.bias1 = round(f1b.max.bias,digits=4)
f1c.max.bias1 = round(f1c.max.bias,digits=4)
f2.max.bias1 = round(f2.max.bias,digits=4)
f2b.max.bias1 = round(f2b.max.bias,digits=4)
f2c.max.bias1 = round(f2c.max.bias,digits=4)

#SD AT CUTOFF TIME
f1sd.end = round(f1sd[ix.grd.end.fin],digits=4)
f2sd.end = round(f2sd[ix.grd.end.fin],digits=4)
f1bsd.end = round(f1bsd[ix.grd.end.fin],digits=4)
f2bsd.end = round(f2bsd[ix.grd.end.fin],digits=4)
f1csd.end = round(f1csd[ix.grd.end.fin],digits=4)
f2csd.end = round(f2csd[ix.grd.end.fin],digits=4)

#CONFIDENCE BAND COVERAGE RATES
nrep.adj = nrep - errtik1
cvrte1 = mean(cvr1,na.rm=T)
cvrte2 = mean(cvr2,na.rm=T)
cvrte1b = mean(cvr1b,na.rm=T)
cvrte2b = mean(cvr2b,na.rm=T)
cvrte1c = mean(cvr1c,na.rm=T)
cvrte2c = mean(cvr2c,na.rm=T)
cvrte1.r = round(cvrte1,digits=3)
cvrte2.r = round(cvrte2,digits=3)
cvrte1b.r = round(cvrte1b,digits=3)
cvrte2b.r = round(cvrte2b,digits=3)
cvrte1c.r = round(cvrte1c,digits=3)
cvrte2c.r = round(cvrte2c,digits=3)

#MEAN CI HALF-WIDTH
d1crit.mean = round(mean(d1crit.vec,na.rm=T),digits=4)
d2crit.mean = round(mean(d2crit.vec,na.rm=T),digits=4)
d1bcrit.mean = round(mean(d1bcrit.vec,na.rm=T),digits=4)
d2bcrit.mean = round(mean(d2bcrit.vec,na.rm=T),digits=4)
d1ccrit.mean = round(mean(d1ccrit.vec,na.rm=T),digits=4)
d2ccrit.mean = round(mean(d2ccrit.vec,na.rm=T),digits=4)

#QUANTILES OF FINAL TOTAL CIF
f.sum.end.b.quant = round(quantile(f.sum.end.b,probs=c(0.01,0.10,0.50,0.90,0.99)),digits=4)
f.sum.end.c.quant = round(quantile(f.sum.end.c,probs=c(0.01,0.10,0.50,0.90,0.99)),digits=4)

our.output1 = cbind(scenario,type,n,rr,znew,censoring,nrep,f1b.max.bias1,f1c.max.bias1,f1.max.bias1,
  f2b.max.bias1,f2c.max.bias1,f2.max.bias1)                  
our.output2 = cbind(scenario,type,n,rr,znew,censoring,nrep,f1bsd.end,f1csd.end,f1sd.end,
  f2bsd.end,f2csd.end,f2sd.end)
our.output3 = cbind(scenario,type,n,rr,znew,censoring,nrep,nboot,cvrte1b.r,cvrte1c.r,cvrte1.r,
  cvrte2b.r,cvrte2c.r,cvrte2.r)
our.output4 = cbind(scenario,type,n,rr,znew,censoring,nrep,nboot,d1bcrit.mean,d1ccrit.mean,d1crit.mean,
  d2bcrit.mean,d2ccrit.mean,d2crit.mean)                                    
our.output5 = cbind(scenario,type,n,rr,znew,nrep,t(f.sum.end.b.quant),t(f.sum.end.c.quant))

if (outflg) {
  write(our.output1, file = "aaa.txt", ncolumns=13, append=T, sep="        ")
  write(our.output2, file = "bbb.txt", ncolumns=13, append=T, sep="        ")
  write(our.output3, file = "ccc.txt", ncolumns=14, append=T, sep="     ")
  write(our.output4, file = "ddd.txt", ncolumns=16, append=T, sep="   ")
  if (!censoring) {write(our.output5, file = "eee.txt", ncolumns=14, append=T, sep="     ")}
}

print(cbind(n,beta1,beta2,znew))
cat('\n')
print(our.output1)
cat('\n')
print(our.output2)
cat('\n')
print(our.output3)
cat('\n')
print(our.output4)
cat('\n')
print(our.output5)
cat('\n')

tim1 = proc.time()
print(tim1-tim0)
play(sin(1:10000/20))
