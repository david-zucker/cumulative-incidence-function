#CUMULATIVE INCIDENCE CALCULATION

#Program for analyzing a dataset using the methods reported in 
#"Revisiting the Cumulative Incidence Function With Competing Risks Data
#by David Zucker and Malka Gorfine

#14 FEB 22


#FUNCTION CALL:
# cif(t,J,delta,z,znew,nboot,covpr,seed,tgrid)

#t = observation times
#J = number of event types
#delta = event type (0, 1, ..., J with 0 for censored)
#z = n x k matrix of covariate values (k = number of covariates)
#znew = vector of length k of covariate values for which we want to estimate the CIF's
#nboot = number of bootstrap replications for confidence band
#cvpr = desired coverage rate of confidence band
#seed = random number seed for bootstrapping
#tgrid = grid of time points for which CIF estimates are to be produced

###############################################################################

library(survival)
options(max.print = 999999)
options(width=200)

#CIF CURVE ESTIMATE COMPUTATION
cifcmp = function(t.o,J,delta.o,z.o,znew,tgrid,wts.o) {

  #observation times are assumed to be ordered

  #SETUP FOR COX MODEL RUNS
  k = ncol(z.o)
  vname = 'z1'
  assign(vname,z.o[,1])
  frmla = 'survobj ~ z1'
  for (kk in 2:k) {
    vname = paste('z',kk,sep='')
    assign(vname,z.o[,kk])
    frmla = paste(frmla,'+')
    frmla = paste(frmla,vname)
  }  
  frmla = as.formula(frmla)

  #RUN COX MODELS AND EXTRACT BETAS
  beta = matrix(0,k,J)
  for (j in 1:J) {
    delcur = (delta.o == j)
    survobj = Surv(t.o,delcur)
    coxreg = coxph(frmla, ties='breslow', weights=wts.o)
    beta[,j] = coef(coxreg)
  }

  #SETUPS
  n = length(t)
  t.o.u = unique(t.o)
  ntou = length(t.o.u)
  rr.o = exp(z.o %*% beta) #n x J matrix
  bresden = matrix(0,ntou,J)
  for (i in 1:ntou) {
    atrsk = which(t.o >= t.o.u[i])
    for (j in 1:J) {
      bresden[i,j] = sum(wts.o[atrsk]*rr.o[atrsk,j])
    }
  }
  thz = exp(t(znew) %*% beta)
  evnt = matrix(0,n,J)
  for (i in 1:n) {
  for (j in 1:J) {
    evnt[i,j] = (delta.o[i] == j)
  }
  }
  cif.a = matrix(0,ntou,J)
  cif.b = matrix(0,ntou,J)
  cif.c = matrix(0,ntou,J)

  #PROPOSED ESTIMATE
  new.wts = matrix(0,ntou,J)
  g = matrix(0,ntou,J)
  for (i in 1:ntou) {
    ix = which(t.o == t.o.u[i])
    for (j in 1:J) {
      new.wts[i,j] = sum(wts.o[ix]*evnt[ix,j])
      numer = sum(wts.o[ix]*evnt[ix,j]*rr.o[ix,j])
      g[i,j] = 1 - (1-numer/bresden[i,j])^(thz[j]/numer)
    }
  }
  gtot = rowSums(g)
  fac.a = cumprod(1-gtot)
  fac.a = c(1,fac.a[1:(ntou-1)])
  for (j in 1:J) {
    cif.a[,j] = cumsum(new.wts[,j]*fac.a*g[,j])
  }

  #ALTERNATIVE ESTIMATE 1
  bres.jmp = rep(0,ntou)
  for (j in 1:J) {
    bres.jmp = bres.jmp + thz[j]*new.wts[,j]/bresden[,j]
  }
  fac.b = exp(-cumsum(bres.jmp))
  fac.b = c(1,fac.b[1:(ntou-1)])
  for (j in 1:J) {
    cif.b[,j] = cumsum(fac.b*thz[j]*new.wts[,j]/bresden[,j])
  }

  #ALTERNATIVE ESTIMATE 2
  bj2 = pmin(bres.jmp,1)
  fac.c = cumprod(1-bj2)
  fac.c = c(1,fac.c[1:(ntou-1)])
  for (j in 1:J) {
    cif.c[,j] = cumsum(fac.c*thz[j]*new.wts[,j]/bresden[,j])
  }

  #CIF ESTIMATES OVER GRID
  ngrid = length(tgrid)
  cif.fin.a = matrix(0,ngrid,J)
  cif.fin.b = matrix(0,ngrid,J)
  cif.fin.c = matrix(0,ngrid,J)
  for (ig in 1:ngrid) {
    ixt = which(t.o.u <= tgrid[ig])
    if (length(ixt)>0) {
      ixt1 = max(ixt)
      for (j in 1:J) {
        cif.fin.a[ig,j] = cif.a[ixt1,j]
        cif.fin.b[ig,j] = cif.b[ixt1,j]
        cif.fin.c[ig,j] = cif.c[ixt1,j]
      }
    }
  }
    
  #RETURN RESULTS
  ans = list(time=t.o.u, beta=beta, cif.a=cif.a, cif.b=cif.b, cif.c=cif.c,
    cif.fin.a=cif.fin.a, cif.fin.b=cif.fin.b, cif.fin.c=cif.fin.c)
  return(ans)

}

#CIF ANALYSIS
cif = function(t,J,delta,z,znew,nboot,covpr,seed,tgrid) {

  #t = observation times
  #J = number of event types
  #delta = event type (0, 1, ..., J with 0 for censored)
  #z = n x k matrix of covariate values (k = number of covariates)
  #znew = vector of length k of covariate values for which we want to estimate the CIF's
  #nboot = number of bootstrap replications for confidence band
  #cvpr = desired coverage rate of confidence band
  #seed = random number seed for bootstrapping
  #tgrid = grid of time points for which CIF estimates are to be produced
  
  #OBSERVATIONS SORTED BY FOLLOW-UP TIME
  ord = order(t)
  t.o = t[ord]
  delta.o = delta[ord]
  z.o = z[ord,]

  #IDENTIFY THE END
  ix.end = max(which(delta.o > 0))
  t.end = t.o[ix.end]
  ix.grd = which(tgrid >= t.end)
  if (length(ix.grd) > 0) {
    ix.grd.end = min(ix.grd)
    ix.range = 1:ix.grd.end
  }
  else {
    ix.range = 1:length(tgrid)
  }   
    
  #COMPUTE CIF ESTIMATES FOR ORIGINAL DATA
  n = length(t)
  cif.orig = cifcmp(t.o,J,delta.o,z.o,znew,tgrid,rep(1,n))
  
  #BOOTSTRAP LOOP
  set.seed(seed)
  sup.a = matrix(0,nboot,J)
  sup.b = matrix(0,nboot,J)
  sup.c = matrix(0,nboot,J)
  for (ib in 1:nboot) {
    wts.o = rexp(n,1)
    wts.o = wts.o / mean(wts.o)
    cif.boot = cifcmp(t.o,J,delta.o,z.o,znew,tgrid,wts.o)
    for (j in 1:J) {
      sup.a[ib,j] = max(abs(cif.boot$cif.fin.a[ix.range,j]-cif.orig$cif.fin.a[ix.range,j]))
      sup.b[ib,j] = max(abs(cif.boot$cif.fin.b[ix.range,j]-cif.orig$cif.fin.b[ix.range,j]))
      sup.c[ib,j] = max(abs(cif.boot$cif.fin.c[ix.range,j]-cif.orig$cif.fin.c[ix.range,j]))
    }
  }

  #CRITICAL VALUES FOR BOOTSTRAP CONFIDENCE BANDS
  dcrit.a = rep(0,J)
  dcrit.b = rep(0,J)
  dcrit.c = rep(0,J)
  for (j in 1:J) {
    dcrit.a[j] = quantile(sup.a[,j],covpr)
    dcrit.b[j] = quantile(sup.b[,j],covpr)
    dcrit.c[j] = quantile(sup.c[,j],covpr)
  }

  #LOWER AND UPPER LIMITS OF CONFIDENCE BANDS
  ngrid = length(tgrid)
  cif.fin.a.lo = matrix(0,ngrid,J)
  cif.fin.b.lo = matrix(0,ngrid,J)
  cif.fin.c.lo = matrix(0,ngrid,J)
  cif.fin.a.hi = matrix(0,ngrid,J)
  cif.fin.b.hi = matrix(0,ngrid,J)
  cif.fin.c.hi = matrix(0,ngrid,J)
  for (i in 1:ngrid) {
  for (j in 1:J) {
    cif.fin.a.lo[i,j] = cif.orig$cif.fin.a[i,j] - dcrit.a[j]
    cif.fin.a.hi[i,j] = cif.orig$cif.fin.a[i,j] + dcrit.a[j]
    cif.fin.b.lo[i,j] = cif.orig$cif.fin.b[i,j] - dcrit.b[j]
    cif.fin.b.hi[i,j] = cif.orig$cif.fin.b[i,j] + dcrit.b[j]
    cif.fin.c.lo[i,j] = cif.orig$cif.fin.c[i,j] - dcrit.c[j]
    cif.fin.c.hi[i,j] = cif.orig$cif.fin.c[i,j] + dcrit.c[j]
   }
  }
  cif.fin.a.lo = pmax(cif.fin.a.lo,0)
  cif.fin.b.lo = pmax(cif.fin.b.lo,0)
  cif.fin.c.lo = pmax(cif.fin.c.lo,0)
  cif.fin.a.hi = pmin(cif.fin.a.hi,1)
  cif.fin.b.hi = pmin(cif.fin.b.hi,1)
  cif.fin.c.hi = pmin(cif.fin.c.hi,1)
  
  cifobj = cif.orig
  cifobj[['cif.fin.a.lo']] = cif.fin.a.lo
  cifobj[['cif.fin.a.hi']] = cif.fin.a.hi
  cifobj[['cif.fin.b.lo']] = cif.fin.b.lo
  cifobj[['cif.fin.b.hi']] = cif.fin.b.hi
  cifobj[['cif.fin.c.lo']] = cif.fin.c.lo
  cifobj[['cif.fin.c.hi']] = cif.fin.c.hi
   
  return(cifobj)

}
