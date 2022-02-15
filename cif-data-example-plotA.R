library(tidyverse)
library(survival)

## "SOURCE" THE CODE WITH THE CIF FUNCTIONS ##
setwd("~/Dropbox/Malka/CIF-newproject")
source('cif-data-analysis-revised.r')

dat <- read.csv("data_snippet103.csv",header = TRUE)

dat$delta <- dat$correct+1
dat$sex <- dat$gender=="male"
dat$st <- ifelse(is.na(dat$selftaught), FALSE, dat$selftaught)
dat$yoe <- ifelse(is.na(dat$yearsofexperience), 3, dat$yearsofexperience)
dat$st <- ifelse(dat$st==TRUE,1,0)

z = cbind(dat$order ,dat$age ,dat$sex  ,dat$yoe )

#SETUP FOR COX MODEL RUNS
k = ncol(z)
vname = 'z1'
assign(vname,z[,1])
frmla = 'survobj ~ z1'
for (kk in 2:k) {
  vname = paste('z',kk,sep='')
  assign(vname,z[,kk])
  frmla = paste(frmla,'+')
  frmla = paste(frmla,vname)
}  
frmla = as.formula(frmla)

#RUN COX MODELS AND EXTRACT BETAS
J=2
beta = matrix(0,k,J)
for (j in 1:J) {
  delcur = (dat$delta == j)
  survobj = Surv(dat$t,delcur)
  coxreg = coxph(frmla, ties='breslow')
  print(coxreg)
  beta[,j] = coef(coxreg)
}

znew1 = c(1 , 35 , 0  , 0 ) 
znew2 = c(1 , 35 , 1  , 0 ) 
znew3 = c(1 , 35 , 0  , 5 ) 
znew4 = c(1 , 35 , 1  , 5 ) 
znew5 = c(10 , 35 , 0  , 0 ) 
znew6 = c(10 , 35 , 1  , 0 ) 
znew7 = c(10 , 35 , 0  , 5 ) 
znew8 = c(10 , 35 , 1  , 5 ) 
znew.mat <- rbind(znew1,znew2,znew3,znew4,znew5,znew6,znew7,znew8)

#max.ev.tim = max(dat$time[which(dat$delta>0)])
tgrid = sort(dat$time)

for (i in 1:8) {
  znew <- znew.mat[i,]
  rottcif = cif(dat$time,2,dat$delta,z,znew,100,0.95,1234,tgrid,0.1)
  
  mean(rottcif$cif.fin.a-rottcif$cif.fin.b )
  mean(rottcif$cif.fin.a-rottcif$cif.fin.c )
  mean(rottcif$cif.fin.b-rottcif$cif.fin.c )
  
  m <- length(tgrid)
  print(c(sum(rottcif$cif.fin.a[m,]),sum(rottcif$cif.fin.b[m,]),sum(rottcif$cif.fin.c[m,])))
  
  table(dat$delta)
  
  ## Method 3
  r1 = cbind(tgrid,rottcif$cif.fin.a[,1],rottcif$cif.fin.a.lo[,1],rottcif$cif.fin.a.hi[,1],
             rottcif$cif.fin.a[,2],rottcif$cif.fin.a.lo[,2],rottcif$cif.fin.a.hi[,2])         
  colnames(r1) = c('time','cif.a.1','ci.a.lo.1','ci.a.hi.1','cif.a.2','ci.a.lo.2','ci.a.hi.2')
  
  cif.data <- data.frame(r1)
  cif.data$anyevent <- cif.data$cif.a.1+cif.data$cif.a.2
  
  p.prop <-  ggplot() + theme_classic() + xlab("Seconds") +ylab("CIF")  +
    scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1)) +
    geom_step(data=cif.data, mapping=aes(x=time, y=cif.a.1,color="Correct")) +
    geom_step(data=cif.data, mapping=aes(x=time, y=cif.a.2,color="Incorrect")) +
    geom_step(data=cif.data, mapping=aes(x=time, y=anyevent,color="Total")) 
  
  ## Method 1
  r1 = cbind(tgrid,rottcif$cif.fin.b[,1],rottcif$cif.fin.b.lo[,1],rottcif$cif.fin.b.hi[,1],
             rottcif$cif.fin.b[,2],rottcif$cif.fin.b.lo[,2],rottcif$cif.fin.b.hi[,2])         
  colnames(r1) = c('time','cif.b.1','ci.b.lo.1','ci.b.hi.1',
                   'cif.b.2','ci.b.lo.2','ci.b.hi.2')
  
  cif.data <- data.frame(r1)
  cif.data$anyevent <- cif.data$cif.b.1+cif.data$cif.b.2
  
  p.asaf <- p.prop + 
    geom_step(data=cif.data, mapping=aes(x=time, y=cif.b.1,color="Correct"), linetype = "dashed", size = 1) +
    geom_step(data=cif.data, mapping=aes(x=time, y=cif.b.2,color="Incorrect"), linetype = "dashed", size = 1) +
    geom_step(data=cif.data, mapping=aes(x=time, y=anyevent,color="Total"), linetype = "dashed", size = 1) 
  
  ## Method 2
  r1 = cbind(tgrid,rottcif$cif.fin.c[,1],rottcif$cif.fin.c.lo[,1],rottcif$cif.fin.c.hi[,1],
             rottcif$cif.fin.c[,2],rottcif$cif.fin.c.lo[,2],rottcif$cif.fin.c.hi[,2])         
  colnames(r1) = c('time','cif.c.1','ci.c.lo.1','ci.c.hi.1',
                   'cif.c.2','ci.c.lo.2','ci.c.hi.2')
  
  cif.data <- data.frame(r1)
  cif.data$anyevent <- cif.data$cif.c.1+cif.data$cif.c.2
  
  p.book <-  p.asaf + 
    geom_step(data=cif.data, mapping=aes(x=time, y=cif.c.1,color="Correct"), linetype="dotted", size = 1) +
    geom_step(data=cif.data, mapping=aes(x=time, y=cif.c.2,color="Incorrect"), linetype="dotted", size = 1) +
    geom_step(data=cif.data, mapping=aes(x=time, y=anyevent,color="Total"), linetype="dotted", size = 1) +
    geom_hline(yintercept=1 ,alpha=0.2) +
    theme(legend.position='none')
  
  print(i)
  assign(paste0("variable_", i),p.book)
}

################
#Combining plots
################
library(multipanelfigure)

figure103A <- multi_panel_figure(width = 210, height = 240, columns = 2, rows = 2, panel_label_type = "none")

figure103A %<>%
  fill_panel(variable_1, column = 1, row = 1) %<>%
  fill_panel(variable_2, column = 1, row = 2) %<>%
  fill_panel(variable_5, column = 2, row = 1) %<>%
  fill_panel(variable_6, column = 2, row = 2) 
figure103A

figure103B <- multi_panel_figure(width = 210, height = 240, columns = 2, rows = 2, panel_label_type = "none")

figure103B %<>%
  fill_panel(variable_3, column = 1, row = 1) %<>%
  fill_panel(variable_4, column = 1, row = 2) %<>%
  fill_panel(variable_7, column = 2, row = 1) %<>%
  fill_panel(variable_8, column = 2, row = 2) 
figure103B















