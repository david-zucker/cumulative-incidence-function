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
znew = c(1 , 35 , 0  , 5 ) 
#max.ev.tim = max(dat$time[which(dat$delta>0)])
tgrid = sort(dat$time)
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
cif.data$ci.a.lo.1 <- ifelse(cif.data$ci.a.lo.1<0,0,cif.data$ci.a.lo.1)
cif.data$ci.a.lo.2 <- ifelse(cif.data$ci.a.lo.2<0,0,cif.data$ci.a.lo.2)

p.prop <-  ggplot() + theme_classic() + xlab("Seconds") +ylab("CIF") + 
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1)) +
  geom_step(data=cif.data, mapping=aes(x=time, y=cif.a.1,color="Correct")) +
  geom_step(data=cif.data, mapping=aes(x=time, y=cif.a.2,color="Incorrect")) +
  geom_step(data=cif.data, mapping=aes(x=time, y=anyevent,color="Total")) 

p.prop <- p.prop + geom_hline(yintercept=1) +
  geom_ribbon(data=cif.data,aes(x=time, ymin = ci.a.lo.1, ymax = ci.a.hi.1),alpha=0.1) +
  geom_ribbon(data=cif.data,aes(x=time, ymin = ci.a.lo.2, ymax = ci.a.hi.2),alpha=0.1) +
  ggtitle("Method 3") +
  theme(legend.position='none')

## Method 1
r1 = cbind(tgrid,rottcif$cif.fin.b[,1],rottcif$cif.fin.b.lo[,1],rottcif$cif.fin.b.hi[,1],
           rottcif$cif.fin.b[,2],rottcif$cif.fin.b.lo[,2],rottcif$cif.fin.b.hi[,2])         
colnames(r1) = c('time','cif.b.1','ci.b.lo.1','ci.b.hi.1',
                 'cif.b.2','ci.b.lo.2','ci.b.hi.2')


cif.data <- data.frame(r1)
cif.data$anyevent <- cif.data$cif.b.1+cif.data$cif.b.2
cif.data$ci.b.lo.1 <- ifelse(cif.data$ci.b.lo.1<0,0,cif.data$ci.b.lo.1)
cif.data$ci.b.lo.2 <- ifelse(cif.data$ci.b.lo.2<0,0,cif.data$ci.b.lo.2)

p.asaf <-  ggplot() + theme_classic() + xlab("Seconds") +ylab("CIF") +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1)) +
  geom_step(data=cif.data, mapping=aes(x=time, y=cif.b.1,color="Correct")) +
  geom_step(data=cif.data, mapping=aes(x=time, y=cif.b.2,color="Incorrect")) +
  geom_step(data=cif.data, mapping=aes(x=time, y=anyevent,color="Total")) 

p.asaf <- p.asaf + geom_hline(yintercept=1) +
  geom_ribbon(data=cif.data,aes(x=time, ymin = ci.b.lo.1, ymax = ci.b.hi.1),alpha=0.1) +
  geom_ribbon(data=cif.data,aes(x=time, ymin = ci.b.lo.2, ymax = ci.b.hi.2),alpha=0.1) +
  ggtitle("Method 1") +
  theme(legend.position='none')


## Method 2
r1 = cbind(tgrid,rottcif$cif.fin.c[,1],rottcif$cif.fin.c.lo[,1],rottcif$cif.fin.c.hi[,1],
           rottcif$cif.fin.c[,2],rottcif$cif.fin.c.lo[,2],rottcif$cif.fin.c.hi[,2])         
colnames(r1) = c('time','cif.c.1','ci.c.lo.1','ci.c.hi.1',
                 'cif.c.2','ci.c.lo.2','ci.c.hi.2')


cif.data <- data.frame(r1)
cif.data$anyevent <- cif.data$cif.c.1+cif.data$cif.c.2
cif.data$ci.c.lo.1 <- ifelse(cif.data$ci.c.lo.1<0,0,cif.data$ci.c.lo.1)
cif.data$ci.c.lo.2 <- ifelse(cif.data$ci.c.lo.2<0,0,cif.data$ci.c.lo.2)

p.book <-  ggplot() + theme_classic() + xlab("Seconds") +ylab("CIF") + 
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1)) +
  geom_step(data=cif.data, mapping=aes(x=time, y=cif.c.1,color="Correct")) +
  geom_step(data=cif.data, mapping=aes(x=time, y=cif.c.2,color="Incorrect")) +
  geom_step(data=cif.data, mapping=aes(x=time, y=anyevent,color="Total")) 

p.book <- p.book + geom_hline(yintercept=1) +
  geom_ribbon(data=cif.data,aes(x=time, ymin = ci.c.lo.1, ymax = ci.c.hi.1),alpha=0.1) +
  geom_ribbon(data=cif.data,aes(x=time, ymin = ci.c.lo.2, ymax = ci.c.hi.2),alpha=0.1) +
  ggtitle("Method 2") +
  theme(legend.position='none')

################
#Combining plots
################
library(multipanelfigure)

figure3 <- multi_panel_figure(width = 120, height = 180, columns = 1, rows = 3, panel_label_type = "none")

figure3 %<>%
  fill_panel(p.asaf, column = 1, row = 1) %<>%
  fill_panel(p.book, column = 1, row = 2) %<>%
  fill_panel(p.prop, column = 1, row = 3) 
figure3



