
#---------------------------------------
# uptake.R
#
# ben arnold (benarnold@berkeley.edu)
#
# summarize measures of uptake / compliance
# by study arm and measurement round
# (enrollment, year 1, year 2)
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-uptake.csv
#
# output files:
# bangladesh-uptake.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(dplyr)

# source the base functions
source("~/WBBpa/src/basefns/washb-base-functions.R")

#---------------------------------------
# load the uptake analysis dataset
#---------------------------------------
d <- read.csv("~/dropbox/wbb-primary-analysis/data/final/ben/washb-bangladesh-uptake.csv")

# re-order the treatment factor for convenience
d$tr <- factor(d$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

#---------------------------------------
# for each uptake indicator, summarize
# the number of obs and the % at each
# measurement round
#---------------------------------------

d.svy <- group_by(d, tr,svy)

# number of observations
ncompounds <- summarise(d.svy,n=n())

# store water
storewat <- summarise(d.svy,mean=mean(storewat,na.rm=T))

# store water with detectable chlorine
freechl <- summarise(d.svy,mean=mean(freechl,na.rm=T))

# Latrine w/ a functional water seal
latseal <- summarise(d.svy,mean=mean(latseal,na.rm=T))

# No visible feces on the latrine slab or floor
latfeces <- summarise(d.svy,mean=mean(latfeces,na.rm=T))

# No human feces in the compound
humfeces <- summarise(d.svy,mean=mean(humfeces,na.rm=T))

# Primary handwashing station has water
hwsw <- summarise(d.svy,mean=mean(hwsw,na.rm=T))

# Primary handwashing station has soap
hwss <- summarise(d.svy,mean=mean(hwss,na.rm=T))

# Primary handwashing station has soap & water
hwsws <- summarise(d.svy,mean=mean(hwsws,na.rm=T))

# Mean sachets of LNS fed in prior week to index child 6-24 mos
rlnsp <- summarise(d.svy,mean=mean(rlnsp,na.rm=T))

#---------------------------------------
# combine estimates into a single matrix
#---------------------------------------
uptake.tab <- as.data.frame(
  rbind(
  ncompounds$n,
  storewat$mean,
  freechl$mean,
  latseal$mean,
  latfeces$mean,
  humfeces$mean,
  hwsw$mean,
  hwss$mean,
  hwsws$mean,
  rlnsp$mean
))
names(uptake.tab) <- paste(rep(levels(d$tr),rep(3,7)),c("0","1","2"))
uptake.tab$label=c(
  "N compounds",
  "Store water",
  "Store water with detectable free chlorine",
  "Latrine with functional water seal",
  "Visible feces on latrine slab or floor",
  "Human feces in house or compound",
  "Primary handwashing station has water",
  "Primary handwashing station has soap",
  "Primary handwashing station has soap & water",
  "LNS sachet consumption"
)
# reorder label
uptake_table_b <- uptake.tab[,c(ncol(uptake.tab),1:(ncol(uptake.tab)-1))]

# print table
uptake_table_b

#---------------------------------------
# Calculate means and influence-curve 
# based 95% CIs by survey round
#---------------------------------------
arms <- levels(d$tr)

# store water
storewat0 <- sapply(arms,tmle.mean.est,Y=d$storewat,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
storewat1 <- sapply(arms,tmle.mean.est,Y=d$storewat,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
storewat2 <- sapply(arms,tmle.mean.est,Y=d$storewat,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# store water with detectable chlorine (not measured at enrollment)
# note: since there were no events in any of the non-water arms (actually, 1 in HW)
# we cannot estimate 95% CIs. Will pad the matrix with zeros and NAs
# freechl0 <- sapply(arms,tmle.mean.est,Y=d$freechl,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
warms <- c("Water","WSH","Nutrition + WSH")
freechl1 <- sapply(warms,tmle.mean.est,Y=d$freechl,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
  # pad with zeros and missings for other treatment arms
  freechl1 <- cbind(c(0,NA,NA),freechl1[,1],c(0,NA,NA),c(0,NA,NA),freechl1[,2],c(0,NA,NA),freechl1[,3])
  colnames(freechl1) <- colnames(storewat1)
  
freechl2 <- sapply(warms,tmle.mean.est,Y=d$freechl,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)
  # pad with zeros and missings for other treatment arms
  freechl2 <- cbind(c(0,NA,NA),freechl2[,1],c(0,NA,NA),c(0,NA,NA),freechl2[,2],c(0,NA,NA),freechl2[,3])
  colnames(freechl2) <- colnames(storewat2)

# Latrine w/ a functional water seal
latseal0 <- sapply(arms,tmle.mean.est,Y=d$latseal,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
latseal1 <- sapply(arms,tmle.mean.est,Y=d$latseal,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
latseal2 <- sapply(arms,tmle.mean.est,Y=d$latseal,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Visible feces on the latrine slab or floor
latfeces0 <- sapply(arms,tmle.mean.est,Y=d$latfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
latfeces1 <- sapply(arms,tmle.mean.est,Y=d$latfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
latfeces2 <- sapply(arms,tmle.mean.est,Y=d$latfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Human feces in the compound
humfeces0 <- sapply(arms,tmle.mean.est,Y=d$humfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
humfeces1 <- sapply(arms,tmle.mean.est,Y=d$humfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
humfeces2 <- sapply(arms,tmle.mean.est,Y=d$humfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Primary handwashing station has water
hwsw0 <- sapply(arms,tmle.mean.est,Y=d$hwsw,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
hwsw1 <- sapply(arms,tmle.mean.est,Y=d$hwsw,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
hwsw2 <- sapply(arms,tmle.mean.est,Y=d$hwsw,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Primary handwashing station has soap
hwss0 <- sapply(arms,tmle.mean.est,Y=d$hwss,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
hwss1 <- sapply(arms,tmle.mean.est,Y=d$hwss,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
hwss2 <- sapply(arms,tmle.mean.est,Y=d$hwss,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Primary handwashing station has soap & water
hwsws0 <- sapply(arms,tmle.mean.est,Y=d$hwsws,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
hwsws1 <- sapply(arms,tmle.mean.est,Y=d$hwsws,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
hwsws2 <- sapply(arms,tmle.mean.est,Y=d$hwsws,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Mean sachets of LNS fed in prior week to index child 6-24 mos (not measured at enrollment, only measured in nutrition arms)
narms <- arms[grep("Nutrition",arms)]
rlnsp1 <- sapply(narms,tmle.mean.est,Y=d$rlnsp,tr=d$tr,svy=d$svy,id=d$clusterid,s=1,family="gaussian")
  # pad with missings for other treatment arms
  rlnsp1 <- cbind(matrix(NA,nrow=3,ncol=5),rlnsp1)
  colnames(rlnsp1) <- colnames(hwsws1)

rlnsp2 <- sapply(narms,tmle.mean.est,Y=d$rlnsp,tr=d$tr,svy=d$svy,id=d$clusterid,s=2,family="gaussian")
  # pad with missings for other treatment arms
  rlnsp2 <- cbind(matrix(NA,nrow=3,ncol=5),rlnsp2)
  colnames(rlnsp2) <- colnames(hwsws2)



#---------------------------------------
# save the objects
#---------------------------------------
rm(d)
save.image(file="~/dropbox/wbb-primary-analysis/results/raw/ben/bangladesh-uptake.RData")




