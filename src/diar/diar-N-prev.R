
#---------------------------------------
# diar-N-prev.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate the number of diarrhea
# cases and the prevalence by
# survey round and treatment arm
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-diar-public.csv
#
# output files:
#	bangladesh-diar-N-prev-ben.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls()); library(here)
library(tidyverse)
library(tmle)

# source the base functions
source("src/basefns/washb-base-functions.R")


#---------------------------------------
# Load the analysis dataset
#---------------------------------------

d <- read.csv(here("data/washb-bangladesh-diar-public.csv"))

# merge in the treatment assignments
d_tr    <- read.csv(here('data/washb-bangladesh-tr-public.csv'))
d <- left_join(d,d_tr,by=c("clusterid","block"))


#---------------------------------------
# Exclude:
# * siblings who were born after enrollment
# * siblings who were >36 mos at enrollment
# * children with missing outcome data
#---------------------------------------

table(d$svy,d$sibnewbirth)
table(d$svy,d$gt36mos)
table(d$svy,is.na(d$diar7d))


table(d$sibnewbirth)
ad <- subset(d,sibnewbirth==0)
dim(ad)

table(ad$gt36mos)
ad <- subset(ad,gt36mos==0)
dim(ad)

table(is.na(ad$diar7d))
ad <- subset(ad,!is.na(ad$diar7d))


# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

#---------------------------------------
# Calculate N and 7-day diarrhea cases
# by survey round and by treatment arm
#---------------------------------------

Nchild <- tapply(ad$diar7d,list(ad$svy,ad$tr),function(x) length(x))
diar7d <- tapply(ad$diar7d,list(ad$svy,ad$tr),function(x) sum(x))

#---------------------------------------
# calculate Ns and prevalence for
# all of follow-up (surveys 1 and 2)
#
# Lancet is requesting we include the
# SD of the prevalence for the manuscript
# table, so add that
#---------------------------------------
Nchildfu <- tapply(ad$diar7d[ad$svy==1|ad$svy==2],ad$tr[ad$svy==1|ad$svy==2],function(x) length(x))
diar7dfu <- tapply(ad$diar7d[ad$svy==1|ad$svy==2],ad$tr[ad$svy==1|ad$svy==2],function(x) sum(x))
diar7dfu_sd <- tapply(ad$diar7d[ad$svy==1|ad$svy==2],ad$tr[ad$svy==1|ad$svy==2],function(x) sd(x))

#---------------------------------------
# Create randomization block indicators
# (need to condition on these for 
# clusters to be indepdent)  The tmle()
# function can't handle factors, so this
# is necessary
#
# NOTE: tmle() returns identical CIs
# irrespective of whether we include
# these indicators
# might want to double-check the
# estimates using Stata as well?
#---------------------------------------

# create indicators for the randomization strata
blocks <- model.matrix(~as.factor(ad$block))[,-c(1)]
colnames(blocks) <- paste("block",2:90,sep="")

#---------------------------------------
# Calculate unadjusted prevalences
# and 95% CIs by arm
#---------------------------------------

# Calculate means and influence-curve based 95% CIs by survey round
set.seed(12345)
arms <- levels(ad$tr)
mu0 <- sapply(arms,tmle.mean.est,Y=ad$diar7d,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=0)
mu1 <- sapply(arms,tmle.mean.est,Y=ad$diar7d,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=1)
mu2 <- sapply(arms,tmle.mean.est,Y=ad$diar7d,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=2)

# compute means by all follow-up rounds
ad$svy2 <- ifelse(ad$svy==1|ad$svy==2,1,0)
mu12 <- sapply(arms,tmle.mean.est,Y=ad$diar7d,tr=ad$tr,svy=ad$svy2,id=ad$clusterid,s=1)

#---------------------------------------
# Create final objects
#---------------------------------------

# condense results into pre-specfied objects
diar_t0_n <- cbind(Nchild[1,],diar7d[1,])
diar_t1_n <- cbind(Nchild[2,],diar7d[2,])
diar_t2_n <- cbind(Nchild[3,],diar7d[3,])
diar_t12_n <- cbind(Nchildfu,diar7dfu)
colnames(diar_t0_n) <- colnames(diar_t1_n) <-  colnames(diar_t2_n) <- colnames(diar_t12_n) <- c("N","n")

diar_t0_prev <- t(mu0)
diar_t1_prev <- t(mu1)
diar_t2_prev <- t(mu2)
diar_t12_prev <- t(mu12)

diar_t12_prev_sd <- diar7dfu_sd

#---------------------------------------
# Internal consistency check:
# the n/N and tmle prevalence estimates
# should be the same
#---------------------------------------
round(cbind(diar_t0_n[,2]/diar_t0_n[,1],diar_t0_prev[,1]),4)
round(cbind(diar_t1_n[,2]/diar_t1_n[,1],diar_t1_prev[,1]),4)
round(cbind(diar_t2_n[,2]/diar_t2_n[,1],diar_t2_prev[,1]),4)

round(cbind(diar_t12_n[,2]/diar_t12_n[,1],diar_t12_prev[,1]),4)

#---------------------------------------
# Print and save results
#---------------------------------------
diar_t0_n
diar_t1_n
diar_t2_n
diar_t12_n
round(diar_t0_prev,4)
round(diar_t1_prev,4)
round(diar_t2_prev,4)
round(diar_t12_prev,4)

round(diar_t12_prev_sd,4)


# add 'b' suffix for comparison w/ jade
# diar_t0_n_b <- diar_t0_n
# diar_t1_n_b <- diar_t1_n
# diar_t2_n_b <- diar_t2_n
# diar_t12_n_b <- diar_t12_n
# diar_t0_prev_b <- diar_t0_prev
# diar_t1_prev_b <- diar_t1_prev
# diar_t2_prev_b <- diar_t2_prev
# diar_t12_prev_b <- diar_t12_prev
# diar_t12_prev_sd_b <- diar_t12_prev_sd

save(diar_t0_n_b,diar_t1_n_b,diar_t2_n_b,diar_t12_n_b,diar_t0_prev_b,diar_t1_prev_b,diar_t2_prev_b,diar_t12_prev_b,diar_t12_prev_sd_b,file="results/bangladesh-diar-N-prev.RData")





	