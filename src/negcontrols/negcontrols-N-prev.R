
#---------------------------------------
# negcontrols-N-prev.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate the number of bruise
# cases and the prevalence by
# survey round and treatment arm
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-diar-public.csv
#
# output files:
#	bangladesh-negcontrols-N-prev.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
source(here::here("src/0-config.R"))


# source the base functions
source(here("src/basefns/washb-base-functions.R"))


#---------------------------------------
# Load the analysis dataset
#---------------------------------------

d <- read.csv(here("data/washb-bangladesh-diar-public.csv"))

# merge in the treatment assignments
tr    <- read.csv(here("data/washb-bangladesh-tr-public.csv"))
d <- left_join(d,tr,by=c("clusterid","block"))

#---------------------------------------
# Exclude:
# * siblings who were born after enrollment
# * siblings who were >36 mos at enrollment
# * children with missing outcome data
#---------------------------------------

table(d$svy,d$sibnewbirth)
table(d$svy,d$gt36mos)
table(d$svy,is.na(d$bruise7d))


table(d$sibnewbirth)
ad <- subset(d,sibnewbirth==0)
dim(ad)

table(ad$gt36mos)
ad <- subset(ad,gt36mos==0)
dim(ad)

table(is.na(ad$bruise7d))
ad <- subset(ad,!is.na(ad$bruise7d))


# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

#---------------------------------------
# Calculate N and 7-day diarrhea cases
# by survey round and by treatment arm
#---------------------------------------

Nchild <- tapply(ad$bruise7d,list(ad$svy,ad$tr),function(x) length(x))
bruise7d <- tapply(ad$bruise7d,list(ad$svy,ad$tr),function(x) sum(x))

#---------------------------------------
# calculate Ns and prevalence for
# all of follow-up (surveys 1 and 2)
#---------------------------------------
Nchildfu <- tapply(ad$bruise7d[ad$svy==1|ad$svy==2],ad$tr[ad$svy==1|ad$svy==2],function(x) length(x))
bruise7dfu <- tapply(ad$bruise7d[ad$svy==1|ad$svy==2],ad$tr[ad$svy==1|ad$svy==2],function(x) sum(x))


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
mu0 <- sapply(arms,tmle.mean.est,Y=ad$bruise7d,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=0)
mu1 <- sapply(arms,tmle.mean.est,Y=ad$bruise7d,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=1)
mu2 <- sapply(arms,tmle.mean.est,Y=ad$bruise7d,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=2)

# compute means by all follow-up rounds
ad$svy2 <- ifelse(ad$svy==1|ad$svy==2,1,0)
mu12 <- sapply(arms,tmle.mean.est,Y=ad$bruise7d,tr=ad$tr,svy=ad$svy2,id=ad$clusterid,s=1)

#---------------------------------------
# Create final objects
#---------------------------------------

# condense results into pre-specfied objects
bruise_t0_n <- cbind(Nchild[1,],bruise7d[1,])
bruise_t1_n <- cbind(Nchild[2,],bruise7d[2,])
bruise_t2_n <- cbind(Nchild[3,],bruise7d[3,])
bruise_t12_n <- cbind(Nchildfu,bruise7dfu)
colnames(bruise_t0_n) <- colnames(bruise_t1_n) <-  colnames(bruise_t2_n) <- colnames(bruise_t12_n) <- c("N","n")

bruise_t0_prev <- t(mu0)
bruise_t1_prev <- t(mu1)
bruise_t2_prev <- t(mu2)
bruise_t12_prev <- t(mu12)



#---------------------------------------
# Internal consistency check:
# the n/N and tmle prevalence estimates
# should be the same
#---------------------------------------
round(cbind(bruise_t0_n[,2]/bruise_t0_n[,1],bruise_t0_prev[,1]),4)
round(cbind(bruise_t1_n[,2]/bruise_t1_n[,1],bruise_t1_prev[,1]),4)
round(cbind(bruise_t2_n[,2]/bruise_t2_n[,1],bruise_t2_prev[,1]),4)

round(cbind(bruise_t12_n[,2]/bruise_t12_n[,1],bruise_t12_prev[,1]),4)

#---------------------------------------
# Print and save results
#---------------------------------------
bruise_t0_n
bruise_t1_n
bruise_t2_n
bruise_t12_n
round(bruise_t0_prev,4)
round(bruise_t1_prev,4)
round(bruise_t2_prev,4)
round(bruise_t12_prev,4)

# add 'b' suffix for comparison w/ jade
# bruise_t0_n_b <- bruise_t0_n
# bruise_t1_n_b <- bruise_t1_n
# bruise_t2_n_b <- bruise_t2_n
# bruise_t12_n_b <- bruise_t12_n
# bruise_t0_prev_b <- bruise_t0_prev
# bruise_t1_prev_b <- bruise_t1_prev
# bruise_t2_prev_b <- bruise_t2_prev
# bruise_t12_prev_b <- bruise_t12_prev


save(bruise_t0_n,bruise_t1_n,bruise_t2_n,bruise_t12_n,bruise_t0_prev,bruise_t1_prev,bruise_t2_prev,bruise_t12_prev,file=here("results/bangladesh-bruise-N-prev.RData"))





	