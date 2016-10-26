
#---------------------------------------
# hcz-unadj-permute.R
#
# ben arnold (benarnold@berkeley.edu)
#
# unadjusted permutation tests for
# differences in hcz at year 1 and year 2
# follow-ups
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(plyr)
library(coin)

# source the base functions
# which includes the permutation test function used below
source("~/WBBpa/src/basefns/washb-base-functions.R")


#---------------------------------------
# load the anthropometry analysis data
#---------------------------------------
d <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-anthro.csv")

d$block <- as.factor(d$block)

# subset the anthropometry to target children (excluding siblings)
dim(d)
d <- subset(d,tchild=="Target child")
dim(d)

# Drop children with extreme hcz values
table(d$hcz_x)
d <- subset(d,hcz_x!=1)

# Exclude children with missing data (none)
table(is.na(d$hcz))

#---------------------------------------
# YEAR 1
#---------------------------------------
# subset to the relevant measurement
table(d$svy)
ad <- subset(d,svy==1)
dim(ad)

# Hypothesis 1 permutation tests
set.seed(2342)
permute.1.C.W <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Water"))
permute.1.C.S <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Sanitation"))
permute.1.C.H <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Handwashing"))
permute.1.C.WSH <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","WSH"))
permute.1.C.N   <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Nutrition"))
permute.1.C.NWSH <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Nutrition + WSH"))

# Hypothesis 3 permutation tests
set.seed(4524)
permute.1.N.WSHN   <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Nutrition","Nutrition + WSH"))
permute.1.WSH.NWSH <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("WSH","Nutrition + WSH"))

#---------------------------------------
# YEAR 2
#---------------------------------------
# subset to the relevant measurement
table(d$svy)
ad <- subset(d,svy==2)
dim(ad)

# Hypothesis 1 permutation tests
set.seed(523423)
permute.2.C.W <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Water"))
permute.2.C.S <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Sanitation"))
permute.2.C.H <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Handwashing"))
permute.2.C.WSH <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","WSH"))
permute.2.C.N   <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Nutrition"))
permute.2.C.NWSH <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Control","Nutrition + WSH"))

# Hypothesis 3 permutation tests
set.seed(23624)
permute.2.N.WSHN   <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("Nutrition","Nutrition + WSH"))
permute.2.WSH.NWSH <- washb.permute(Y=ad$hcz,tr=ad$tr,block=ad$block,c("WSH","Nutrition + WSH"))


#---------------------------------------
# put objects in the standard format
#---------------------------------------
t1h1res <- list(permute.1.C.W,permute.1.C.S,permute.1.C.H,permute.1.C.WSH,permute.1.C.N,permute.1.C.NWSH)
hcz_t1_h1_pval_unadj <- as.matrix(sapply(t1h1res,function(x) x$p.value),nrow=6)
rownames(hcz_t1_h1_pval_unadj) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

t1h3res <- list(permute.1.N.WSHN,permute.1.WSH.NWSH)
hcz_t1_h3_pval_unadj <- as.matrix(sapply(t1h3res,function(x) x$p.value),nrow=6)
rownames(hcz_t1_h3_pval_unadj) <- c("Nutrition + WSH v Nutrition","Nutrition + WSH v WSH")

t2h1res <- list(permute.2.C.W,permute.2.C.S,permute.2.C.H,permute.2.C.WSH,permute.2.C.N,permute.2.C.NWSH)
hcz_t2_h1_pval_unadj <- as.matrix(sapply(t2h1res,function(x) x$p.value),nrow=6)
rownames(hcz_t2_h1_pval_unadj) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

t2h3res <- list(permute.2.N.WSHN,permute.2.WSH.NWSH)
hcz_t2_h3_pval_unadj <- as.matrix(sapply(t2h3res,function(x) x$p.value),nrow=6)
rownames(hcz_t2_h3_pval_unadj) <-c("Nutrition + WSH v Nutrition","Nutrition + WSH v WSH")


#---------------------------------------
# print results
#---------------------------------------
hcz_t1_h1_pval_unadj

hcz_t1_h3_pval_unadj

hcz_t2_h1_pval_unadj

hcz_t2_h3_pval_unadj

#---------------------------------------
# add suffix for replication
#---------------------------------------
hcz_t1_h1_pval_unadj_b <- hcz_t1_h1_pval_unadj
hcz_t1_h3_pval_unadj_b <- hcz_t1_h3_pval_unadj
hcz_t2_h1_pval_unadj_b <- hcz_t2_h1_pval_unadj
hcz_t2_h3_pval_unadj_b <- hcz_t2_h3_pval_unadj
rm(hcz_t1_h1_pval_unadj,hcz_t1_h3_pval_unadj,hcz_t2_h1_pval_unadj,hcz_t2_h3_pval_unadj)

#---------------------------------------
# save all of the results
# excluding the datasets
#---------------------------------------
rm(d,ad)
rm(t1h1res,t1h3res,t2h1res,t2h3res)
save.image("~/dropbox/WBB-primary-analysis/results/raw/ben/bangladesh-hcz-unadj-permute.RData")


