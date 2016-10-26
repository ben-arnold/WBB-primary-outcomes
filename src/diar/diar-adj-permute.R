
#---------------------------------------
# diar-adj-permute.R
#
# ben arnold (benarnold@berkeley.edu)
#
# adjusted permutation tests for
# differences in diarrhea prevalence
# for the year 1 and year 2 
# visits combined
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(plyr)
library(coin)
library(SuperLearner)

# source the base functions
# which includes the permutation test function used below
source("~/WBBpa/src/basefns/washb-base-functions.R")


#---------------------------------------
# Load the analysis dataset,
# the baseline covariate dataset
#---------------------------------------

bd <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-enrol.csv")

d <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-diar.csv")

# merge the baseline dataset to the follow-up dataset
ad <- merge(bd,d,by=c("dataid","clusterid","block","tr"),all.x=F,all.y=T)
dim(d)
dim(ad)

ad$block <- as.factor(ad$block)

# ensure that month is coded as a factor
ad$month <- factor(ad$month)

#---------------------------------------
# Subset the Data to Follow-up data only
#---------------------------------------
table(ad$svy)
ad <- subset(ad,svy>0)

#---------------------------------------
# Exclude:
# * siblings who were born after enrollment
# * siblings who were >36 mos at enrollment
# * children with missing outcome data
#---------------------------------------
table(ad$sibnewbirth)
table(ad$gt36mos)
table(is.na(ad$diar7d))

ad <- subset(ad,sibnewbirth==0)
dim(ad)

table(ad$gt36mos)
ad <- subset(ad,gt36mos==0)
dim(ad)

table(is.na(ad$diar7d))
ad <- subset(ad,!is.na(ad$diar7d))
dim(ad)

# ensure that month is coded as a factor
ad$month <- factor(ad$month)

# lump field investigators with <100 measurements into a single code (arbitrarily choosing N08002 from the list). This FRA code list with N<100 is the same for years 1 and 2.
ad$fracode[ad$fracode=="20"    ] <- "N08002"
ad$fracode[ad$fracode=="D03538"] <- "N08002"
ad$fracode[ad$fracode=="N00436"] <- "N08002"
ad$fracode[ad$fracode=="N05267"] <- "N08002"
ad$fracode[ad$fracode=="N05268"] <- "N08002"
ad$fracode[ad$fracode=="N05271"] <- "N08002"
ad$fracode[ad$fracode=="N06485"] <- "N08002"
ad$fracode[ad$fracode=="N08000"] <- "N08002"
ad$fracode <- factor(ad$fracode)

# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Water","Sanitation","Handwashing","Nutrition","WSH","Nutrition + WSH","Control"))

# sort the data for perfect replication with jade on the V-fold cross-validation
ad <- ad[order(ad$block,ad$clusterid,ad$dataid,ad$childid),]

#---------------------------------------
# Permutation tests
#---------------------------------------

# excluding birthord (only measured in index children in year 2)

# pre-specified covariates (but not treatment)
Ws <- subset(ad,select=c("fracode","month","agedays","sex","momage","momedu","momheight","hfiacat","Nlt18","Ncomp","watmin","elec","floor","walls","roof","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile"))

# restrict to complete cases
SLd <- data.frame(id=ad$clusterid,block=ad$block,tr=ad$tr,Y=ad$diar7d,Ws)
dim(SLd)
SLd <- SLd[complete.cases(SLd),]
dim(SLd)


# pre-screen the covariates for those associated with the outcome (LR test P<0.2)
# see washb_prescreen() and design_matrix() in the base functions
Wscreen <- washb_prescreen(Y=SLd$Y,Ws=SLd[,5:ncol(SLd)],family="binomial",pval=0.2)
Wselect <- subset(SLd,select=Wscreen)
Wselect <- design_matrix(Wselect)


# algorithmic fit of the outcome as a function of selected covariates
set.seed(12345)
SLfit <- SuperLearner(Y=SLd$Y,X=Wselect,id=SLd$id,
                       family="binomial",
                       SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet")
)
SLfit
SLd$pY <- as.vector(predict(SLfit)$pred)
SLd$r <- SLd$Y-SLd$pY


# Hypothesis 1 permutation tests
permute.C.W <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Water"),seed=12345)
permute.C.S <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Sanitation"),seed=12345)
permute.C.H <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Handwashing"),seed=12345)
permute.C.WSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","WSH"),seed=12345)
permute.C.N   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Nutrition"),seed=12345)
permute.C.NWSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Nutrition + WSH"),seed=12345)

# Hypothesis 2 permutation tests
permute.W.WSH   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,c("Water","WSH"),seed=12345)
permute.S.WSH   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,c("Sanitation","WSH"),seed=12345)
permute.H.WSH   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,c("Handwashing","WSH"),seed=12345)

#---------------------------------------
# put objects in the standard format
#---------------------------------------
h1res <- list(permute.C.W,permute.C.S,permute.C.H,permute.C.WSH,permute.C.N,permute.C.NWSH)
diar_h1_pval_adj <- as.matrix(sapply(h1res,function(x) x$p.value),nrow=6)
rownames(diar_h1_pval_adj) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

h2res <- list(permute.W.WSH,permute.S.WSH,permute.H.WSH)
diar_h2_pval_adj <- as.matrix(sapply(h2res,function(x) x$p.value),nrow=3)
rownames(diar_h2_pval_adj) <- c("WSH v W","WSH v S","WSH v H")


#---------------------------------------
# print results
#---------------------------------------
diar_h1_pval_adj

diar_h2_pval_adj


#---------------------------------------
# add suffix for replication
#---------------------------------------
diar_h1_pval_adj_b <- diar_h1_pval_adj
diar_h2_pval_adj_b <- diar_h2_pval_adj
rm(diar_h1_pval_adj,diar_h2_pval_adj)

#---------------------------------------
# save all of the results
# excluding the datasets
#---------------------------------------
rm(bd,d,ad)
rm(SLfit,h1res,h2res)
save.image("~/dropbox/WBB-primary-analysis/results/raw/ben/bangladesh-diar-adj-permute.RData")


