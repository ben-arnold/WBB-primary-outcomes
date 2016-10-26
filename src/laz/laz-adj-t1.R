
#---------------------------------------
# laz-diff-adj-t1.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate unadjusted differences
# between treatment arms for H1 and H3
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-anthro.csv
#	washb-bangladesh-enrol.csv
#
# output files:
#	bangladesh-laz-adj-t1-ben.RData
#
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(tmle)
library(SuperLearner)

# source the base functions
source("~/WBBpa/src/basefns/washb-base-functions.R")


#---------------------------------------
# Load the analysis dataset,
# the baseline covariate dataset
#---------------------------------------

bd <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-enrol.csv",colClasses=c("dataid"="character"))

d <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-anthro.csv",colClasses=c("dataid"="character"))

# merge the baseline dataset to the follow-up dataset
ad <- merge(bd,d,by=c("dataid","clusterid","block","tr"),all.x=T,all.y=T)
dim(d)
dim(ad)

#---------------------------------------
# subset to the relevant measurement
# Year 1 or Year 2
#---------------------------------------
table(ad$svy)
ad <- subset(ad,svy==1)
dim(ad)

# subset the anthropometry to target children (excluding siblings)
dim(ad)
ad <- subset(ad,tchild=="Target child")
dim(ad)

# Drop children with extreme LAZ values
table(ad$laz_x)
ad <- subset(ad,laz_x!=1)


# Exclude children with missing data (none)
table(is.na(ad$laz))

# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

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

# sort the data for perfect replication with jade on the V-fold cross-validation
ad <- ad[order(ad$block,ad$clusterid,ad$dataid,ad$childid),]

#---------------------------------------
# Select covariates with univariate
# associations with the outcome of
# P<0.2 based on a liklihood ratio test
#---------------------------------------


# drop due to so many missing values?
# asset_clock
# birthord

Ws <- subset(ad,select=c("fracode","month","aged","sex","momage","momedu","momheight","hfiacat","Nlt18","Ncomp","watmin","elec","floor","walls","roof","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile"))


#---------------------------------------
# Estimate adjusted mean differences
#---------------------------------------
#---------------------------------------
# H1: Each intervention arm vs. Control
#---------------------------------------
h1.contrasts <- list(
  c("Control","Water"),
  c("Control","Sanitation"),
  c("Control","Handwashing"),
  c("Control","WSH"),
  c("Control","Nutrition"),
  c("Control","Nutrition + WSH")
)

# unadjusted estimates (paired t-test)
diff.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ad$laz,tr=ad$tr,strat=ad$block))
rownames(diff.h1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

# adjusted estimates (tmle)
cwfit    <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,family="gaussian",contrast=c("Control","Water"), seed=12345)
csfit    <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,contrast=c("Control","Sanitation"), seed=12345)
chfit    <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,contrast=c("Control","Handwashing"), seed=12345)
cwshfit  <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,contrast=c("Control","WSH"), seed=12345)
cnfit    <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,contrast=c("Control","Nutrition"), seed=12345)
cwshnfit <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,contrast=c("Control","Nutrition + WSH"), seed=12345)

# pull out the estimates from the tmle objects and summarize them in a matrix
tmle.summary <- function(x) {
  res <- c(x$estimates$ATE$psi,x$estimates$ATE$CI[1],x$estimates$ATE$CI[2],x$estimates$ATE$pvalue)
  names(res) <- c("diff","ci.lb","ci.ub","p")
  return(res)
}

tmle.h1 <- list(cwfit,csfit,chfit,cwshfit,cnfit,cwshnfit)
diff.tmle.h1 <- t(sapply(tmle.h1,tmle.summary))
rownames(diff.tmle.h1) <- rownames(diff.h1)

# print results
round(diff.h1,4)
round(diff.tmle.h1,4)

#---------------------------------------
# H2: Combined WSH versus single interventions
# NOT CURRENTLY RUN BECAUSE NOT PRE-SPECIFIED
#---------------------------------------
# h2.contrasts <- list(
#   c("Water","WSH"),
#   c("Sanitation","WSH"),
#   c("Handwashing","WSH")
# )
# 
# # unadjusted estimates (paired t-test)
# diff.h2 <- t(sapply(h2.contrasts,ITT.unadj,Y=ad$laz,tr=ad$tr,strat=ad$block))
# rownames(diff.h2) <- c("WSH v Water","WSH v Sanitation","WSH v Handwashing")
# 
# # adjusted estimates (tmle)
# wshwfit    <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,family="gaussian",contrast=c("Water","WSH"),seed=12345)
# wshhfit    <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,family="gaussian",contrast=c("Handwashing","WSH"),seed=12345)
# wshsfit    <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,family="gaussian",contrast=c("Sanitation","WSH"),seed=12345)
# 
# tmle.h2 <- list(wshwfit,wshhfit,wshsfit)
# diff.tmle.h2 <- t(sapply(tmle.h2,tmle.summary))
# rownames(diff.tmle.h2) <- rownames(diff.h2)

#---------------------------------------
# H3: WSH+Nutrition vs. WSH or Nutrition alone
#---------------------------------------
h3.contrasts <- list(
  c("Nutrition","Nutrition + WSH"),
  c("WSH","Nutrition + WSH")
)

# unadjusted estimates (paired t-test)
diff.h3 <- t(sapply(h3.contrasts,ITT.unadj,Y=ad$laz,tr=ad$tr,strat=ad$block))
rownames(diff.h3) <- c("Nutrition + WSH v Nutrition", "Nutrition + WSH v WSH")

# adjusted estimates (tmle)
nwshnfit   <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,family="gaussian",contrast=c("Nutrition","Nutrition + WSH"),seed=12345)
wshwshnfit <- washb_tmle(Y=ad$laz,tr=ad$tr,W=Ws,id=ad$block,family="gaussian",contrast=c("WSH","Nutrition + WSH"),seed=12345)

tmle.h3 <- list(nwshnfit,wshwshnfit)
diff.tmle.h3 <- t(sapply(tmle.h3,tmle.summary))
rownames(diff.tmle.h3) <- rownames(diff.h3)

# print results
round(diff.h3,4)
round(diff.tmle.h3,4)


#---------------------------------------
# re-name into pre-specified objects
#---------------------------------------

laz_t1_h1_diff_adj <- diff.tmle.h1
# laz_t1_h2_diff_adj <- diff.tmle.h2
laz_t1_h3_diff_adj <- diff.tmle.h3


#---------------------------------------
# Print and save results
#---------------------------------------
laz_t1_h1_diff_adj
laz_t1_h3_diff_adj

# add 'b' suffix for comparison with jade
laz_t1_h1_diff_adj_b <- laz_t1_h1_diff_adj
laz_t1_h3_diff_adj_b <- laz_t1_h3_diff_adj
rm(laz_t1_h1_diff_adj,laz_t1_h3_diff_adj)

# save everything except the datasets themselves
# that way we have all of the block-specific estimates if needed for plotting or additional stats
rm(list=c("d","ad"))
save.image(file="~/dropbox/WBB-primary-analysis/results/raw/ben/bangladesh-laz-adj-t1-ben.RData")




	