
#---------------------------------------
# diar-adj.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate adjusted differences
# between treatment arms for H1 and H2
#
# diarrhea
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-diar-public.csv
#
# output files:
#	bangladesh-diar-adj-ben.RData
#
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls()); library(here)
library(tmle)
library(SuperLearner)

# source the base functions
source(here("src/basefns/washb-base-functions.R"))


#---------------------------------------
# Load the analysis dataset,
# the baseline covariate dataset
#---------------------------------------

bd <- read.csv(here("data/washb-bangladesh-enrol-public.csv"))

  # drop svydate and month because they are superceded in the child level diarrhea data
  bd$svydate <- NULL
  bd$month <- NULL

d <- read.csv(here("data/washb-bangladesh-diar-public.csv"))

# merge in the treatment assignments
tr <- read.csv(here("data/washb-bangladesh-tr-public.csv"))
bd_tr <- left_join(bd,tr,by=c("clusterid","block"))
d_tr <- left_join(d,tr,by=c("clusterid","block"))

# merge the baseline dataset to the follow-up dataset
ad <- merge(bd_tr,d_tr,by=c("dataid","clusterid","block","tr"),all.x=F,all.y=T)
dim(d_tr)
dim(ad)

#---------------------------------------
# subset to the relevant measurement
# Year 1 or Year 2
#---------------------------------------
table(ad$svy)
ad <- subset(ad,svy==1|svy==2)
dim(ad)

# subset the diarrhea to children <36 mos at enrollment
# (exlude new births that are not target children)
dim(ad)
table(ad$sibnewbirth)
ad <- subset(ad,sibnewbirth==0)
dim(ad)
table(ad$gt36mos)
ad <- subset(ad,gt36mos==0)


# Exclude children with missing data
table(ad$tchild,is.na(ad$diar7d),ad$svy)
ad <- subset(ad,!is.na(ad$diar7d))

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

# Excluded because only measured in target children:
# birthord

# drop due to so many missing values?
# asset_clock
# momheight?

Ws <- subset(ad,select=c("fracode","month","agedays","sex","momage","momedu","momheight","hfiacat","Nlt18","Ncomp","watmin","elec","floor","walls","roof","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile"))


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

# unadjusted estimates (mantel-haenszel)
diff.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ad$diar7d,tr=ad$tr,strat=ad$block,binomial=TRUE,measure="RD"))
rownames(diff.h1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

# adjusted estimates (tmle)
cwfit    <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Control","Water"),seed=12345)
csfit    <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Control","Sanitation"),seed=12345)
chfit    <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Control","Handwashing"),seed=12345)
cwshfit  <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Control","WSH"),seed=12345)
cnfit    <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Control","Nutrition"),seed=12345)
cwshnfit <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Control","Nutrition + WSH"),seed=12345)

# pull out the estimates from the tmle objects and summarize them in a matrix
tmle.summary.rd <- function(x) {
  res <- c(x$estimates$ATE$psi,x$estimates$ATE$CI[1],x$estimates$ATE$CI[2],x$estimates$ATE$pvalue)
  names(res) <- c("rd","ci.lb","ci.ub","p")
  return(res)
}
tmle.summary.pr <- function(x) {
  res <- c(x$estimates$RR$psi,x$estimates$RR$CI[1],x$estimates$RR$CI[2],x$estimates$RR$pvalue)
  names(res) <- c("pr","ci.lb","ci.ub","p")
  return(res)
}

tmle.h1 <- list(cwfit,csfit,chfit,cwshfit,cnfit,cwshnfit)
rd.tmle.h1 <- t(sapply(tmle.h1,tmle.summary.rd))
pr.tmle.h1 <- t(sapply(tmle.h1,tmle.summary.pr))
rownames(rd.tmle.h1) <- rownames(pr.tmle.h1) <- rownames(diff.h1)

# print results
round(diff.h1,4)
round(rd.tmle.h1,4)
round(pr.tmle.h1,4)


#---------------------------------------
# H2: Combined WSH versus single interventions
#---------------------------------------
h2.contrasts <- list(
  c("Water","WSH"),
  c("Sanitation","WSH"),
  c("Handwashing","WSH")
)

# unadjusted estimates (paired t-test)
diff.h2 <- t(sapply(h2.contrasts,ITT.unadj,Y=ad$diar7d,tr=ad$tr,strat=ad$block,binomial=TRUE,measure="RD"))
rownames(diff.h2) <- c("WSH v Water","WSH v Sanitation","WSH v Handwashing")

# adjusted estimates (tmle)
wshwfit    <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Water","WSH"),seed=12345)
wshhfit    <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Handwashing","WSH"),seed=12345)
wshsfit    <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Sanitation","WSH"),seed=12345)

tmle.h2 <- list(wshwfit,wshsfit,wshhfit)
rd.tmle.h2 <- t(sapply(tmle.h2,tmle.summary.rd))
pr.tmle.h2 <- t(sapply(tmle.h2,tmle.summary.pr))
rownames(rd.tmle.h2) <- rownames(pr.tmle.h2) <- rownames(diff.h2)

# print results
round(diff.h2,4)
round(rd.tmle.h2,4)
round(pr.tmle.h2,4)

#---------------------------------------
# H3: WSH+Nutrition vs. WSH or Nutrition alone
# CURRENTLY COMMENTED OUT TO SAVE COMPUTING
# TIME B/C NOT PRE-SPECIFIED for diarrhea
#---------------------------------------
# h3.contrasts <- list(
#   c("WSH","Nutrition + WSH"),
#   c("Nutrition","Nutrition + WSH")
# )
# 
# # unadjusted estimates (paired t-test)
# diff.h3 <- t(sapply(h3.contrasts,ITT.unadj,Y=ad$diar7d,tr=ad$tr,strat=ad$block))
# rownames(diff.h3) <- c("Nutrition + WSH v Nutrition","Nutrition + WSH v WSH")
# 
# # adjusted estimates (tmle)
# wshwshnfit <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("WSH","Nutrition + WSH"),seed=12345)
# nwshnfit   <- washb_tmle(Y=ad$diar7d,tr=ad$tr,W=Ws,id=ad$block,family="binomial",contrast=c("Nutrition","Nutrition + WSH"),seed=12345)
# 
# tmle.h3 <- list(nwshnfit,wshwshnfit)
# rd.tmle.h3 <- t(sapply(tmle.h3,tmle.summary.rd))
# rownames(rd.tmle.h3) <- rownames(diff.h3)


#---------------------------------------
# re-name into pre-specified objects
#---------------------------------------

diar_h1_rd_adj <- rd.tmle.h1
diar_h2_rd_adj <- rd.tmle.h2
# diar_h3_rd_adj <- rd.tmle.h3

diar_h1_pr_adj <- pr.tmle.h1
diar_h2_pr_adj <- pr.tmle.h2
# diar_h3_pr_adj <- pr.tmle.h3

#---------------------------------------
# Print and save results
#---------------------------------------
diar_h1_rd_adj

diar_h2_rd_adj

diar_h1_pr_adj

diar_h2_pr_adj



# add 'b' suffix for comparison with jade
# diar_h1_rd_adj_b <- diar_h1_rd_adj
# diar_h2_rd_adj_b <- diar_h2_rd_adj
# diar_h1_pr_adj_b <- diar_h1_pr_adj
# diar_h2_pr_adj_b <- diar_h2_pr_adj
# rm(diar_h1_rd_adj,diar_h2_rd_adj,diar_h1_pr_adj,diar_h2_pr_adj)

# save everything except the datasets themselves
# that way we have all of the block-specific estimates if needed for plotting or additional stats
rm(list=c("d","ad","d_tr"))
save.image(file="results/bangladesh-diar-adj.RData")




	