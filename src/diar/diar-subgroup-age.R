
#---------------------------------------
# diar-subgroup-age.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate unadjusted differences 
# in prevalence for diarrhea, stratified
# by index children vs. others
#
# using the Mantel-Hanzel estimator
#
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-diar-public.csv
#
# output files:
# bangladesh-diar-subgroup-age.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls()); library(here)
library(metafor)

# source the base functions
source(here("src/basefns/washb-base-functions.R"))


#---------------------------------------
# Load the analysis dataset
#---------------------------------------
d <- read.csv(here("data/washb-bangladesh-diar-public.csv"))

# merge in the treatment assignments
d_tr    <- read.csv(here('data/washb-bangladesh-tr-public.csv'))
d <- left_join(d,d_tr,by=c("clusterid","block"))

#---------------------------------------
# Subset the Data to Follow-up data only
#---------------------------------------
table(d$svy)
ad <- subset(d,svy>0)

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

# re-order the tr factor for convenience
# careful: this order needs to correspond to labels in the results objects
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))


# stratify the data into index children and siblings/others
adi <- subset(ad,tchild=="Target child")
ads <- subset(ad,tchild=="Sibling")

#---------------------------------------
# cross-tabs of final observations
# in the analysis, by survey round
#---------------------------------------
iN  <- table(adi$tr,adi$diar7d)
iNy <- table(adi$tr,adi$diar7d,adi$svy)
iN
iNy

sN <- table(ads$tr,ads$diar7d)
sNy <- table(ads$tr,ads$diar7d,ads$svy)
sN
sNy

#---------------------------------------
# Calculate unadjusted prevalences
# and 95% CIs by arm
#---------------------------------------

# Calculate means and influence-curve based 95% CIs
# stratified by year of measurement
arms <- levels(ad$tr)
imu.y1 <- sapply(arms,tmle.mean.est,Y=adi$diar7d,tr=adi$tr,svy=adi$svy,id=adi$clusterid,s=1)
imu.y2 <- sapply(arms,tmle.mean.est,Y=adi$diar7d,tr=adi$tr,svy=adi$svy,id=adi$clusterid,s=2)

smu.y1 <- sapply(arms,tmle.mean.est,Y=ads$diar7d,tr=ads$tr,svy=ads$svy,id=ads$clusterid,s=1)
smu.y2 <- sapply(arms,tmle.mean.est,Y=ads$diar7d,tr=ads$tr,svy=ads$svy,id=ads$clusterid,s=2)

# and over the entire follow-up period
adi$svy2 <- 1
ads$svy2 <- 1
imu <- sapply(arms,tmle.mean.est,Y=adi$diar7d,tr=adi$tr,svy=adi$svy2,id=adi$clusterid,s=1)
smu <- sapply(arms,tmle.mean.est,Y=ads$diar7d,tr=ads$tr,svy=ads$svy2,id=ads$clusterid,s=1)

# quick checks of n/N vs. means obtained from tmle
round(cbind(iN[,2]/rowSums(iN),imu[1,],iN[,2]/rowSums(iN)-imu[1,]),4)
round(cbind(iNy[,2,1]/rowSums(iNy[,,1]),imu.y1[1,],iNy[,2,1]/rowSums(iNy[,,1])-imu.y1[1,]),4)
round(cbind(iNy[,2,2]/rowSums(iNy[,,2]),imu.y2[1,],iNy[,2,2]/rowSums(iNy[,,2])-imu.y2[1,]),4)

round(cbind(sN[,2]/rowSums(sN),smu[1,],sN[,2]/rowSums(sN)-smu[1,]),4)
round(cbind(sNy[,2,1]/rowSums(sNy[,,1]),smu.y1[1,],sNy[,2,1]/rowSums(sNy[,,1])-smu.y1[1,]),4)
round(cbind(sNy[,2,2]/rowSums(sNy[,,2]),smu.y2[1,],sNy[,2,2]/rowSums(sNy[,,2])-smu.y2[1,]),4)

#---------------------------------------
# Mantel-Haenszel PR and RD estimates
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

# Index children
ipr.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=adi$diar7d,tr=adi$tr,strat=adi$block,binomial=T))
ird.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=adi$diar7d,tr=adi$tr,strat=adi$block,binomial=T,measure="RD"))
rownames(ipr.h1) <- rownames(ird.h1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

  # print results (note exponentiated PR and CI)
  round(ird.h1,3)
  round(cbind(exp(ipr.h1[,c(1,3,4)]),ipr.h1[,c(5,6)]),3)

# non-index children
spr.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ads$diar7d,tr=ads$tr,strat=ads$block,binomial=T))
srd.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ads$diar7d,tr=ads$tr,strat=ads$block,binomial=T,measure="RD"))
rownames(spr.h1) <- rownames(srd.h1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

  # print results (note exponentiated PR and CI)
  round(srd.h1,3)
  round(cbind(exp(spr.h1[,c(1,3,4)]),spr.h1[,c(5,6)]),3)

#---------------------------------------
# Mantel-Haenszel PR and RD estimates
# H1: Each intervention arm vs. Control
# further stratified by year of measurement
#---------------------------------------
# Index children, years 1 and 2
ipr.h1.y1 <- t(sapply(h1.contrasts,ITT.unadj,Y=adi$diar7d[adi$svy==1],tr=adi$tr[adi$svy==1],strat=adi$block[adi$svy==1],binomial=T))
ird.h1.y1 <- t(sapply(h1.contrasts,ITT.unadj,Y=adi$diar7d[adi$svy==1],tr=adi$tr[adi$svy==1],strat=adi$block[adi$svy==1],binomial=T,measure="RD"))
rownames(ipr.h1.y1) <- rownames(ird.h1.y1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

  # print results (note exponentiated PR and CI)
  round(ird.h1.y1,3)
  round(cbind(exp(ipr.h1.y1[,c(1,3,4)]),ipr.h1.y1[,c(5,6)]),3)

ipr.h1.y2 <- t(sapply(h1.contrasts,ITT.unadj,Y=adi$diar7d[adi$svy==2],tr=adi$tr[adi$svy==2],strat=adi$block[adi$svy==2],binomial=T))
ird.h1.y2 <- t(sapply(h1.contrasts,ITT.unadj,Y=adi$diar7d[adi$svy==2],tr=adi$tr[adi$svy==2],strat=adi$block[adi$svy==2],binomial=T,measure="RD"))
rownames(ipr.h1.y2) <- rownames(ird.h1.y2) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

  # print results (note exponentiated PR and CI)
  round(ird.h1.y2,3)
  round(cbind(exp(ipr.h1.y2[,c(1,3,4)]),ipr.h1.y2[,c(5,6)]),3)

# non-index children, years 1 and 2
spr.h1.y1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ads$diar7d[ads$svy==1],tr=ads$tr[ads$svy==1],strat=ads$block[ads$svy==1],binomial=T))
srd.h1.y1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ads$diar7d[ads$svy==1],tr=ads$tr[ads$svy==1],strat=ads$block[ads$svy==1],binomial=T,measure="RD"))
rownames(spr.h1.y1) <- rownames(srd.h1.y1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

  # print results (note exponentiated PR and CI)
  round(srd.h1.y1,3)
  round(cbind(exp(spr.h1.y1[,c(1,3,4)]),spr.h1.y1[,c(5,6)]),3)

spr.h1.y2 <- t(sapply(h1.contrasts,ITT.unadj,Y=ads$diar7d[ads$svy==2],tr=ads$tr[ads$svy==2],strat=ads$block[ads$svy==2],binomial=T))
srd.h1.y2 <- t(sapply(h1.contrasts,ITT.unadj,Y=ads$diar7d[ads$svy==2],tr=ads$tr[ads$svy==2],strat=ads$block[ads$svy==2],binomial=T,measure="RD"))
rownames(spr.h1.y2) <- rownames(srd.h1.y2) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

  # print results (note exponentiated PR and CI)
  round(srd.h1.y2,3)
  round(cbind(exp(spr.h1.y2[,c(1,3,4)]),spr.h1.y2[,c(5,6)]),3)


#---------------------------------------
# Create final objects 
# combine N/n and prevalence estimates with
# RD estimates
#---------------------------------------

diar_h1_rd_age <- cbind(rowSums(iN),iN[,2],imu[1,],
                        rbind(rep(NA,3),ird.h1[,c(1,3,4)]),
                        rowSums(sN),sN[,2],smu[1,],
                        rbind(rep(NA,3),srd.h1[,c(1,3,4)])
                        )

diar_h1_rd_age_y1 <- cbind(rowSums(iNy[,,1]),iNy[,2,1],imu.y1[1,],
                        rbind(rep(NA,3),ird.h1.y1[,c(1,3,4)]),
                        rowSums(sNy[,,1]),sNy[,2,1],smu.y1[1,],
                        rbind(rep(NA,3),srd.h1.y1[,c(1,3,4)])
                        )

diar_h1_rd_age_y2 <- cbind(rowSums(iNy[,,2]),iNy[,2,2],imu.y2[1,],
                           rbind(rep(NA,3),ird.h1.y2[,c(1,3,4)]),
                           rowSums(sNy[,,2]),sNy[,2,2],smu.y2[1,],
                           rbind(rep(NA,3),srd.h1.y2[,c(1,3,4)])
)

colnames(diar_h1_rd_age) <- colnames(diar_h1_rd_age_y1) <- colnames(diar_h1_rd_age_y2) <- c(paste(rep(c("index","other"),c(6,6)),c("N","n","prev","RD","RDlb","RDub"),sep="."))


#---------------------------------------
# Create final objects 
# combine N/n and prevalence estimates with
# PR estimates
#---------------------------------------

diar_h1_pr_age <- cbind(rowSums(iN),iN[,2],imu[1,],
                        rbind(rep(NA,3),exp(ipr.h1[,c(1,3,4)])),
                        rowSums(sN),sN[,2],smu[1,],
                        rbind(rep(NA,3),exp(spr.h1[,c(1,3,4)]))
)

diar_h1_pr_age_y1 <- cbind(rowSums(iNy[,,1]),iNy[,2,1],imu.y1[1,],
                           rbind(rep(NA,3),exp(ipr.h1.y1[,c(1,3,4)])),
                           rowSums(sNy[,,1]),sNy[,2,1],smu.y1[1,],
                           rbind(rep(NA,3),exp(spr.h1.y1[,c(1,3,4)]))
)

diar_h1_pr_age_y2 <- cbind(rowSums(iNy[,,2]),iNy[,2,2],imu.y2[1,],
                           rbind(rep(NA,3),exp(ipr.h1.y2[,c(1,3,4)])),
                           rowSums(sNy[,,2]),sNy[,2,2],smu.y2[1,],
                           rbind(rep(NA,3),exp(spr.h1.y2[,c(1,3,4)]))
)

colnames(diar_h1_pr_age) <- colnames(diar_h1_pr_age_y1) <- colnames(diar_h1_pr_age_y2) <- c(paste(rep(c("index","other"),c(6,6)),c("N","n","prev","PR","PRlb","PRub"),sep="."))



#---------------------------------------
# Print and save results
#---------------------------------------
round(diar_h1_rd_age[,1:6],3)
round(diar_h1_rd_age[,7:12],3)

round(diar_h1_rd_age_y1[,1:6],3)
round(diar_h1_rd_age_y1[,7:12],3)

round(diar_h1_rd_age_y2[,1:6],3)
round(diar_h1_rd_age_y2[,7:12],3)

round(diar_h1_pr_age[,1:6],3)
round(diar_h1_pr_age[,7:12],3)

round(diar_h1_pr_age_y1[,1:6],3)
round(diar_h1_pr_age_y1[,7:12],3)

round(diar_h1_pr_age_y2[,1:6],3)
round(diar_h1_pr_age_y2[,7:12],3)

# save everything except the datasets themselves
# that way we have all of the block-specific estimates if needed for plotting or additional stats
rm(list=c("d","ad","adi","ads"))
save.image(file=here("results/bangladesh-diar-subgroup-age.RData"))





	