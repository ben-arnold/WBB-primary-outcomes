
#---------------------------------------
# whzminus2-unadj-t2.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate unadjusted prevalence ratios
# and differences in prevalence for whzminus2
#
# using the Mantel-Hanzel estimator
#
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-anthro.csv
#
# output files:
# bangladesh-whzminus2-unadj-t2.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(metafor)

# source the base functions
source("~/WBBpa/src/basefns/washb-base-functions.R")


#---------------------------------------
# Load the analysis dataset
#---------------------------------------
d <- read.csv("~/dropbox/wbb-primary-analysis/data/final/ben/washb-bangladesh-anthro.csv")


#---------------------------------------
# subset to the relevant measurement
# Year 1 or Year 2
#---------------------------------------
table(d$svy)
ad <- subset(d,svy==2)
dim(ad)

#---------------------------------------
# Drop children with extreme LAZ values
#---------------------------------------
table(ad$whz_x)
ad <- subset(ad,whz_x!=1)

#---------------------------------------
# Exclude children with missing data
# and subset the data to target children
#---------------------------------------
table(is.na(ad$whzminus2))
ad <- subset(ad,!is.na(ad$whzminus2))

table(ad$tchild)
ad <- subset(ad,tchild=="Target child")
dim(ad)

# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

#---------------------------------------
# cross-tabs of final observations
# in the analysis, by group
#---------------------------------------
table(ad$tr,ad$whzminus2)


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
pr.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ad$whzminus2,tr=ad$tr,strat=ad$block,binomial=T))
rd.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ad$whzminus2,tr=ad$tr,strat=ad$block,binomial=T,measure="RD"))
rownames(pr.h1) <- rownames(rd.h1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

#---------------------------------------
# Mantel-Haenszel PR and RD estimates
# H3: WSH+Nutrition vs. WSH or Nutrition alone
#---------------------------------------
h3.contrasts <- list(
  c("Nutrition","Nutrition + WSH"),
  c("WSH","Nutrition + WSH")
)
pr.h3 <- t(sapply(h3.contrasts,ITT.unadj,Y=ad$whzminus2,tr=ad$tr,strat=ad$block,binomial=T))
rd.h3 <- t(sapply(h3.contrasts,ITT.unadj,Y=ad$whzminus2,tr=ad$tr,strat=ad$block,binomial=T,measure="RD"))
rownames(pr.h3) <- rownames(rd.h3) <- c("Nutrition + WSH v Nutrition","Nutrition + WSH v WSH")


#---------------------------------------
# Create final objects (pre-specified names)
#---------------------------------------
wast_t2_h1_pr_unadj <- pr.h1[,c(1,3,4,6)]
wast_t2_h1_rd_unadj <- rd.h1[,c(1,3,4,6)]
wast_t2_h3_pr_unadj <- pr.h3[,c(1,3,4,6)]
wast_t2_h3_rd_unadj <- rd.h3[,c(1,3,4,6)]

# exponentiate the PR and 95% estimates
wast_t2_h1_pr_unadj[,c(1,2,3)] <- exp(wast_t2_h1_pr_unadj[,c(1,2,3)])
wast_t2_h3_pr_unadj[,c(1,2,3)] <- exp(wast_t2_h3_pr_unadj[,c(1,2,3)])
colnames(wast_t2_h1_pr_unadj)[1] <- colnames(wast_t2_h3_pr_unadj)[1] <- "PR"

#---------------------------------------
# Print and save results
#---------------------------------------
round(wast_t2_h1_pr_unadj,4)

round(wast_t2_h1_rd_unadj,4)

round(wast_t2_h3_pr_unadj,4)

round(wast_t2_h3_rd_unadj,4)


# add a 'b' suffix for comparison with jade
wast_t2_h1_pr_unadj_b <- wast_t2_h1_pr_unadj
wast_t2_h1_rd_unadj_b <- wast_t2_h1_rd_unadj
wast_t2_h3_pr_unadj_b <- wast_t2_h3_pr_unadj
wast_t2_h3_rd_unadj_b <- wast_t2_h3_rd_unadj
rm(wast_t2_h1_pr_unadj, wast_t2_h1_rd_unadj, wast_t2_h3_pr_unadj, wast_t2_h3_rd_unadj)

# save everything except the datasets themselves
# that way we have all of the block-specific estimates if needed for plotting or additional stats
rm(list=c("d","ad"))
save.image(file="~/dropbox/wbb-primary-analysis/results/raw/ben/bangladesh-whzminus2-unadj-t2-ben.RData")





	
