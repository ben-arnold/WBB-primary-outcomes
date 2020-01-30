
#---------------------------------------
# diar-unadj.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate unadjusted prevalence ratios
# and differences in prevalence for diarrhea
#
# using the Mantel-Hanzel estimator
#
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-diar-public.csv
#
# output files:
# bangladesh-diar-unadj.RData
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
d_tr    <- read.csv(here("data/washb-bangladesh-tr-public.csv"))
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
ad$tr <- factor(ad$tr,levels=c("Water","Sanitation","Handwashing","Nutrition","WSH","Nutrition + WSH","Control"))


#---------------------------------------
# cross-tabs of final observations
# in the analysis, by survey round
#---------------------------------------
table(ad$tr,ad$diar7d,ad$svy)


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
pr.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ad$diar7d,tr=ad$tr,strat=ad$block,binomial=T))
rd.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ad$diar7d,tr=ad$tr,strat=ad$block,binomial=T,measure="RD"))
rownames(pr.h1) <- rownames(rd.h1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

#---------------------------------------
# Mantel-Haenszel PR and RD estimates
# H2: Combined WSH versus single interventions
#---------------------------------------
h2.contrasts <- list(
  c("Water","WSH"),
  c("Sanitation","WSH"),
  c("Handwashing","WSH")
)
pr.h2 <- t(sapply(h2.contrasts,ITT.unadj,Y=ad$diar7d,tr=ad$tr,strat=ad$block,binomial=T))
rd.h2 <- t(sapply(h2.contrasts,ITT.unadj,Y=ad$diar7d,tr=ad$tr,strat=ad$block,binomial=T,measure="RD"))
rownames(pr.h2) <- rownames(rd.h2) <- c("WSH v Water","WSH v Sanitation","WSH v Handwashing")


#---------------------------------------
# Create final objects (pre-specified names)
#---------------------------------------
diar_h1_pr_unadj <- pr.h1[,c(1,3,4,6)]
diar_h1_rd_unadj <- rd.h1[,c(1,3,4,6)]
diar_h2_pr_unadj <- pr.h2[,c(1,3,4,6)]
diar_h2_rd_unadj <- rd.h2[,c(1,3,4,6)]

# exponentiate the PR and 95% estimates
diar_h1_pr_unadj[,c(1,2,3)] <- exp(diar_h1_pr_unadj[,c(1,2,3)])
diar_h2_pr_unadj[,c(1,2,3)] <- exp(diar_h2_pr_unadj[,c(1,2,3)])
colnames(diar_h1_pr_unadj)[1] <- colnames(diar_h2_pr_unadj)[1] <- "PR"

#---------------------------------------
# Print and save results
#---------------------------------------
round(diar_h1_pr_unadj,4)

round(diar_h1_rd_unadj,4)

round(diar_h2_pr_unadj,4)

round(diar_h2_rd_unadj,4)


# add a 'b' suffix for comparison with jade
# diar_h1_pr_unadj_b <- diar_h1_pr_unadj
# diar_h1_rd_unadj_b <- diar_h1_rd_unadj
# diar_h2_pr_unadj_b <- diar_h2_pr_unadj
# diar_h2_rd_unadj_b <- diar_h2_rd_unadj
# rm(diar_h1_pr_unadj, diar_h1_rd_unadj, diar_h2_pr_unadj, diar_h2_rd_unadj)

# save everything except the datasets themselves
# that way we have all of the block-specific estimates if needed for plotting or additional stats
rm(list=c("d","ad"))
save.image(file=here("results/bangladesh-diar-unadj.RData"))





	