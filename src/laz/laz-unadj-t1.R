
#---------------------------------------
# laz-diff-unadj-t1.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate unadjusted differences
# between treatment arms for H1 and H3
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-anthro.csv
#
# output files:
#	bangladesh-laz-unadj-t1-ben.RData
#
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())

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
ad <- subset(d,svy==1)
dim(ad)

#---------------------------------------
# Drop children with extreme LAZ values
#---------------------------------------
table(ad$laz_x)
ad <- subset(ad,laz_x!=1)

#---------------------------------------
# Exclude children with missing data
# and subset the data to target children
#---------------------------------------
table(is.na(ad$laz))
ad <- subset(ad,!is.na(ad$laz))

table(ad$tchild)
ad <- subset(ad,tchild=="Target child")
dim(ad)

# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))


#---------------------------------------
# Calculate N, mean and SD of LAZ
#---------------------------------------
laz.n <- tapply(ad$laz,ad$tr,function(x) length(x))
laz.mu <- tapply(ad$laz,ad$tr,function(x) mean(x))
laz.sd <- tapply(ad$laz,ad$tr,function(x) sd(x))

# condense results into pre-specfied objects
laz_t1_n <- cbind(laz.n, laz.mu, laz.sd)

# final labeling
colnames(laz_t1_n) <- c("N","Mean","SD")

# print
laz_t1_n

# add 'b' suffix for comparison with jade
laz_t1_n_b <- laz_t1_n
rm(laz_t1_n)

#---------------------------------------
# Estimate paired t-tests for differences
# in means at the randomization block level
#---------------------------------------

#---------------------------------------
# paired T-test for differences in LAZ
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
diff.h1 <- t(sapply(h1.contrasts,ITT.unadj,Y=ad$laz,tr=ad$tr,strat=ad$block))
rownames(diff.h1) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

#---------------------------------------
# paired T-test for differences in LAZ
# H2: Combined WSH versus single interventions
#---------------------------------------
h2.contrasts <- list(
  c("Water","WSH"),
  c("Sanitation","WSH"),
  c("Handwashing","WSH")
)
diff.h2 <- t(sapply(h2.contrasts,ITT.unadj,Y=ad$laz,tr=ad$tr,strat=ad$block))
rownames(diff.h2) <- c("WSH v Water","WSH v Sanitation","WSH v Handwashing")

#---------------------------------------
# paired T-test for differences in LAZ
# H3: WSH+Nutrition vs. WSH or Nutrition alone
#---------------------------------------
h3.contrasts <- list(
  c("Nutrition","Nutrition + WSH"),
  c("WSH","Nutrition + WSH")
)
diff.h3 <- t(sapply(h3.contrasts,ITT.unadj,Y=ad$laz,tr=ad$tr,strat=ad$block))
rownames(diff.h3) <- c("Nutrition + WSH v Nutrition","Nutrition + WSH v WSH")


#---------------------------------------
# re-name into pre-specified objects
#---------------------------------------

laz_t1_h1_diff_unadj <- diff.h1
laz_t1_h2_diff_unadj <- diff.h2
laz_t1_h3_diff_unadj <- diff.h3


#---------------------------------------
# Print and save results
#---------------------------------------
round(laz_t1_h1_diff_unadj,4)

round(laz_t1_h2_diff_unadj,4)

round(laz_t1_h3_diff_unadj,4)

# add 'b' suffix for comparison with jade
laz_t1_h1_diff_unadj_b <- laz_t1_h1_diff_unadj
laz_t1_h2_diff_unadj_b <- laz_t1_h2_diff_unadj
laz_t1_h3_diff_unadj_b <- laz_t1_h3_diff_unadj
rm(laz_t1_h1_diff_unadj,laz_t1_h2_diff_unadj,laz_t1_h3_diff_unadj)

# save everything except the datasets themselves
# that way we have all of the block-specific estimates if needed for plotting or additional stats
rm(list=c("d","ad"))
save.image(file="~/dropbox/wbb-primary-analysis/results/raw/ben/bangladesh-laz-unadj-t1-ben.RData")




	