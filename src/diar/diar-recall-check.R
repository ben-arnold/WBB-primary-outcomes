
#---------------------------------------
# diar-recall-check.R
#
# ben arnold (benarnold@berkeley.edu)
#
# Check for differential recall errors
# in diarrhea by comparing current cases
# (2d recall) with terminated cases
# (those who answered "Yes" to 7d recall
# but "No" to 2d recall) using the C/T ratio.
#
# This is the approach used in Boerma et al 1991
# and recommended by Arnold et al. 2013 (appendix 4)
#
# Arnold et al. 2013. Optimal Recall Period for Caregiver-Reported Illness 
# in Risk Factor and Intervention Studies: A Multicountry Study.
# American Journal of Epidemiology 177 (4): 361–70.
#
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-diar-public.csv
#
# output files:
# bangladesh-diar-recall-check.RData
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
ad$tr <- factor(ad$tr,levels=c("Water","Sanitation","Handwashing","Nutrition","WSH","Nutrition + WSH","Control"))

# create a variable for terminated diarrhea cases between 7d and 2d recall periods
ad$diart <- ifelse(ad$diar2d==0 & ad$diar7d==1,1,0)
table(ad$diar2d,ad$diart)  # only 146 terminated cases

#---------------------------------------
# Current:Terminated ratio difference function
#---------------------------------------
ctratio <- function(dcurr,dterm,tr) {
  # dcurr: indicator of a current diarrhea case
  # dterm: indicator of a terminated diarrhea case
  # tr   : assigned treatment group
  curr <- tapply(dcurr,tr,sum)
  term <- tapply(dterm,tr,sum)
  ctratio <- curr/term
  return(ctratio)
}

ctratio(ad$diar2d,ad$diart,ad$tr)



#---------------------------------------
# cross-tabs of final observations
# in the analysis, by survey round
#---------------------------------------
table(ad$tr,ad$diar2d)
table(ad$tr,ad$diar7d)
table(ad$tr,ad$diart)


#---------------------------------------
# bootstrap the CT ratio by re-sampling
# randomization blocks with replacement
#---------------------------------------

set.seed(1349175)
nreps <- 1000
bsamp <- matrix(sample(unique(ad$block),size=length(unique(ad$block))*nreps,replace=TRUE),ncol=nreps)
ctratios <- matrix(rep(NA,nreps*7),ncol=7)
for(i in 1:nreps) {
  bd <- merge(ad,data.frame(block=bsamp[,i]),by="block",all.x=FALSE)
  ctratios[i,] <- ctratio(bd$diar2d,bd$diart,bd$tr)
}

# compute differences between each arm and control, then the mean and percentile 95% CIs
muctratio <- apply(ctratios,2,mean)[c(7,1:6)]
ctdiff <- ctratios[,1:6]-ctratios[,7]
ctmeans <- apply(ctdiff,2,mean)
ct95ci  <- apply(ctdiff,2,function(x) quantile(x,probs=c(0.025,0.975)))

res <- cbind(muctratio,c(NA,ctmeans),t(cbind(c(NA,NA),ct95ci)))
rownames(res) <- levels(ad$tr)[c(7,1:6)]
colnames(res) <- c("CTratio","CTratio diff","Lower 95% CI","Upper 95% CI")


# print results
res


# save everything except the datasets themselves
# that way we have all of the block-specific estimates if needed for plotting or additional stats
rm(list=c("d","ad"))
save.image(file=here("results/bangladesh-diar-recall-check.RData"))





	