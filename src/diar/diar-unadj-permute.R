
#---------------------------------------
# diar-unadj-permute.R
#
# ben arnold (benarnold@berkeley.edu)
#
# unadjusted permutation tests for
# differences in diarrhea prevalence
# for the year 1 and year 2 
# visits combined
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
# load all packages
source(here::here("src/0-config.R"))

# source the base functions
# which includes the permutation test function used below
source(here("src/basefns/washb-base-functions.R"))


#---------------------------------------
# load the diarrhea analysis data
#---------------------------------------
d <- read.csv(here("data/washb-bangladesh-diar-public.csv"))

# merge in the treatment assignments
d_tr    <- read.csv(here('data/washb-bangladesh-tr-public.csv'))
d <- left_join(d,d_tr,by=c("clusterid","block"))

d$block <- as.factor(d$block)

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
# Permutation tests
#---------------------------------------

# Hypothesis 1 permutation tests
set.seed(242524)
permute.C.W <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Control","Water"))
permute.C.S <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Control","Sanitation"))
permute.C.H <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Control","Handwashing"))
permute.C.WSH <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Control","WSH"))
permute.C.N   <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Control","Nutrition"))
permute.C.NWSH <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Control","Nutrition + WSH"))

# Hypothesis 2 permutation tests
set.seed(35234)
permute.W.WSH   <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Water","WSH"))
permute.S.WSH   <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Sanitation","WSH"))
permute.H.WSH   <- washb.permute(Y=ad$diar7d,tr=ad$tr,block=ad$block,c("Handwashing","WSH"))

#---------------------------------------
# put objects in the standard format
#---------------------------------------
h1res <- list(permute.C.W,permute.C.S,permute.C.H,permute.C.WSH,permute.C.N,permute.C.NWSH)
diar_h1_pval_unadj <- as.matrix(sapply(h1res,function(x) x$p.value),nrow=6)
rownames(diar_h1_pval_unadj) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

h2res <- list(permute.W.WSH,permute.S.WSH,permute.H.WSH)
diar_h2_pval_unadj <- as.matrix(sapply(h2res,function(x) x$p.value),nrow=3)
rownames(diar_h2_pval_unadj) <- c("WSH v W","WSH v S","WSH v H")


#---------------------------------------
# print results
#---------------------------------------
diar_h1_pval_unadj

diar_h2_pval_unadj


#---------------------------------------
# add suffix for replication
#---------------------------------------
# diar_h1_pval_unadj_b <- diar_h1_pval_unadj
# diar_h2_pval_unadj_b <- diar_h2_pval_unadj
# rm(diar_h1_pval_unadj,diar_h2_pval_unadj)

#---------------------------------------
# save all of the results
# excluding the datasets
#---------------------------------------
rm(d,ad)
save.image(here("results/bangladesh-diar-unadj-permute.RData"))


