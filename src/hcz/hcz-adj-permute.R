
#---------------------------------------
# hcz-adj-permute.R
#
# ben arnold (benarnold@berkeley.edu)
#
# adjusted permutation tests for
# differences in hcz at year 1 and year 2
# follow-ups
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls()); library(here)
library(plyr)
library(coin)
library(SuperLearner)

# source the base functions
# which includes the permutation test function used below
source(here("src/basefns/washb-base-functions.R"))


#---------------------------------------
# Load the analysis dataset,
# the baseline covariate dataset
#---------------------------------------

bd <- read.csv(here("data/washb-bangladesh-enrol-public.csv"))

d <- read.csv(here("data/washb-bangladesh-anthro-public.csv"))

# merge in the treatment assignments
tr    <- read.csv(here('data/washb-bangladesh-tr-public.csv'))

d <- left_join(d,tr,by=c("clusterid","block"))
bd <- left_join(bd,tr,by=c("clusterid","block"))

# merge the baseline dataset to the follow-up dataset
ad <- merge(bd,d,by=c("dataid","clusterid","block","tr"),all.x=F,all.y=T)
dim(d)
dim(ad)

#---------------------------------------
# load the anthropometry analysis data
#---------------------------------------
ad$block <- as.factor(ad$block)

# ensure that month is coded as a factor
ad$month <- factor(ad$month)

# subset the anthropometry to target children (excluding siblings)
dim(ad)
ad <- subset(ad,tchild=="Target child")
dim(ad)

# Drop children with extreme hcz values
table(ad$hcz_x)
ad <- subset(ad,hcz_x!=1)

# Exclude children with missing data (none)
table(is.na(ad$hcz))

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

#---------------------------------------
# YEAR 1
#---------------------------------------
# subset to the relevant measurement
table(ad$svy)
ad1 <- subset(ad,svy==1)
dim(ad1)

# sort the data for perfect replication with jade on the V-fold cross-validation
ad1 <- ad1[order(ad1$block,ad1$clusterid,ad1$dataid,ad1$childid),]

# excluding: birth order (missing for 241 children in year 1)
# excluding: asset_clock (missing for 1000s at enrollment)

# pre-specified covariates (but not treatment)
Ws <- subset(ad1,select=c("fracode","month","aged","sex","momage","momedu","momheight","hfiacat","Nlt18","Ncomp","watmin","elec","floor","walls","roof","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile"))

# restrict to complete cases
SLd <- data.frame(id=ad1$clusterid,block=ad1$block,tr=ad1$tr,Y=ad1$hcz,Ws)
dim(SLd)
SLd <- SLd[complete.cases(SLd),]
dim(SLd)
table(SLd$tr)

# pre-screen the covariates for those associated with the outcome (LR test P<0.2)
# see washb_prescreen() and design_matrix() in the base functions
Wscreen <- washb_prescreen(Y=SLd$Y,Ws=SLd[,5:ncol(SLd)],family="gaussian",pval=0.2)
Wselect <- subset(SLd,select=Wscreen)
Wselect <- design_matrix(Wselect)

# algorithmic fit of the outcome as a function of selected covariates
set.seed(12345)
SLfit1 <- SuperLearner(Y=SLd$Y,X=Wselect,id=SLd$block,
                       family="gaussian",
                       SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet")
                       )
SLfit1
SLd$pY <- as.vector(predict(SLfit1)$pred)
SLd$r <- SLd$Y-SLd$pY


# Hypothesis 1 permutation tests
permute.1.C.W <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Water"),seed=12345)
permute.1.C.S <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Sanitation"),seed=12345)
permute.1.C.H <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Handwashing"),seed=12345)
permute.1.C.WSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","WSH"),seed=12345)
permute.1.C.N   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Nutrition"),seed=12345)
permute.1.C.NWSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Nutrition + WSH"),seed=12345)
t1h1res <- list(permute.1.C.W,permute.1.C.S,permute.1.C.H,permute.1.C.WSH,permute.1.C.N,permute.1.C.NWSH)

# Hypothesis 3 permutation tests
permute.1.N.WSHN   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Nutrition","Nutrition + WSH"),seed=12345)
permute.1.WSH.NWSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("WSH","Nutrition + WSH"),seed=12345)
t1h3res <- list(permute.1.N.WSHN,permute.1.WSH.NWSH)

#---------------------------------------
# YEAR 2
#---------------------------------------
# subset to the relevant measurement
table(ad$svy)
ad2 <- subset(ad,svy==2)
dim(ad2)

# sort the data for perfect replication with jade on the V-fold cross-validation
ad2 <- ad2[order(ad2$block,ad2$clusterid,ad2$dataid,ad2$childid),]

# pre-specified covariates (but not treatment)
Ws <- subset(ad2,select=c("fracode","month","aged","sex","birthord","momage","momedu","momheight","hfiacat","Nlt18","Ncomp","watmin","elec","floor","walls","roof","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile"))

# restrict to complete cases
SLd <- data.frame(id=ad2$clusterid,block=ad2$block,tr=ad2$tr,Y=ad2$hcz,Ws)
dim(SLd)
SLd <- SLd[complete.cases(SLd),]
dim(SLd)
table(SLd$tr)

# pre-screen the covariates for those associated with the outcome (LR test P<0.2)
# see washb_prescreen() and design_matrix() in the base functions
Wscreen <- washb_prescreen(Y=SLd$Y,Ws=SLd[,5:ncol(SLd)],family="gaussian",pval=0.2)
Wselect <- subset(SLd,select=Wscreen)
Wselect <- design_matrix(Wselect)

# algorithmic fit of the outcome as a function of selected covariates
set.seed(12345)
SLfit2 <- SuperLearner(Y=SLd$Y,X=Wselect,id=SLd$block,
                       family="gaussian",
                       SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet")
)
SLfit2
SLd$pY <- as.vector(predict(SLfit2)$pred)
SLd$r <- SLd$Y-SLd$pY

# Hypothesis 1 permutation tests
permute.2.C.W <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Water"),seed=12345)
permute.2.C.S <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Sanitation"),seed=12345)
permute.2.C.H <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Handwashing"),seed=12345)
permute.2.C.WSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","WSH"),seed=12345)
permute.2.C.N   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Nutrition"),seed=12345)
permute.2.C.NWSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Control","Nutrition + WSH"),seed=12345)
t2h1res <- list(permute.2.C.W,permute.2.C.S,permute.2.C.H,permute.2.C.WSH,permute.2.C.N,permute.2.C.NWSH)

# Hypothesis 3 permutation tests
permute.2.N.WSHN   <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("Nutrition","Nutrition + WSH"),seed=12345)
permute.2.WSH.NWSH <- washb.permute(Y=SLd$r,tr=SLd$tr,block=SLd$block,contrast=c("WSH","Nutrition + WSH"),seed=12345)
t2h3res <- list(permute.2.N.WSHN,permute.2.WSH.NWSH)


#---------------------------------------
# put objects in the standard format
#---------------------------------------
hcz_t1_h1_pval_adj <- as.matrix(sapply(t1h1res,function(x) x$p.value),nrow=6)
rownames(hcz_t1_h1_pval_adj) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

hcz_t1_h3_pval_adj <- as.matrix(sapply(t1h3res,function(x) x$p.value),nrow=6)
rownames(hcz_t1_h3_pval_adj) <- c("Nutrition + WSH v Nutrition","Nutrition + WSH v WSH")

hcz_t2_h1_pval_adj <- as.matrix(sapply(t2h1res,function(x) x$p.value),nrow=6)
rownames(hcz_t2_h1_pval_adj) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

hcz_t2_h3_pval_adj <- as.matrix(sapply(t2h3res,function(x) x$p.value),nrow=6)
rownames(hcz_t2_h3_pval_adj) <- c("Nutrition + WSH v Nutrition","Nutrition + WSH v WSH")


#---------------------------------------
# print results
#---------------------------------------
hcz_t1_h1_pval_adj

hcz_t1_h3_pval_adj

hcz_t2_h1_pval_adj

hcz_t2_h3_pval_adj

#---------------------------------------
# add suffix for replication
#---------------------------------------
# hcz_t1_h1_pval_adj_b <- hcz_t1_h1_pval_adj
# hcz_t1_h3_pval_adj_b <- hcz_t1_h3_pval_adj
# hcz_t2_h1_pval_adj_b <- hcz_t2_h1_pval_adj
# hcz_t2_h3_pval_adj_b <- hcz_t2_h3_pval_adj
# rm(hcz_t1_h1_pval_adj,hcz_t1_h3_pval_adj,hcz_t2_h1_pval_adj,hcz_t2_h3_pval_adj)

#---------------------------------------
# save all of the results
# excluding the datasets
#---------------------------------------
rm(bd,d,ad,ad1,ad2,SLd,Ws,Wselect)
rm(SLfit1,SLfit2,t1h1res,t1h3res,t2h1res,t2h3res)
save.image(here("results/bangladesh-hcz-adj-permute.RData"))


