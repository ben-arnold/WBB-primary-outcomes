
#---------------------------------------
# whzminus2-N-prev.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate the number of whzminus2
# cases and the prevalence by
# survey round and treatment arm
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-anthro.csv
#
# output files:
#	bangladesh-whzminus2-N-prev-ben.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(tmle)

#---------------------------------------
# Load the analysis dataset
#---------------------------------------
d <- read.csv("~/dropbox/wbb-primary-analysis/data/final/ben/washb-bangladesh-anthro.csv")

#---------------------------------------
# Drop children with extreme LAZ values
#---------------------------------------
table(d$whz_x)
ad <- subset(d,whz_x!=1)

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
# Calculate N and whzminus2 cases
# by survey round and by treatment arm
#---------------------------------------

Nchild    <- tapply(ad$whzminus2,list(ad$svy,ad$tr),function(x) length(x))
whzminus2 <- tapply(ad$whzminus2,list(ad$svy,ad$tr),function(x) sum(x))

#---------------------------------------
# Calculate unadjusted prevalences
# and 95% CIs by arm
#---------------------------------------
# wrapper function to call for each treatment and survey round
# has to be applied to the ad data frame with diar7d and tr variables
tmle.mean.est <- function(Y,tr,svy,id,group="Control",s=0) {
	# Y : outcome variable
	# tr: treatment indicator variable
  # svy  : measurment round variable
  # id: cluster ID variable
  # group : string. treatment factor level to compute mean
  # s     : survey round to compute mean. 0, 1, or 2 
  tmledat <- data.frame(id=id[tr==group & svy==s],
                        svy=svy[tr==group & svy==s],
                        Y=Y[tr==group & svy==s],
                        tr=tr[tr==group & svy==s])
  mu.fit <- tmle(Y=tmledat$Y,A=NULL,W=as.matrix(rep(1,nrow(tmledat)),nrow=length(Y)),id=tmledat$id,Q.SL.library="SL.mean")
	print(mu.fit)
	mu <- mu.fit$estimates$EY1$psi
	se <- sqrt(mu.fit$estimates$EY1$var.psi)
	ci <- mu.fit$estimates$EY1$CI
	res <- c(mu,ci[1],ci[2])
	names(res) <- c("mean","ci.lb","ci.ub")
	return(res)
}


# Calculate means and influence-curve based 95% CIs by survey round
set.seed(12345)
arms <- levels(ad$tr)
mu1 <- sapply(arms,tmle.mean.est,Y=ad$whzminus2,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=1)
mu2 <- sapply(arms,tmle.mean.est,Y=ad$whzminus2,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=2)

#---------------------------------------
# Create final objects
#---------------------------------------

# condense results into pre-specfied objects
wast_t1_n <- cbind(Nchild[1,],whzminus2[1,])
wast_t2_n <- cbind(Nchild[2,],whzminus2[2,])
colnames(wast_t1_n) <-  colnames(wast_t2_n) <- c("N","n")

wast_t1_prev <- t(mu1)
wast_t2_prev <- t(mu2)

#---------------------------------------
# Internal consistency check:
# the n/N and tmle prevalence estimates
# should be the same
#---------------------------------------
round(cbind(wast_t1_n[,2]/wast_t1_n[,1],wast_t1_prev[,1]),4)
round(cbind(wast_t2_n[,2]/wast_t2_n[,1],wast_t2_prev[,1]),4)

#---------------------------------------
# Print and save results
#---------------------------------------

wast_t1_n
wast_t2_n
round(wast_t1_prev,4)
round(wast_t2_prev,4)

# add 'b' suffix for comparison w/ jade
wast_t1_n_b <- wast_t1_n
wast_t2_n_b <- wast_t2_n
wast_t1_prev_b <- wast_t1_prev
wast_t2_prev_b <- wast_t2_prev


save(wast_t1_n_b,wast_t2_n_b,wast_t1_prev_b,wast_t2_prev_b,file="~/dropbox/wbb-primary-analysis/results/raw/ben/bangladesh-whzminus2-N-prev-ben.RData")





	