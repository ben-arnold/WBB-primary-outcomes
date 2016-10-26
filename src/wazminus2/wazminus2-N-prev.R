
#---------------------------------------
# wazminus2-N-prev.R
#
# ben arnold (benarnold@berkeley.edu)
#
# calculate the number of wazminus2
# cases and the prevalence by
# survey round and treatment arm
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-anthro.csv
#
# output files:
#	bangladesh-wazminus2-N-prev-ben.RData
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
table(d$waz_x)
ad <- subset(d,waz_x!=1)

#---------------------------------------
# Exclude children with missing data
# and subset the data to target children
#---------------------------------------
table(is.na(ad$wazminus2))
ad <- subset(ad,!is.na(ad$wazminus2))

table(ad$tchild)
ad <- subset(ad,tchild=="Target child")
dim(ad)

# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

#---------------------------------------
# Calculate N and wazminus2 cases
# by survey round and by treatment arm
#---------------------------------------

Nchild    <- tapply(ad$wazminus2,list(ad$svy,ad$tr),function(x) length(x))
wazminus2 <- tapply(ad$wazminus2,list(ad$svy,ad$tr),function(x) sum(x))

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
mu1 <- sapply(arms,tmle.mean.est,Y=ad$wazminus2,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=1)
mu2 <- sapply(arms,tmle.mean.est,Y=ad$wazminus2,tr=ad$tr,svy=ad$svy,id=ad$clusterid,s=2)

#---------------------------------------
# Create final objects
#---------------------------------------

# condense results into pre-specfied objects
underwt_t1_n <- cbind(Nchild[1,],wazminus2[1,])
underwt_t2_n <- cbind(Nchild[2,],wazminus2[2,])
colnames(underwt_t1_n) <-  colnames(underwt_t2_n) <- c("N","n")

underwt_t1_prev <- t(mu1)
underwt_t2_prev <- t(mu2)

#---------------------------------------
# Internal consistency check:
# the n/N and tmle prevalence estimates
# should be the same
#---------------------------------------
round(cbind(underwt_t1_n[,2]/underwt_t1_n[,1],underwt_t1_prev[,1]),4)
round(cbind(underwt_t2_n[,2]/underwt_t2_n[,1],underwt_t2_prev[,1]),4)

#---------------------------------------
# Print and save results
#---------------------------------------

underwt_t1_n
underwt_t2_n
round(underwt_t1_prev,4)
round(underwt_t2_prev,4)

# add 'b' suffix for comparison w/ jade
underwt_t1_n_b <- underwt_t1_n
underwt_t2_n_b <- underwt_t2_n
underwt_t1_prev_b <- underwt_t1_prev
underwt_t2_prev_b <- underwt_t2_prev


save(underwt_t1_n_b,underwt_t2_n_b,underwt_t1_prev_b,underwt_t2_prev_b,file="~/dropbox/wbb-primary-analysis/results/raw/ben/bangladesh-wazminus2-N-prev-ben.RData")





	