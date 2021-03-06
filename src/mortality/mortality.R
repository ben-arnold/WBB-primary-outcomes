
#---------------------------------------
# mortality.R
#
# ben arnold (bfarnold@gmail.com)
#
# calculate unadjusted comparisons
# between treatment arms in cumulative
# all cause mortality during the trial
#---------------------------------------

#---------------------------------------
# input files:
# washb-bangladesh-track-compound.csv
# washb-bangladesh-enrol-tr.csv
#	washb-bangladesh-anthro.csv
#
# output files:
#	bangladesh-mortality.RData
#
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
library(here)
here()
library(tmle)
library(metafor)

# source the base functions
source(here::here("src/basefns/","washb-base-functions.R"))


#---------------------------------------
# Load the analysis dataset,
# the compound tracking dataset
#---------------------------------------

trd <- read.csv(here::here("data","washb-bangladesh-tr-public.csv"),colClasses=c("clusterid"="character"))

td <- read.csv(here::here("data","washb-bangladesh-track-compound-public.csv"),colClasses=c("dataid"="character"))

d <- read.csv(here::here("data","washb-bangladesh-anthro-public.csv"),colClasses=c("dataid"="character"))

# merge treatment assigments to the tracking dataset
ad <- merge(td,trd,by=c("block","clusterid"),all.x=T,all.y=T)
dim(td)
dim(ad)

# merge the tracking dataset to the follow-up dataset
ad <- merge(ad,d,by=c("dataid","clusterid","block"),all.x=T,all.y=T)
dim(d)
dim(ad)

# re-order the tr factor for convenience
ad$tr <- factor(ad$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

# restrict the dataset to a single observation per child 
# (there were 32 cases of 2 index children in 1 compound)
md <- ad[!duplicated(ad[,c('clusterid','dataid','childid')]),]

# drop cases where there was no live birth to get the correct denominator
md <- subset(md,miss1reason!='No live birth')

# identify child deaths
md$death <- ifelse(md$miss1reason=='Child death'|md$miss2reason=='Child death',1,0)

#---------------------------------------
# Cross-tab of deaths by treatment arm
# calculate means and 95% CIs
#---------------------------------------

death.xtab <- table(md$tr,md$death)
death.xtab


arms <- levels(md$tr)
md$svy <- 0
death.mu <- sapply(arms,tmle.mean.est,Y=md$death,tr=md$tr,svy=md$svy,id=md$clusterid,s=0)
death.mu <- t(death.mu)

death.tab <- cbind(rowSums(death.xtab),death.xtab[,2],death.mu)
colnames(death.tab)[1:3] <- c("N at risk","N deaths","Cum Incidence")

#---------------------------------------
# Mantel-Haenszel CIR and RD estimates
# Each intervention arm vs. Control
#---------------------------------------

h1.contrasts <- list(
  c("Control","Water"),
  c("Control","Sanitation"),
  c("Control","Handwashing"),
  c("Control","WSH"),
  c("Control","Nutrition"),
  c("Control","Nutrition + WSH")
)
death.cir <- t(sapply(h1.contrasts,ITT.unadj,Y=md$death,tr=md$tr,strat=md$block,binomial=T))
death.rd<- t(sapply(h1.contrasts,ITT.unadj,Y=md$death,tr=md$tr,strat=md$block,binomial=T,measure="RD"))
rownames(death.cir) <- rownames(death.rd) <- c("Water v C","Sanitation v C","Handwashing v C","WSH v C","Nutrition v C","Nutrition + WSH v C")

# compare the MH RD estimates with the crude estimates to ensure there is nothing to whacky going on with
# the relatively rare outcome and highly stratified data
death.rd.crude <- death.mu[2:7,1]-death.mu[1,1]
round(cbind(death.rd.crude,death.rd[,1], death.rd.crude-death.rd[,1]),5)

# exponentiate the CIR and 95% estimates
death.cir[,c(1,3,4)] <- exp(death.cir[,c(1,3,4)])
colnames(death.cir)[1:2] <- c("CIR","se.logCIR")


#---------------------------------------
# print results
#---------------------------------------

round(death.tab,4)

round(death.rd,4)

round(death.cir,4)


# save everything except the datasets themselves
rm(list=c("d","td","ad","md"))
save.image(file=here("results/raw/","bangladesh-mortality.RData"))




	