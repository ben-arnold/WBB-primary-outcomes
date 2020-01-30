
#---------------------------------------
# uptake.R
#
# ben arnold (benarnold@berkeley.edu)
#
# summarize measures of uptake / compliance
# by study arm and measurement round
# (enrollment, year 1, year 2)
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-uptake-public.csv
#
# output files:
# bangladesh-uptake.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
source(here::here("src/0-config.R"))

# source the base functions
source(here("src/basefns/washb-base-functions.R"))

#---------------------------------------
# load the uptake analysis dataset
#---------------------------------------
d <- read.csv(here("data/washb-bangladesh-uptake-public.csv"))

# merge in the treatment assignments
d_tr    <- read.csv(here("data/washb-bangladesh-tr-public.csv"))
d <- left_join(d,d_tr,by=c("clusterid","block"))


# re-order the treatment factor for convenience
d$tr <- factor(d$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

#---------------------------------------
# for each uptake indicator, summarize
# the number of obs and the % at each
# measurement round
#---------------------------------------

# list the indicators to include
inds <- c("storewat","freechl","latseal","latfeces","humfeces","hwsw","hwss","hwsws","rlnsp")

# number of observations
ncompounds <- group_by(d,tr,svy)  %>%
  summarise(n=n())

# reshape long to calculate means and Ns for selected indicators
dlong <- d %>%
  gather(indicator,value,-dataid,-clusterid,-block,-tr,-svy) %>%
  filter(indicator %in% inds ) %>%
  group_by(indicator,tr,svy)

dsum <- dlong %>%
  summarize(mean=mean(value,na.rm=T),N=sum(!is.na(value)),n=sum(value,na.rm=T)) %>%
  ungroup() %>%
  mutate(indicator = factor(indicator,levels=inds)) %>%
  arrange(indicator,tr,svy)

# note that the "n" does not make sense for the LNS sachet consumption
# since that is % of expected consumption
dsum <- dsum %>%
  mutate(n=ifelse(indicator=="rlnsp",NA,n))


# format the indicator labels
dsum$indlab <- ""
dsum$indlab[dsum$indicator %in% "storewat"] <-  "Store water"
dsum$indlab[dsum$indicator %in% "freechl"] <- "Store water with detectable free chlorine"
dsum$indlab[dsum$indicator %in% "latseal"] <- "Latrine with functional water seal"
dsum$indlab[dsum$indicator %in% "latfeces"] <- "Visible feces on latrine slab or floor"
dsum$indlab[dsum$indicator %in% "humfeces"] <- "Human feces in house or compound"
dsum$indlab[dsum$indicator %in% "hwsw"] <- "Primary handwashing station has water"
dsum$indlab[dsum$indicator %in% "hwss"] <- "Primary handwashing station has soap"
dsum$indlab[dsum$indicator %in% "hwsws"] <- "Primary handwashing station has soap & water"
dsum$indlab[dsum$indicator %in% "rlnsp"] <- "LNS sachet consumption"

# now spread the data wide by arm
# for final tables

uptake_tab_n <- dsum %>%
  select(-mean,-N) %>%
  spread(tr,n)

uptake_tab_mean <- dsum %>%
  select(-N,-n) %>%
  spread(tr,mean)

# print tables
data.frame(uptake_tab_n)
data.frame(uptake_tab_mean)


#---------------------------------------
# Calculate means and influence-curve 
# based 95% CIs by survey round
#---------------------------------------
arms <- levels(d$tr)

# store water
storewat0 <- sapply(arms,tmle.mean.est,Y=d$storewat,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
storewat1 <- sapply(arms,tmle.mean.est,Y=d$storewat,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
storewat2 <- sapply(arms,tmle.mean.est,Y=d$storewat,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# store water with detectable chlorine (not measured at enrollment)
# note: since there were no events in any of the non-water arms (actually, 1 in HW)
# we cannot estimate 95% CIs. Will pad the matrix with zeros and NAs
# freechl0 <- sapply(arms,tmle.mean.est,Y=d$freechl,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
warms <- c("Water","WSH","Nutrition + WSH")
freechl1 <- sapply(warms,tmle.mean.est,Y=d$freechl,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
  # pad with zeros and missings for other treatment arms
  freechl1 <- cbind(c(0,NA,NA),freechl1[,1],c(0,NA,NA),c(0,NA,NA),freechl1[,2],c(0,NA,NA),freechl1[,3])
  colnames(freechl1) <- colnames(storewat1)
  
freechl2 <- sapply(warms,tmle.mean.est,Y=d$freechl,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)
  # pad with zeros and missings for other treatment arms
  freechl2 <- cbind(c(0,NA,NA),freechl2[,1],c(0,NA,NA),c(0,NA,NA),freechl2[,2],c(0,NA,NA),freechl2[,3])
  colnames(freechl2) <- colnames(storewat2)

# Latrine w/ a functional water seal
latseal0 <- sapply(arms,tmle.mean.est,Y=d$latseal,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
latseal1 <- sapply(arms,tmle.mean.est,Y=d$latseal,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
latseal2 <- sapply(arms,tmle.mean.est,Y=d$latseal,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Visible feces on the latrine slab or floor
latfeces0 <- sapply(arms,tmle.mean.est,Y=d$latfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
latfeces1 <- sapply(arms,tmle.mean.est,Y=d$latfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
latfeces2 <- sapply(arms,tmle.mean.est,Y=d$latfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Human feces in the compound
humfeces0 <- sapply(arms,tmle.mean.est,Y=d$humfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
humfeces1 <- sapply(arms,tmle.mean.est,Y=d$humfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
humfeces2 <- sapply(arms,tmle.mean.est,Y=d$humfeces,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Primary handwashing station has water
hwsw0 <- sapply(arms,tmle.mean.est,Y=d$hwsw,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
hwsw1 <- sapply(arms,tmle.mean.est,Y=d$hwsw,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
hwsw2 <- sapply(arms,tmle.mean.est,Y=d$hwsw,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Primary handwashing station has soap
hwss0 <- sapply(arms,tmle.mean.est,Y=d$hwss,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
hwss1 <- sapply(arms,tmle.mean.est,Y=d$hwss,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
hwss2 <- sapply(arms,tmle.mean.est,Y=d$hwss,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Primary handwashing station has soap & water
hwsws0 <- sapply(arms,tmle.mean.est,Y=d$hwsws,tr=d$tr,svy=d$svy,id=d$clusterid,s=0)
hwsws1 <- sapply(arms,tmle.mean.est,Y=d$hwsws,tr=d$tr,svy=d$svy,id=d$clusterid,s=1)
hwsws2 <- sapply(arms,tmle.mean.est,Y=d$hwsws,tr=d$tr,svy=d$svy,id=d$clusterid,s=2)

# Mean sachets of LNS fed in prior week to index child 6-24 mos (not measured at enrollment, only measured in nutrition arms)
narms <- arms[grep("Nutrition",arms)]
rlnsp1 <- sapply(narms,tmle.mean.est,Y=d$rlnsp,tr=d$tr,svy=d$svy,id=d$clusterid,s=1,family="gaussian")
  # pad with missings for other treatment arms
  rlnsp1 <- cbind(matrix(NA,nrow=3,ncol=5),rlnsp1)
  colnames(rlnsp1) <- colnames(hwsws1)

rlnsp2 <- sapply(narms,tmle.mean.est,Y=d$rlnsp,tr=d$tr,svy=d$svy,id=d$clusterid,s=2,family="gaussian")
  # pad with missings for other treatment arms
  rlnsp2 <- cbind(matrix(NA,nrow=3,ncol=5),rlnsp2)
  colnames(rlnsp2) <- colnames(hwsws2)



#---------------------------------------
# save the objects
#---------------------------------------
rm(d)
save.image(file=here("results/bangladesh-uptake.RData"))




