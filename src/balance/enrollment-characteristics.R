


#---------------------------------------
# enrollment-characteristics.R
#
# ben arnold (benarnold@berkeley.edu)
#
# summarize enrollment characteristics
# by treatment arm
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-enrol.csv
#
# output files:
#	bangladesh-enrol-characteristics-ben.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(dplyr)

#---------------------------------------
# CONSORT: compounds
#---------------------------------------

d <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-enrol.csv")

# re-order the tr factor for convenience
d$tr <- factor(d$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

# number of observations
d.tr <- group_by(d, tr)
ncompounds <- summarise(d.tr,n=n())

# calculate N obs and means for each variable
vlist <- c("momage","momeduy","dadeduy","dadagri","Nhh","elec","cement","landacre","tubewell","storewat","treatwat","odmen","odwom","odch815","odch38","odchu3","latown","latslab","latseal","latfeces","potty","humfeces","humfecesch","hwlatwat","hwlatsoap","hwkitwat","hwkitsoap")

ns <- sapply(d[vlist],function(x) tapply(x,d$tr,function(y) length(y[!is.na(y)])))
mus <- sapply(d[vlist],function(x) tapply(x,d$tr,function(y) mean(y,na.rm=T)))

#---------------------------------------
# combine results into a single dataframe
#---------------------------------------
balance.tab.n <- t(ns)

balance.tab.mu <- t(mus)


table1 <- data.frame(
   variable=rownames(balance.tab.n),
   balance.tab.n[,1],
    balance.tab.mu[,1],
    balance.tab.n[,2],
    balance.tab.mu[,2],
    balance.tab.n[,3],
    balance.tab.mu[,3],
    balance.tab.n[,4],
    balance.tab.mu[,4],
    balance.tab.n[,5],
    balance.tab.mu[,5],
    balance.tab.n[,6],
    balance.tab.mu[,6],
    balance.tab.n[,7],
    balance.tab.mu[,7],
   stringsAsFactors = F
)
names(table1) <- c("variable",paste(rep(colnames(balance.tab.n),rep(2,7)),rep(c(".n",".mu"),7),sep=""))

# print table
table1

table1_b <- table1
rm(table1)

#---------------------------------------
# save objects 
#(drop datasets to save space)
#---------------------------------------
rm(d)
save(table1_b,file="~/dropbox/WBB-primary-analysis/results/raw/ben/bangladesh-enrol-characteristics-ben.RData")




