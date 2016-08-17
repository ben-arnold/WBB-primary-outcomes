

#-------------------------------------
# 1-re-randomize-tr-assignments.R
#
# ben arnold
#
# scramble the treatment assignments
# to make a dataset that other team
# members can use for blinded analyses
# of the WASH Benefits Bangladesh trial
#
# this scrambled assignment will be
# included in the washb R package
#-------------------------------------

rm(list=ls())
library(foreign)
d <- read.dta('~/dropbox/wbb-primary-analysis/data/final/ben/washb-bang-tr.dta')

# randomly re-order treatment within each block
set.seed(6435345)
ru <- runif(nrow(d))
trr <- d$tr[order(d$block,ru)]

# ensure there is no correlation between actual treatment assignments and re-randomized ones
chisq.test(d$tr,trr)

# make a re-randomized dataset
blind_tr <- d
blind_tr$tr <- trr
attr(blind_tr,'var.labels')[3] <- "Randomized treatment assignment (scrambled/blinded!!)"
attr(blind_tr,'datalabel') <- 'WASH Benefits Bangladesh cluster level treatment assignments (scrambled/blinded!!)'

# save R version of the TRUE tr assignments
save(d,file='/Volumes/0-Treatment-assignments/washb-bang-tr.RData')

# save a re-randomized dataset
save(blind_tr,file='~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bang-blind-tr.RData')
write.dta(blind_tr,file='~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bang-blind-tr.dta')
write.csv(blind_tr,file='~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bang-blind-tr.csv')
