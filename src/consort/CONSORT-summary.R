


#---------------------------------------
# CONSORT-summary.R
#
# ben arnold (benarnold@berkeley.edu)
#
# summarize participant flow by 
# treatment arm
#---------------------------------------

#---------------------------------------
# input files:
#	washb-bangladesh-track-compound.csv
# washb-bangladesh-anthro.csv
# washb-bangladesh-diar.csv
#
# output files:
#	bangladesh-CONSORT-ben.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())


#---------------------------------------
# CONSORT: compounds
#---------------------------------------

d <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-track-compound.csv")

# re-order the tr factor for convenience
d$tr <- factor(d$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

# re-order the miss2reason factor for convenience
#d$miss2reason <- factor(d$miss2reason,levels=c("Not lost","Moved away","Absent","Withdrew","No live birth","Child death"))
d$miss2reason <- factor(d$miss2reason,levels=c("Not lost","No live birth","Child death","Moved away","Withdrew","Absent"))

# summarize compounds by survey round
c.enrol <- table(d$tr)
c.year1 <- table(d$tr[d$miss1==0])
c.year2 <- table(d$tr[d$miss2==0])
compound_tracking <- rbind(c.enrol,c.year1,c.year2)


# summarize reasons for dropout by the year 2 measurment
c.lost <- table(d$miss2reason,d$tr)
compound_lost <- rbind(colSums(c.lost[2:6,]),c.lost[2:6,])
rownames(compound_lost)[1] <- "Compounds Lost"

# confirm that the reasons compounds were lost sum to the differences
(c.enrol - compound_lost[1,] ) == c.year2



#---------------------------------------
# CONSORT: index children
# based on anthropometry dataset
#---------------------------------------
chd <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-anthro.csv")

# re-order the tr factor for convenience
chd$tr <- factor(chd$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

# identify missing observations for anthro measurements
chd$miss.laz <- ifelse(is.na(chd$laz),1,0)
chd$miss.waz <- ifelse(is.na(chd$waz),1,0)
chd$miss.whz <- ifelse(is.na(chd$whz),1,0)
chd$miss.hcz <- ifelse(is.na(chd$hcz),1,0)

# Summary of analyzed and missing values by year of follow-up
ch.year1.laz <- table(chd$miss.laz[chd$svy==1],chd$tr[chd$svy==1])
ch.year2.laz <- table(chd$miss.laz[chd$svy==2],chd$tr[chd$svy==2])

# summary object
tchild_tracking <- rbind(ch.year1.laz[1,],
                         ch.year2.laz[1,],
                         ch.year1.laz[2,],
                         ch.year2.laz[2,])
rownames(tchild_tracking) <- c("Year 1 analyzed", "Year 2 analyzed","Year 1 missing", "Year 2 missing")


#---------------------------------------
# CONSORT:  children <36 mos at enrollment
# based on diarrhea dataset
#---------------------------------------

d36 <- read.csv("~/dropbox/WBB-primary-analysis/data/final/ben/washb-bangladesh-diar.csv")

# re-order the tr factor for convenience
d36$tr <- factor(d36$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

# drop sibling newbirths
table(d$sibnewbirth)
d36 <- subset(d36,sibnewbirth==0)

# drop any children over 36 mos at enrollment
table(d36$gt36mos)
d36 <- subset(d36,gt36mos==0)

# identify missing observations for diarrhea measurements
d36$miss.diar <- ifelse(is.na(d36$diar7d),1,0)

# Summary of analyzed and missing values by year of follow-up
ch36.year0.diar <- table(d36$miss.diar[d36$svy==0],d36$tr[d36$svy==0])
ch36.year1.diar <- table(d36$miss.diar[d36$svy==1],d36$tr[d36$svy==1])
ch36.year2.diar <- table(d36$miss.diar[d36$svy==2],d36$tr[d36$svy==2])

# summary object
ch36_tracking <- rbind(ch36.year0.diar[1,],
                       ch36.year1.diar[1,],
                       ch36.year2.diar[1,],
                       ch36.year0.diar[2,],
                       ch36.year1.diar[2,],
                       ch36.year2.diar[2,])
rownames(ch36_tracking) <- c("Enroll analyzed","Year 1 analyzed", "Year 2 analyzed","Enroll missing","Year 1 missing", "Year 2 missing")

#---------------------------------------
# Print results
#---------------------------------------
compound_tracking

compound_lost

tchild_tracking

ch36_tracking

#---------------------------------------
# add "_b" suffix to compare with jade
#---------------------------------------
compound_tracking_b <- compound_tracking
compound_lost_b     <- compound_lost
tchild_tracking_b   <- tchild_tracking
ch36_tracking_b     <- ch36_tracking

rm(compound_tracking,compound_lost,tchild_tracking,ch36_tracking)


#---------------------------------------
# save objects 
#(drop datasets to save space)
#---------------------------------------
rm(d,chd,d36)
save.image(file="~/dropbox/WBB-primary-analysis/results/raw/ben/bangladesh-CONSORT-ben.RData")




