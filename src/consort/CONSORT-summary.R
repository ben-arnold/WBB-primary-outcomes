


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
#	washb-bangladesh-track-compound-public.csv
# washb-bangladesh-anthro-public.csv
# washb-bangladesh-diar-public.csv
#
# output files:
#	bangladesh-CONSORT.RData
#
#---------------------------------------


#---------------------------------------
# preamble
#---------------------------------------
source(here::here("src/0-config.R"))

#---------------------------------------
# CONSORT: compounds
#---------------------------------------

# read in the compound tracking file
d <- read.csv(here("data/washb-bangladesh-track-compound-public.csv"))

# merge in the treatment assignments
d_tr    <- read.csv(here('data/washb-bangladesh-tr-public.csv'))
d <- left_join(d,d_tr,by=c("clusterid","block"))

# re-order the tr factor for convenience
d$tr <- factor(d$tr,levels=c("Control","Water","Sanitation","Handwashing","WSH","Nutrition","Nutrition + WSH"))

# re-order the miss1reason and miss2reason factors for convenience
d$miss1reason <- factor(d$miss1reason,levels=c("Not lost","Moved away","Absent","Refused","No live birth","Child death"))
d$miss2reason <- factor(d$miss2reason,levels=c("Not lost","Moved away","Absent","Refused","No live birth","Child death"))


# summarize compounds by survey round
c.enrol <- table(d$tr)
c.year1 <- table(d$tr[d$miss1==0])
c.year2 <- table(d$tr[d$miss2==0])
compound_tracking <- rbind(c.enrol,c.year1,c.year2)


# summarize reasons for dropout by the year 1 measurment
c.lost1 <- table(d$miss1reason,d$tr)
compound_lost1 <- rbind(colSums(c.lost1[2:6,]),c.lost1[2:6,])
rownames(compound_lost1)[1] <- "Compounds Lost"

# summarize reasons for dropout by the year 2 measurment
c.lost2 <- table(d$miss2reason,d$tr)
compound_lost2 <- rbind(colSums(c.lost2[2:6,]),c.lost2[2:6,])
rownames(compound_lost2)[1] <- "Compounds Lost"

#---------------------------------------
# CONSORT: index children
# based on anthropometry dataset
#---------------------------------------
chd <- read.csv(here("data/washb-bangladesh-anthro-public.csv"))

chd <- left_join(chd,d_tr,by=c("clusterid","block"))

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
laz_tracking <- rbind(colSums(ch.year1.laz),ch.year1.laz,
                         colSums(ch.year2.laz),ch.year2.laz)
rownames(laz_tracking) <- c("Year 1 total","Year 1 analyzed", "Year 1 missing","Year 2 total","Year 2 analyzed", "Year 2 missing")


#---------------------------------------
# CONSORT:  children <36 mos at enrollment
# based on diarrhea dataset
#---------------------------------------

d36 <- read.csv(here("data/washb-bangladesh-diar-public.csv"))

d36 <- left_join(d36,d_tr,by=c("clusterid","block"))

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
diar_tracking <- rbind(colSums(ch36.year0.diar),ch36.year0.diar,
                       colSums(ch36.year1.diar),ch36.year1.diar,
                       colSums(ch36.year2.diar),ch36.year2.diar)
rownames(diar_tracking) <- c("Enroll total", "Enroll analyzed", "Enroll missing",
                             "Year 1 total", "Year 1 analyzed", "Year 1 missing",
                             "Year 2 total","Year 2 analyzed", "Year 2 missing")

# identify the number of index children by arm and survey
dindex <- d36 %>%
  filter(tchild=="Target child")

n_index <- table(dindex$svy,dindex$tr)


# identify the number of twin pairs by arm
# twins are identified by childid=="T2"
n_twins <- table(d36$svy[d36$childid=="T2"],d36$tr[d36$childid=="T2"])


#---------------------------------------
# Print results
#---------------------------------------
compound_tracking

compound_lost1

compound_lost2

laz_tracking

diar_tracking

n_index

n_twins


#---------------------------------------
# add "_b" suffix to compare with jade
#---------------------------------------
# compound_tracking_b <- compound_tracking
# compound_lost1_b    <- compound_lost1
# compound_lost2_b    <- compound_lost2
# 
# laz_tracking_b   <- laz_tracking
# diar_tracking_b  <- diar_tracking
# n_index_b <- n_index
# n_twins_b <- n_twins
# rm(compound_tracking,compound_lost1,compound_lost2,laz_tracking,diar_tracking,n_index,n_twins)
# 

#---------------------------------------
# save objects 
#(drop datasets to save space)
#---------------------------------------
rm(d,chd,d36,dindex)
save.image(file=here("results/bangladesh-CONSORT.RData"))




