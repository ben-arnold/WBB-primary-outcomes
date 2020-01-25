#---------------------------------
# diar-calendar-time-figure.R
#
# summarize diarrhea prevalence
# by calendar time
#---------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls()); library(here)
library(washb)
library(lubridate)
library(scales)

#---------------------------------------
# Load the analysis dataset
#---------------------------------------
d <- read.csv(here("data/washb-bangladesh-diar-public.csv"))

# merge in the treatment assignments
d_tr    <- read.csv(here("data/washb-bangladesh-tr-public.csv"))
d <- left_join(d,d_tr,by=c("clusterid","block"))

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
# create a calendar time aggregator
# by month
#---------------------------------------
ad$caldate <- as.Date(as.character(ad$svydate),format="%d%b%Y")
ad$my <- floor_date(ad$caldate,"month")


#---------------------------------------
# control monthly means w/ robust 95% SEs
#---------------------------------------
mys <- unique(ad$my)
diar_c <- sapply(mys,
                 function(x) {
                   washb_mean(ad$diar7d[ad$tr=="Control" & ad$my==x],
                              id=ad$block[ad$tr=="Control" & ad$my==x],
                              print=F
                              )
                 }
                 )
diar_c <- data.frame(t(diar_c))
diar_c$month <- as.Date(mys)
colnames(diar_c) <- c("n","mean","sd","se","lb","ub","month")
diar_c <- diar_c[order(diar_c$month),c("month","n","mean","sd","se","lb","ub")]

# summarize dist of obs per month
summary(diar_c$n)

#---------------------------------------
# intervention (all arms) monthly means w/ robust 95% SEs
#---------------------------------------
mys <- unique(ad$my)
diar_i <- sapply(mys,
                 function(x) {
                   washb_mean(ad$diar7d[ad$tr!="Control" & ad$my==x],
                              id=ad$block[ad$tr!="Control" & ad$my==x],
                              print=F
                   )
                 }
)
diar_i <- data.frame(t(diar_i))
diar_i$month <- as.Date(mys)
colnames(diar_i) <- c("n","mean","sd","se","lb","ub","month")
diar_i <- diar_i[order(diar_i$month),c("month","n","mean","sd","se","lb","ub")]

# summarize dist of obs per month
summary(diar_i$n)

#---------------------------------------
# make plot w/ color (no CIs)
#---------------------------------------
# brighter color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"

pdf("results/figs/wbb-diar-calendar-month.pdf",width=6,height=3)
# general plotting parameters
op <- par(mar=c(4,4,1,1)+0.1)
ytics <- c(0,0.05,0.10,0.15)
xtics <- unique(diar_c$month)
xtics2 <- as.Date(c("2013-10-15","2014-05-01","2015-05-01"))
xtics3 <- as.Date(c("2014-03-01","2015-05-01"))
cols <- c(cblack,cblack)
# empty plot
plot(diar_c$month[diar_c$month<="2014-09-01"],
     diar_c$mean[diar_c$month<="2014-09-01"],
     type="n",
     ylim=range(ytics),ylab="Diarrhoea prevalence (%)",yaxt="n",
     xlim=range(xtics),xlab="",xaxt="n",
     bty="n",las=1)
# axis(1,at=xtics,las=2,labels=format(xtics,"%m"),cex.axis=0.6)
axis(1,at=xtics,las=2,labels=NA,cex.axis=0.6)
mtext(month(xtics),side=1,line=0.5,at=xtics,las=1,cex=0.6)

axis(2,at=ytics,las=1,labels=sprintf("%1.0f",ytics*100))
mtext(format(xtics2,"%Y"),at=xtics2,side=1,line=1.5,cex=0.75)
mtext(c("Year 1 follow-up","Year 2 follow-up"),at=xtics3,side=1,line=2.5,cex=0.75)
legend("topright",c("Control","Intervention (all arms)"),ncol=1,lty=c(1,2),col=cols,bty="n",cex=0.75)

# control years 1 and 2 lines + points
lines(diar_c$month[diar_c$month<="2014-09-01"],diar_c$mean[diar_c$month<="2014-09-01"],lty=1,col=cols[1])
# points(diar_c$month[diar_c$month<="2014-09-01"],diar_c$mean[diar_c$month<="2014-09-01"],col=cols[1],pch=19,cex=0.5)

lines(diar_c$month[diar_c$month>"2014-09-01"],diar_c$mean[diar_c$month>"2014-09-01"],lty=1,col=cols[1])
# points(diar_c$month[diar_c$month>"2014-09-01"],diar_c$mean[diar_c$month>"2014-09-01"],col=cols[1],pch=19,cex=0.5)

# intervention years 1 and 2 lines + points
lines(diar_i$month[diar_i$month<="2014-09-01"],diar_i$mean[diar_i$month<="2014-09-01"],lty=2,col=cols[2])
# points(diar_i$month[diar_i$month<="2014-09-01"],diar_i$mean[diar_i$month<="2014-09-01"],col=cols[2],pch=19,cex=0.5)

lines(diar_i$month[diar_i$month>"2014-09-01"],diar_i$mean[diar_i$month>"2014-09-01"],lty=2,col=cols[2])
# points(diar_i$month[diar_i$month>"2014-09-01"],diar_i$mean[diar_i$month>"2014-09-01"],col=cols[2],pch=19,cex=0.5)


par(op)
dev.off()


#---------------------------------------
# make plot w/ color + shaded 95% CIs
#---------------------------------------
pdf("results/figs/wbb-diar-calendar-month-ci.pdf",width=6,height=3)
# general plotting parameters
op <- par(mar=c(4,4,2,1)+0.1)
ytics <- c(0,0.05,0.10,0.15,0.2,0.25)
xtics <- unique(diar_c$month)
xtics2 <- as.Date(c("2013-10-01","2014-05-01","2015-05-01"))
xtics3 <- as.Date(c("2014-03-01","2015-05-01"))
cols <- c(corange,cblue)
# empty plot
plot(diar_c$month[diar_c$month<="2014-09-01"],
     diar_c$mean[diar_c$month<="2014-09-01"],
     type="n",
     ylim=range(ytics),ylab="Diarrhoea prevalence (%)",yaxt="n",
     xlim=range(xtics),xlab="",xaxt="n",
     bty="n",las=1)
  axis(1,at=xtics,las=2,labels=format(xtics,"%m"),cex.axis=0.6)
  axis(2,at=ytics,las=1,labels=sprintf("%1.0f",ytics*100))
  mtext(format(xtics2,"%Y"),at=xtics2,side=1,line=1.5,cex=0.75)
  mtext(c("Year 1 follow-up","Year 2 follow-up"),at=xtics3,side=1,line=2.5,cex=0.75)
  legend("topright",c("Control","Intervention (all arms)"),ncol=1,lty=c(1,2),col=cols,bty="n",cex=0.75)
  
# control shaded 95% CIs
  polygon(x=c(diar_c$month[diar_c$month<="2014-09-01"],
              rev(diar_c$month[diar_c$month<="2014-09-01"])),
          y=c(diar_c$lb[diar_c$month<="2014-09-01"],
              rev(diar_c$ub[diar_c$month<="2014-09-01"])),
          col=alpha(cols[1],alpha=0.2),border=NA)
  polygon(x=c(diar_c$month[diar_c$month>"2014-09-01"],
              rev(diar_c$month[diar_c$month>"2014-09-01"])),
          y=c(diar_c$lb[diar_c$month>"2014-09-01"],
              rev(diar_c$ub[diar_c$month>"2014-09-01"])),
          col=alpha(cols[1],alpha=0.2),border=NA)

# intervention shaded 95% CIs
  polygon(x=c(diar_i$month[diar_i$month<="2014-09-01"],
              rev(diar_i$month[diar_i$month<="2014-09-01"])),
          y=c(diar_i$lb[diar_i$month<="2014-09-01"],
              rev(diar_i$ub[diar_i$month<="2014-09-01"])),
          col=alpha(cols[2],alpha=0.2),border=NA)
  polygon(x=c(diar_i$month[diar_i$month>"2014-09-01"],
              rev(diar_i$month[diar_i$month>"2014-09-01"])),
          y=c(diar_i$lb[diar_i$month>"2014-09-01"],
              rev(diar_i$ub[diar_i$month>"2014-09-01"])),
          col=alpha(cols[2],alpha=0.2),border=NA)

# control years 1 and 2 lines + points
lines(diar_c$month[diar_c$month<="2014-09-01"],diar_c$mean[diar_c$month<="2014-09-01"],lty=1,col=cols[1])
points(diar_c$month[diar_c$month<="2014-09-01"],diar_c$mean[diar_c$month<="2014-09-01"],col=cols[1],pch=19,cex=0.5)

lines(diar_c$month[diar_c$month>"2014-09-01"],diar_c$mean[diar_c$month>"2014-09-01"],lty=1,col=cols[1])
points(diar_c$month[diar_c$month>"2014-09-01"],diar_c$mean[diar_c$month>"2014-09-01"],col=cols[1],pch=19,cex=0.5)

# intervention years 1 and 2 lines + points
lines(diar_i$month[diar_i$month<="2014-09-01"],diar_i$mean[diar_i$month<="2014-09-01"],lty=2,col=cols[2])
points(diar_i$month[diar_i$month<="2014-09-01"],diar_i$mean[diar_i$month<="2014-09-01"],col=cols[2],pch=19,cex=0.5)

lines(diar_i$month[diar_i$month>"2014-09-01"],diar_i$mean[diar_i$month>"2014-09-01"],lty=2,col=cols[2])
points(diar_i$month[diar_i$month>"2014-09-01"],diar_i$mean[diar_i$month>"2014-09-01"],col=cols[2],pch=19,cex=0.5)

par(op)
dev.off()

