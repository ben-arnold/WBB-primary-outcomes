# --------------------------------------
# lazminus2-figure.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot the WASH B bangladesh stunting
# results
#
#
# --------------------------------------

# --------------------------------------
# input files:
#	bangladesh-lazminus2-N-prev.RData
#   bangladesh-lazminus2-unadj.RData
#
# output files:
#	  bangladesh-lazminus2.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)

# --------------------------------------
# load the analysis output files
# --------------------------------------

load('~/Dropbox/wbb-primary-analysis/results/raw/ben/bangladesh-lazminus2-N-prev-ben.RData')
load('~/Dropbox/wbb-primary-analysis/results/raw/ben/bangladesh-lazminus2-unadj-t2-ben.RData')

# --------------------------------------
# format the objects for plotting
# --------------------------------------

glab <- c("Control\n","Water\n","Sanitation\n","Handwashing\n","Combined\nWSH","Nutrition\n","Combined\nNutrion+WSH")
glab2 <- c("C","W","S","H","WSH","N","WSHN")


fuPrev <- stunt_t2_prev_b
fuN <- stunt_t2_n_b

h1rd <- stunt_t2_h1_rd_unadj_b
h3rd <- stunt_t2_h3_rd_unadj_b


pdf("~/dropbox/wbb-primary-analysis/results/figs/bangladesh-lazminus2.pdf",width=10,height=4)

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cols <- c("gray30",cbPalette[c(2:4,6:8)])
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
cols=c(cblack,cblue,cteal,cgreen,corange,cred,cmagent)

ytics <- seq(0,50,by=10)
op <- par(mar=c(1,9,9,0)+0.1,xpd=TRUE)
# set up an empty plot
MidPts <- barplot(1:7,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
	)
	segments(x0=0,x1=max(MidPts+0.5),y0=ytics,lty=1,lwd=1,col="gray90")
	axis(2,at=ytics,las=1)
	mtext("LAZ < -2\nPrevalence\nat Year 2\n(%)",side=2,line=3,las=1)

	# plot estimates
	arrows(x0=MidPts, y0=fuPrev[,2]*100, y1=fuPrev[,3]*100, col=cols,lwd=2,length=0.05,angle=90,code=3)
	points(MidPts,fuPrev[,1]*100,pch=21,cex=2,lwd=2,col=cols,bg="white")
	points(MidPts,fuPrev[,1]*100,pch=21,cex=2,lwd=0,col=cols,bg=alpha(cols,alpha=0.5))
	text(x=MidPts+0.05,y=fuPrev[,1]*100,sprintf("%1.1f",fuPrev[,1]*100),pos=4,cex=1,col=cols,font=1)
	
	# print header and footer labels
	mtext(glab,at=MidPts,side=3,line=6,col=cols,font=1  )
	hx <- MidPts[1]-2.5
	hx2 <- MidPts[1]-2.6
	rdform <- function(rd,lb,ub) {
		paste(sprintf("%1.1f",rd*100)," (",sprintf("%1.1f",lb*100),", ",sprintf("%1.1f",ub*100),")",sep="")
	}
	
	# print Ns in the footer
	mtext(c("N =",fuN[,1]),side=1,line=0,at=c(0,MidPts),cex=0.8,col=c(cols[1],cols))
	
	# print header table - RDs for H1	
	mtext("Difference (95% CI)",side=3,line=5.5,at=hx2,adj=0,cex=1)
	mtext("Intervention v. Control",side=3,line=4.5,at=hx,adj=0,cex=0.8,col="gray30")
	mtext(c("ref",rdform(h1rd[,1],h1rd[,2],h1rd[,3])),side=3,line=4.5,at=MidPts,cex=0.8,col="gray30")
	
	# print header table - RDs for H3
	mtext("Nutrition + WSH v. WSH",side=3,line=3,at=hx,adj=0,cex=0.8,col="gray30")
	mtext(c("ref",rdform(h3rd[2,1],h3rd[2,2],h3rd[2,3])),side=3,line=3,at=MidPts[c(5,7)],cex=0.8,col=c(cols[5],"gray30"))
	
	mtext("Nutrition + WSH v. Nutrition",side=3,line=2,at=hx,adj=0,cex=0.8,col="gray30")
	mtext(c("ref",rdform(h3rd[1,1],h3rd[1,2],h3rd[1,3])),side=3,line=2,at=MidPts[c(6,7)],cex=0.8,col=c(cols[6],"gray30"))
	
	
par(op)
dev.off()




