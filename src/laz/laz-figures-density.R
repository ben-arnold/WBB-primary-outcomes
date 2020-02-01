# --------------------------------------
# laz-figures-density.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot distributions by study arm for
# LAZ at follow-up 1 and 2

# --------------------------------------

# --------------------------------------
# input files:
#	bangladesh-laz-unadj-t1.RData
#	bangladesh-laz-unadj-t2.RData
# bangladesh-laz-unadj-permute.RData
# 
#	bangladesh-laz-adj-t1.RData
#	bangladesh-laz-adj-t2.RData
# bangladesh-laz-adj-permute.RData
#
# output files:
#	bangladesh-laz-unadj-t1.pdf
# bangladesh-laz-unadj-t2.pdf
#	bangladesh-laz-adj-t1.pdf
# bangladesh-laz-adj-t2.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

source(here::here("src/0-config.R"))



# --------------------------------------
# load the analysis output
# --------------------------------------
load(here("results/bangladesh-laz-unadj-t1.RData"))
load(here("results/bangladesh-laz-unadj-t2.RData"))
load(here("results/bangladesh-laz-unadj-permute.RData"))


# --------------------------------------
# rename analysis output objects
# for convenience
# --------------------------------------
t1sum  <- laz_t1_n
t1diffh1 <- laz_t1_h1_diff_unadj
t1pvalh1 <- laz_t1_h1_pval_unadj
t1diffh3 <- laz_t1_h3_diff_unadj
t1pvalh3 <- laz_t1_h3_pval_unadj

t2sum  <- laz_t2_n
t2diffh1 <- laz_t2_h1_diff_unadj
t2pvalh1 <- laz_t2_h1_pval_unadj
t2diffh3 <- laz_t2_h3_diff_unadj
t2pvalh3 <- laz_t2_h3_pval_unadj

#---------------------------------------
# load the anthropometry analysis data
#---------------------------------------
d <- read.csv(here("data/washb-bangladesh-anthro-public.csv"))

# merge in the treatment assignments
tr    <- read.csv(here("data/washb-bangladesh-tr-public.csv"))
d <- left_join(d,tr,by=c("clusterid","block"))

# subset the anthropometry to target children (excluding siblings)
dim(d)
d <- subset(d,tchild=="Target child")
dim(d)

# Drop children with extreme LAZ values
table(d$laz_x)
d <- subset(d,laz_x!=1)

# Exclude children with missing data (none)
table(is.na(d$laz))

#---------------------------------------
# YEARS 1 and 2
#---------------------------------------
# subset to the relevant measurement
table(d$svy)
ad1 <- subset(d,svy==1)
dim(ad1)
ad2 <- subset(d,svy==2)
dim(ad2)



# --------------------------------------
# plotting and analysis function 
# for comparison of 2 groups using
# kernel densities
# --------------------------------------
lazdenplot <- function(x,y,namex,namey,dstats,ppval,main,letter,cols,ylab=TRUE,mulab=TRUE) {
	# x : a series of LAZ data from the first group
	# y : a series of LAZ data from the second group
	# namex : string. name for x group
	# namey : string. name for y group
	# dstats : a vector of length 5 with t-test statistics (diff, min, max, t-stat, P)
  # ppval : permutation test p-value
	# main : string. plot title
	# letter : letter for multipanel plots
	# cols  : colors used for plotting
	# ylab : Label the Y-axis?
  # mulab : logical. label group mean dots?
  
	
	# format the difference and 95% CI
	diff <- paste(sprintf("%1.2f",dstats[1])," (",sprintf("%1.2f",dstats[2]),", ",sprintf("%1.2f",dstats[3]),")",sep="")
	
	
	# make the empty plot
	op <- par(xpd=FALSE,mar=c(5,5,6,2)+0.1)
	ytics <- seq(0,0.4,by=0.1)
	xtics <- seq(-5,2,by=1)
	plot(density(x),type="n",
		main="",
		ylim=c(0,0.5),yaxt="n",ylab="",
		xlim=c(-5.1,2),xaxt="n",xlab="",
		las=1,bty="n"
		)
		axis(1,at=xtics,las=1,cex.axis=1.5)
		mtext("Length for Age Z-score (LAZ)",side=1,line=3)
		if(ylab==TRUE) {
		  axis(2,at=ytics,las=1,cex.axis=1.5,lwd=0,lwd.ticks=1)
		  mtext("Kernel Density",side=2,line=3.5)
		}
		mtext(main,side=3,line=2,cex=1.5,adj=0)
		mtext(letter,side=3,line=2,cex=2,adj=0,at=-6,font=2)
	
		# draw shaded regions for stunted and severely stunted (not used)
		minx <- min(xtics)
		miny <- min(ytics)
		maxy <- max(ytics)
		# polygon(x=c(-3,-2,-2,-3),y=c(miny,miny,maxy,maxy),border=NA,col="gray95")
		# polygon(x=c(minx,-3,-3,minx),y=c(miny,miny,maxy,maxy),border=NA,col="gray90")
		# mtext(c("Severe","Stunted"),at=c(-4,-2.5),side=3,line=-0.5,cex=0.8,col="gray30")
	
		# add a marker at -2
		#segments(x0=-2,y0=0,y1=0.25,lty=1,col="gray40")
	
		# draw kernal density distributions
		dx <- density(x)
		dy <- density(y)
		lines(dx,col=cols[1],lwd=3,lty=2)
		lines(dy,col=cols[2],lwd=3)
		
		# label group means
		if(mulab==TRUE) {
			segments(x0=mean(x),y0=miny-0.015,y1=0.05,col="gray40",lty=2)
			segments(x0=mean(y),y0=miny-0.015,y1=0.05,col="gray40",lty=2)
			text(x=mean(c(x,y)),y=0.07,"Group Means",col="gray40",cex=1)
			
		}
		# segments(x0=mean(x),y0=miny-0.015,y1=max(dx$y)+0.05,col="gray40",lty=2)
		# segments(x0=mean(x),x1=mean(x)+0.2,y0=max(dx$y)+0.05,col="gray40",lty=2)
		op <- par(xpd=TRUE)
		points(mean(x),miny-0.015,pch=21,cex=1.75,col=cols[1],bg=alpha(cols[1],alpha=0.5))
		
		# segments(x0=mean(y),y0=miny-0.015,y1=max(dy$y)+0.02,col="gray40",lty=2)
		# segments(x0=mean(y),x1=mean(y)+0.2,y0=max(dy$y)+0.02,col="gray40",lty=2)
		points(mean(y),miny-0.015,pch=21,cex=1.75,col=cols[2],bg=alpha(cols[2],alpha=0.5))
		
		
		# # add labels
		# text(x=mean(x)+0.3,y=max(dx$y)+0.05,paste("Mean (SD) in Control: ",muform(x)),adj=0,cex=0.8,col="gray40")
		# text(x=mean(y)+0.3,y=max(dy$y)+0.02,paste("Mean (SD) in WSH: ",muform(y)),adj=0,cex=0.8,col="gray40")
		
		
		# draw a small table in the upper right
		#txs <- c(-1,0,1.1,2)
		txs <- c(-3.4,-2.6,-1.5,-0.6)
		txs2 <- c(2.2)
		# mtext("LAZ",side=3,line=0,at=txs[1],adj=1)
			mtext(c("","N","Mean","SD"),side=3,line=0,at=txs,cex=1,adj=1)
		mtext(namex,side=3,line=-1.2,at=txs[1],adj=1,col=cols[1])
			mtext(format(length(x),big.mark=","),side=3,line=-1.2,at=txs[2],adj=1,cex=0.9,col=cols[1])	
			mtext(sprintf("%1.2f",mean(x)),side=3,line=-1.2,at=txs[3],adj=1,cex=0.9,col=cols[1])	
			mtext(sprintf("%1.2f",sd(x))  ,side=3,line=-1.2,at=txs[4],adj=1,cex=0.9,col=cols[1])	
		mtext(namey,side=3,line=-2.4,at=txs[1],adj=1,col=cols[2])
			mtext(format(length(y),big.mark=","),side=3,line=-2.4,at=txs[2],adj=1,cex=0.9,col=cols[2])	
			mtext(sprintf("%1.2f",mean(y)),side=3,line=-2.4,at=txs[3],adj=1,cex=0.9,col=cols[2])	
			mtext(sprintf("%1.2f",sd(y))  ,side=3,line=-2.4,at=txs[4],adj=1,cex=0.9,col=cols[2])	
		
		mtext("Diff. (95% CI)",side=3,line=0,at=txs2,adj=1)
			mtext(diff,side=3,line=-2.4,at=txs2,adj=1,col="gray20",cex=0.9)
			
		mtext(paste("t-test p =",sprintf("%1.3f",dstats[5])),side=3,line=-4,at=txs2,adj=1,col="gray20",cex=0.9)
		mtext(paste("permute p =",sprintf("%1.3f",ppval)),side=3,line=-5.5,at=txs2,adj=1,col="gray20",cex=0.9)

		par(op)
}		

# --------------------------------------
#  make a multi-panel density plot - year 1
# --------------------------------------
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

pdf(here("results/figs/bangladesh-laz-unadj-t1.pdf"),width=10,height=20)
# set up a layout
lo <- layout(mat=matrix(1:8,nrow=4,ncol=2,byrow=T))

lazdenplot(x=ad1$laz[ad1$tr=="Control"],y=ad1$laz[ad1$tr=="Water"],namex="Control",namey="Water", dstats=t1diffh1[1,], ppval=t1pvalh1[1], main="Water v. Control",letter="a",cols=cols[c(1,2)],mulab=TRUE)

lazdenplot(x=ad1$laz[ad1$tr=="Control"],y=ad1$laz[ad1$tr=="Sanitation"], namex="Control",namey="Sanitation", dstats=t1diffh1[2,], ppval=t1pvalh1[2], main="Sanitation v. Control",letter="b",cols=cols[c(1,3)],mulab=TRUE,ylab=FALSE)

lazdenplot(x=ad1$laz[ad1$tr=="Control"],y=ad1$laz[ad1$tr=="Handwashing"], namex="Control",namey="Handwashing", dstats=t1diffh1[3,], ppval=t1pvalh1[3], main="Handwashing v. Control",letter="c",cols=cols[c(1,4)],mulab=TRUE)

lazdenplot(x=ad1$laz[ad1$tr=="Control"],y=ad1$laz[ad1$tr=="WSH"],namex="Control",namey="WSH",dstats=t1diffh1[4,],ppval=t1pvalh1[4], main="Combined WSH v. Control",letter="d",cols=cols[c(1,5)],mulab=TRUE,ylab=FALSE)

lazdenplot(x=ad1$laz[ad1$tr=="Control"],y=ad1$laz[ad1$tr=="Nutrition"], namex="Control",namey="Nutrition", dstats=t1diffh1[5,],ppval=t1pvalh1[5], main="Nutrition v. Control",letter="e",cols=cols[c(1,6)],mulab=TRUE)

lazdenplot(x=ad1$laz[ad1$tr=="Control"],y=ad1$laz[ad1$tr=="Nutrition + WSH"], namex="Control",namey="Nutrition + WSH", dstats=t1diffh1[6,], ppval=t1pvalh1[6],main="Nutrition + WSH v. Control",letter="f",cols=cols[c(1,7)],mulab=TRUE,ylab=FALSE)

lazdenplot(x=ad1$laz[ad1$tr=="Nutrition"],y=ad1$laz[ad1$tr=="Nutrition + WSH"], namex="Nutrition",namey="Nutrition + WSH", dstats=t1diffh3[1,],ppval=t1pvalh3[1], main="Nutrition + WSH v. Nutrition",letter="g",cols=cols[c(6,7)],mulab=TRUE)

lazdenplot(x=ad1$laz[ad1$tr=="WSH"],y=ad1$laz[ad1$tr=="Nutrition + WSH"], namex="WSH",namey="Nutrition + WSH", dstats=t1diffh3[2,],ppval=t1pvalh3[2], main="Nutrition + WSH v. WSH",letter="h",cols=cols[c(5,7)],mulab=TRUE,ylab=FALSE)

dev.off()

# --------------------------------------
#  make a multi-panel density plot - year 2
# --------------------------------------
pdf(here("results/figs/bangladesh-laz-unadj-t2.pdf"),width=10,height=20)
# set up a layout
lo <- layout(mat=matrix(1:8,nrow=4,ncol=2,byrow=T))

lazdenplot(x=ad2$laz[ad2$tr=="Control"],y=ad2$laz[ad2$tr=="Water"],namex="Control",namey="Water", dstats=t2diffh1[1,], ppval=t2pvalh1[1], main="Water v. Control",letter="a",cols=cols[c(1,2)],mulab=TRUE)

lazdenplot(x=ad2$laz[ad2$tr=="Control"],y=ad2$laz[ad2$tr=="Sanitation"], namex="Control",namey="Sanitation", dstats=t2diffh1[2,], ppval=t2pvalh1[2], main="Sanitation v. Control",letter="b",cols=cols[c(1,3)],mulab=TRUE,ylab=FALSE)

lazdenplot(x=ad2$laz[ad2$tr=="Control"],y=ad2$laz[ad2$tr=="Handwashing"], namex="Control",namey="Handwashing", dstats=t2diffh1[3,], ppval=t2pvalh1[3], main="Handwashing v. Control",letter="c",cols=cols[c(1,4)],mulab=TRUE)

lazdenplot(x=ad2$laz[ad2$tr=="Control"],y=ad2$laz[ad2$tr=="WSH"],namex="Control",namey="WSH",dstats=t2diffh1[4,],ppval=t2pvalh1[4], main="Combined WSH v. Control",letter="d",cols=cols[c(1,5)],mulab=TRUE,ylab=FALSE)

lazdenplot(x=ad2$laz[ad2$tr=="Control"],y=ad2$laz[ad2$tr=="Nutrition"], namex="Control",namey="Nutrition", dstats=t2diffh1[5,],ppval=t2pvalh1[5], main="Nutrition v. Control",letter="e",cols=cols[c(1,6)],mulab=TRUE)

lazdenplot(x=ad2$laz[ad2$tr=="Control"],y=ad2$laz[ad2$tr=="Nutrition + WSH"], namex="Control",namey="Nutrition + WSH", dstats=t2diffh1[6,], ppval=t2pvalh1[6],main="Nutrition + WSH v. Control",letter="f",cols=cols[c(1,7)],mulab=TRUE,ylab=FALSE)

lazdenplot(x=ad2$laz[ad2$tr=="Nutrition"],y=ad2$laz[ad2$tr=="Nutrition + WSH"], namex="Nutrition",namey="Nutrition + WSH", dstats=t2diffh3[1,],ppval=t2pvalh3[1], main="Nutrition + WSH v. Nutrition",letter="g",cols=cols[c(6,7)],mulab=TRUE)

lazdenplot(x=ad2$laz[ad2$tr=="WSH"],y=ad2$laz[ad2$tr=="Nutrition + WSH"], namex="WSH",namey="Nutrition + WSH", dstats=t2diffh3[2,],ppval=t2pvalh3[2], main="Nutrition + WSH v. WSH",letter="h",cols=cols[c(5,7)],mulab=TRUE,ylab=FALSE)

dev.off()

