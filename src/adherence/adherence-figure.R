#---------------------------------------
# bangladesh-uptake-plot.R
#
# ben arnold (benarnold@berkeley.edu)
#
# plot uptake measures
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls()); library(here)
library(scales)

#---------------------------------------
# load the uptake estimates
#---------------------------------------

load(here("results/bangladesh-uptake.RData"))



#---------------------------------------
# general plotting function for uptake
#---------------------------------------

uplot <- function(x,varlab="",svy,cols,yaxis=F) {
  
  # empty plot
  plot(1:7,1:7,type="n",
       xaxt="n",xlab="",xlim=c(0.5,7.5),
       yaxt="n",ylab="",bty="n",ylim=c(0,1)
  )
  
  # Y-axis
  if(yaxis==TRUE) mtext(seq(0,100,by=20),side=2,line=0,at=seq(0,1,by=0.2),las=1,col="gray20")
  segments(x0=0,x1=8,y0=seq(0,1,by=0.2),col="gray80",lty=2)
  
  # X-axis labels
  mtext(c("C","W","S","H","WSH","N","WSHN"),side=1,line=0,at=c(1:7),col=cols,cex=0.8,las=1)
  #mtext(c("C","W","S","H","WSH","N","WSHN")[c(2,4,6)],side=1,line=1,at=c(2,4,6),col=cols[c(2,4,6)],cex=0.8,las=1)
  mtext(svy,side=1,line=3,col="gray20",cex=1)
  
  # plot points
  arrows(x0=1:7, y0=x[2,], y1=x[3,], col=cols,lwd=1,length=0.05,angle=90,code=3)
  points(1:7,x[1,],pch=21,cex=1.5,lwd=1,col=cols,bg="white")
  points(1:7,x[1,],pch=21,cex=1.5,lwd=0,col=cols,bg=alpha(cols,alpha=0.5))

}

# general label plot
ulabplot <- function(title) {
  plot(1,1,type="n",
       xaxt="n",xlab="",xlim=c(0,1),
       yaxt="n",ylab="",bty="n",ylim=c(0,1)
  )
  text(1,0.5,title,adj=1,cex=1.5)
}



#---------------------------------------
# put the result objects into a list
# to make them easier to plot in a loop
#---------------------------------------

freechl0 <- matrix(NA,nrow=nrow(freechl1),ncol=ncol(freechl1))
rlnsp0 <- matrix(NA,nrow=nrow(rlnsp1),ncol=ncol(rlnsp1))
uptakelist <- list(list(storewat0,storewat1,storewat2),
                   list(freechl0,freechl1,freechl2),
                   list(latseal0,latseal1,latseal2),
                   list(latfeces0,latfeces1,latfeces2),
                   list(hwss0,hwss1,hwss2),
                   list(rlnsp0,rlnsp1,rlnsp2))
# labels for the lists
uptakelabs <- c(
  "Stored drinking\nwater\n(%)",
  "Stored drinking\nwater has\ndetectable\nfree chlorine\n(%)",
  "Latrine with \na functional\nwater seal\n(%)",
  "No visible feces\non latrine slab\nor floor\n(%)",
  "Handwashing\nlocation\nhas soap\n(%)",
  "LNS sachets\nconsumed\n(% of expected)"
)

pdf(here("results/figs/bangladesh-uptake.pdf"),width=10,height=14)
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
lo <- layout(mat=matrix(1:24,ncol=4,nrow=6,byrow=T),widths=c(0.5,1,1,1))
op <- par(mar=c(4,2.5,3,0.5)+0.1)

for(i in 1:length(uptakelist)) {
  if(i!=length(uptakelist)){
    op <- par(mar=c(3,1,1,0.5)+0.1)
    ulabplot(uptakelabs[i])
    op <- par(mar=c(3,2.5,1,0.5)+0.1)
    uplot(x=uptakelist[[i]][[1]],cols=cols,svy="Enrollment",yaxis=T)
    uplot(x=uptakelist[[i]][[2]],cols=cols,svy="Year 1")
    uplot(x=uptakelist[[i]][[3]],cols=cols,svy="Year 2")
  } else{
    op <- par(mar=c(3,1,1,0.5)+0.1)
    ulabplot(uptakelabs[i])
    op <- par(mar=c(3,2.5,1,0.5)+0.1)
    uplot(x=uptakelist[[i]][[1]],cols=cols,svy="",yaxis=T)
    uplot(x=uptakelist[[i]][[2]],cols=cols,svy="")
    uplot(x=uptakelist[[i]][[3]],cols=cols,svy="")
  }
  
}
par(op)
dev.off()








