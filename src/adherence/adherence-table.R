#---------------------------------------
# adherence-table.R
#
# ben arnold (benarnold@berkeley.edu)
#
# make a table of adherence measures
# at enrollment, year 1, and year 2
#---------------------------------------

#---------------------------------------
# preamble
#---------------------------------------
source(here::here("src/0-config.R"))

#---------------------------------------
# load the uptake estimates
#---------------------------------------

load(here("results/bangladesh-uptake.RData"))



#---------------------------------------
# idiosyncratic estimate formatting function
# takes a matrix with rows for arms and cols for mean, lb, ub
#---------------------------------------

makeci <- function(x) {
  paste(sprintf("%1.0f",x[,1]*100), " (",
        sprintf("%1.0f",x[,2]*100),", ",
        sprintf("%1.0f",x[,3]*100),")",sep="")
}
estformat <- function(x) {
  fx <- makeci(t(x))
  fx[grep("NA",fx)] <- "-"
  names(fx) <- colnames(x)
  fx
}

#---------------------------------------
# format estimates for each indicator
#---------------------------------------

# two adherence measures without measurements in all arms
freechl0 <- matrix(NA,nrow=nrow(freechl1),ncol=ncol(freechl1))
rlnsp0 <- matrix(NA,nrow=nrow(rlnsp1),ncol=ncol(rlnsp1))

# stored water
tab_storewat <- sapply(list(storewat0,storewat1,storewat2),estformat)

# free chlorine
tab_freechl <- sapply(list(freechl0,freechl1,freechl2),estformat)

# latrine with water seal
tab_latseal <- sapply(list(latseal0,latseal1,latseal2),estformat)

# latrine w/ no visible feces
tab_latfeces <- sapply(list(latfeces0,latfeces1,latfeces2),estformat)

# handwasthing location has soap
tab_hwss <- sapply(list(hwss0,hwss1,hwss2),estformat)

# LNS consumption
tab_rlns <-  sapply(list(rlnsp0,rlnsp1,rlnsp2),estformat)

#---------------------------------------
# grab N compounds measured at each time
#---------------------------------------
Ncomp <- matrix(format(uptake_table_b[1,2:22],big.mark=" "),nrow=3,ncol=7)
#---------------------------------------
# collate estimates into a 
# single adherence matrix
#---------------------------------------
adtab <- rbind(
  rep(NA,7),
  Ncomp,
  rep(NA,7),
  t(tab_storewat),
  rep(NA,7),
  t(tab_freechl),
  rep(NA,7),
  t(tab_latseal),
  rep(NA,7),
  t(tab_latfeces),
  rep(NA,7),
  t(tab_hwss),
  rep(NA,7),
  t(tab_rlns)
)
# labels for the adherence indicators
adlabs <- c(
  "Stored drinking water (%)",
  "Stored drinking water has detectable free chlorine (%)",
  "Latrine with a functional water seal (%)",
  "No visible feces on latrine slab or floor (%)",
  "Handwashing location has soap (%)",
  "LNS sachets consumed (% of expected)"
)

rownames(adtab) <- c(
  "Number of compounds measured (N)",c("Enrollment","Year 1","Year 2"),
  adlabs[1],c("Enrollment","Year 1","Year 2"),
  adlabs[2],c("Enrollment","Year 1","Year 2"),
  adlabs[3],c("Enrollment","Year 1","Year 2"),
  adlabs[4],c("Enrollment","Year 1","Year 2"),
  adlabs[5],c("Enrollment","Year 1","Year 2"),
  adlabs[6],c("Enrollment","Year 1","Year 2")
)

#---------------------------------------
# write the matrix to a file
#---------------------------------------
write.csv2(adtab,file=here("results/table-adherence-public.csv"),quote=TRUE)

# having trouble w/ excel opening ; delimited file
# so just output to html
adtab2 <- cbind(rownames(adtab),adtab)
print(xtable(adtab2),
      file=here('results/raw/tables/table-adherence.xls'),
      type='html',include.rownames=FALSE
  
)




