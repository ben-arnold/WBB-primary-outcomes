
capture log close
set more off
clear all

log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/1-format-tr-assignments.log", text replace

*--------------------------------------------
* 1-format-tr-assignments.do
*
* ben arnold (benarnold@berkeley.edu)
*
*
* input REAL treatment assignments from the
* bangladesh trial and save a formatted
* stata dataset and csv file
*
* need to mount the encrypted disk to run
* this script
*--------------------------------------------

*--------------------------------------------
* input files:
*  washb-bangladesh-raw-tr-assignments.csv
*
* output files:
*  washb-bangladesh-tr.dta /.csv
* 
*--------------------------------------------
insheet using "/Volumes/0-Treatment-assignments/washb-bangladesh-raw-tr-assignments.csv", comma clear
gen clusterid = string(cluster_id)
	replace clusterid = "00"+string(cluster_id) if cluster_id<10
	replace clusterid = "0"+string(cluster_id) if cluster_id>=10 & cluster_id<100
	label var clusterid "Cluster ID"
label var block "Randomization block"
gen tr = 0
	label define tr  1 "Control" 2 "Water" 3 "Sanitation" 4 "Handwashing" 5 "WSH" 6 "Nutrition" 7  "Nutrition + WSH"
	label values tr tr
	label var tr "Randomized treatment assignment"
	replace tr = 1 if tr_label == "Control"
	replace tr = 2 if tr_label == "Water"
	replace tr = 3 if tr_label == "Sanitation"
	replace tr = 4 if tr_label == "Handwashing"
	replace tr = 5 if tr_label == "WSH"
	replace tr = 6 if tr_label == "Nutrition"
	replace tr = 7 if tr_label == "Nutrition + WSH"
	assert tr!=. 
	tab tr
keep block clusterid tr
order block clusterid tr
sort block clusterid

label data "WASH Benefits Bangladesh cluster level treatment assignments"
saveold "/Volumes/0-Treatment-assignments/washb-bangladesh-tr.dta", v(12) replace
outsheet using "/Volumes/0-Treatment-assignments/washb-bangladesh-tr.csv", comma replace

log close
exit
