capture log close
set more off
clear all


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/5-bangladesh-dm-consort.log", text replace

*--------------------------------------------
* 5-bangladesh-dm-consort.do
*
* ben arnold (benarnold@berkeley.edu)
*
* create summary datasets that track
* compounds and children through the trial
*
* this will enable us to make a CONSORT flow
* diagram of enrollment, randomization, and loss
* and carefully identify the elidible
* population for each analysis
*
*--------------------------------------------


*--------------------------------------------
* input files:
*   3_Endline/01. WASHB_Midline_Endline_data_count_cleaned.dta
*	washb-bangladesh-enrol.dta
*
* output files:
*	washb-bangladesh-track-compound.dta
*--------------------------------------------


*--------------------------------------------
* summarize reasons for withrawal
* at the compound level
*--------------------------------------------
use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/01. WASHB_Midline_Endline_data_count_cleaned.dta", clear

* recode withdrawal reasons to standard categories
* if a household was absent at year 1 but present in year 2
* then code them as not lost
gen miss1 = status_midline!=1
	label var miss1 "Not measured at year 1"

gen miss2 = status_endline!=1
	label var miss2 "Not measured at year 2"

* summary reasons for compounds missing measurements in years 1 and 2
gen miss1reason = 0
	label var miss1reason "Reason for no measurement in year 1"
	label define missreason 0 "Not lost" 1 "No live birth" 2 "Withdrew" 3 "Moved away" 4 "Child death" 5 "Absent"
	label values miss1reason missreason

	replace miss1reason = 1 if (miss1==1) & inlist(reason_midline,"FALSE PREGNANCY","MISCARRIAGE", "STILL BIRTH","ABORTION")
	replace miss1reason = 2 if (miss1==1) & inlist(reason_midline,"REFUSE")
	replace miss1reason = 3 if (miss1==1) & inlist(reason_midline,"MIGRATION OUT")
	replace miss1reason = 4 if (miss1==1) & inlist(reason_midline,"CHILD DEATH")
	replace miss1reason = 5 if (miss1==1) & inlist(reason_midline,"ABSENT")

	
* summary reason for a compound lost from the study population
gen miss2reason = .
	label var miss2reason "Reason for no measurement in year 2"
	label values miss2reason missreason
	
	* start with year 2 reasons
	replace miss2reason = 1 if (miss2==1) & inlist(reason_endline,"FALSE PREGNANCY","MISCARRIAGE", "STILL BIRTH","ABORTION")
	replace miss2reason = 2 if (miss2==1) & inlist(reason_endline,"REFUSE")
	replace miss2reason = 3 if (miss2==1) & inlist(reason_endline,"MIGRATION OUT")
	replace miss2reason = 4 if (miss2==1) & inlist(reason_endline,"CHILD DEATH")
	replace miss2reason = 5 if (miss2==1) & inlist(reason_endline,"ABSENT")
	
	* back fill with year 1 reasons
	replace miss2reason = 1 if (miss2==1) & (miss2reason==.) & inlist(reason_midline,"FALSE PREGNANCY","MISCARRIAGE", "STILL BIRTH","ABORTION")
	replace miss2reason = 2 if (miss2==1) & (miss2reason==.) & inlist(reason_midline,"REFUSE")
	replace miss2reason = 3 if (miss2==1) & (miss2reason==.) & inlist(reason_midline,"MIGRATION OUT")
	replace miss2reason = 4 if (miss2==1) & (miss2reason==.) & inlist(reason_midline,"CHILD DEATH")
	replace miss2reason = 5 if (miss2==1) & (miss2reason==.) & inlist(reason_midline,"ABSENT")
	
	replace miss2reason = 0 if (miss2==0)
	assert miss2reason!=.
	
* there are 2 compounds with identified measurements at year 1, but are absent at year 2.
* final determination of their status is pending, but in the interim set them to not lost at year 1 and absent at year 2
replace miss1 = 0 if inlist(dataid,"02604","03704")
replace miss1reason = 0 if inlist(dataid,"02604","03704")
replace miss2reason = 5 if inlist(dataid,"02604","03704")

	
keep dataid miss1 miss1reason miss2 miss2reason
order dataid miss1 miss1reason miss2 miss2reason
sort dataid
tempfile withd
save `withd'

*--------------------------------------------
* CONSORT flow: compounds
*--------------------------------------------
use "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-enrol.dta", clear
keep dataid clusterid block tr
sort dataid
merge 1:1 dataid using `withd'
assert _merge==3
drop _merge
compress
sort dataid

*** DROP TREATMENT ASSIGNMENTS just to keep this easier
drop tr

label data "WASH Benefits Bangladesh tracking file, compounds"
saveold "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-track-compound.dta", version(12) replace
outsheet using "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-track-compound.csv", comma replace

log close
exit



