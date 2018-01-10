capture log close
set more off
clear all


log using "~/WBB-primary-outcomes/src/dm/6-bangladesh-dm-uptake.log", text replace

*--------------------------------------------
* 6-bangladesh-dm-uptake.do
*
* ben arnold (benarnold@berkeley.edu)
*
* process the Bangladesh trial uptake
* data for analysis
*
* note: unlike other dm scripts, this script
* relies on an earlier dm script (3-bangladesh-dm-diar.do)
* to grab child ages at survey measurements
* rather than re-calculate them here
*
*--------------------------------------------




*--------------------------------------------
* input files:
*
*  Treatment assignments 
*  washb-bangladesh-tr.dta (or washb-bangladesh-blind-tr.dta if blinded still)
*
*  1. WASHB_Baseline_main_survey.dta
*  1. WASHB_Midline_main_survey_cleaned.dta
* 04. WASHB_Endline_main_survey_cleaned.dta
* 11. WASHB_Endline_LNS_cleaned.dta
*
*  washb-bangladesh-diar.dta
*
* output files:
*  final/washb-bangladesh-uptake.dta / .csv
*--------------------------------------------


*--------------------------------------------
* grab the treatment assignments
*--------------------------------------------
* blinded treatment assignments:
*use "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-blind-tr.dta", clear

* real treatment assignments:
use "/Volumes/0-Treatment-assignments/washb-bangladesh-tr.dta", clear

sort clusterid
tempfile trdata
save `trdata'


*--------------------------------------------
* load index child ages at each survey
* to subset nutrition uptake measures to 
* children 6-24 months old
*--------------------------------------------
use "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-diar.dta", clear
* restrict to index children and drop twins
keep if childid=="T1"
keep dataid svy agedays ageyrs
sort dataid 
tempfile aged
save `aged'


*--------------------------------------------
* load the enrollment dataset
*--------------------------------------------

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/1_Baseline/1. WASHB_Baseline_main_survey.dta", clear
sort dataid 

* format survey dates
gen svydate = q4002
	format svydate %d
	label var svydate "Survey date"
	codebook svydate

* identify survey round
gen svy = 0
	label var svy "Survey round (0,1,2)"
	
* restrict to the relevant variables for uptake

# delimit;
keep dataid svy svydate
q1003_*
q809_9a q809_9 q809_9a q809_16 
q4201 q4203 q4205
q704_1 q704_2 q704_3 q704_4 q704_5 q704_6
q1009_* q1026 q1027*
;
# delimit cr
order dataid svy svydate
sort dataid

tempfile svy0
save `svy0'


*--------------------------------------------
* load the year 1 dataset
*--------------------------------------------

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/2_Midline/1. WASHB_Midline_main_survey_cleaned.dta", clear

* format survey dates
gen svydate = q4002
	format svydate %d
	label var svydate "Survey date"
	codebook svydate

* identify survey round
gen svy = 1
	label var svy "Survey round (0,1,2)"
	
* restrict to the relevant variables for uptake

# delimit;
keep dataid svy svydate
q1003_*
q809_9a q809_9 q809_9a q809_16 
q4201 q4203 q4205
q704_1 q704_2 q704_3 q704_4 q704_5 q704_6
q1009_* q1026 q1027*
n1402* n1403_* n1404 n1405 n1406 n1407 n1407a n1408 n1409
;
# delimit cr
order dataid svy svydate
sort dataid

* convert LNS indicators from string to real (for clean append, below)
local vlist "n1404 n1405 n1406 n1407 n1408 n1409"
foreach var of local vlist {
	gen r`var' = real(`var')
	drop `var'
	rename r`var' `var'
}


tempfile svy1
save `svy1'



*--------------------------------------------
* load the year 2 LNS uptake dataset
*--------------------------------------------

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/11. WASHB_Endline_LNS_cleaned.dta", clear

* drop 9 twins (seem to have identical data)
drop if childid=="T2"
keep dataid n1402* n1403_* n1404 n1405 n1406 n1407 n1407a n1408 n1409

sort dataid
tempfile y2lns
save `y2lns'


*--------------------------------------------
* load the year 2 main survey dataset
*--------------------------------------------

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/04. WASHB_Endline_main_survey_cleaned.dta", clear

* format survey dates
gen svydate = date(q4002,"DMY")
	format svydate %d
	label var svydate "Survey date"
	codebook svydate

* identify survey round
gen svy = 2
	label var svy "Survey round (0,1,2)"
	
* variable q4205 changed to q4206 in year 2
rename q4206 q4205
	
* restrict to the relevant variables for uptake
# delimit;
keep dataid svy svydate
q1003_*
q809_9a q809_9 q809_9a q809_16 
q4201 q4203 q4205
q704_1 q704_2 q704_3 q704_4 q704_5 q704_6
q1009_* q1026 q1027*
;
# delimit cr

* for handwashing station indicators, missing values should be 0
for any q704_1 q704_2 q704_3 q704_4 q704_5 q704_6: replace X = 0 if X==.

order dataid svy svydate
sort dataid

* merge in the LNS uptake variables
merge 1:1 dataid using `y2lns'
assert _merge !=2
drop _merge


*--------------------------------------------
* Append the enrollment, year 1, year 2 data
*--------------------------------------------
append using `svy0'
append using `svy1'

*--------------------------------------------
* merge in child ages at each measurement
*--------------------------------------------
sort dataid svy
merge 1:1 dataid svy using `aged'
assert _merge !=2
drop _merge

*--------------------------------------------
* merge in the treatment assignment info
*--------------------------------------------
gen clusterid = substr(dataid,1,3)
sort clusterid
capture drop  block tr
merge m:1 clusterid using `trdata'
assert _merge == 3
drop _merge


*--------------------------------------------
* Updtake indicators
*--------------------------------------------


* has stored drinking water
recode q1003_6 q1003_7 q1003_8 (.=0)
gen byte storewat = (q1003_6==1 | q1003_7==1 | q1003_8==1)
	replace storewat = . if (q1003_6==. & q1003_7==. & q1003_8==.)
	label var storewat "Store drinking water (q10001)"
	
* store water with detectable free chlorine (not at baseline)
gen byte freechl = q1027level>0.1 & q1027level<.
	replace freechl = . if (svy==0) | (q1027==999)
	label var freechl "Free chlorine detected in stored water (>0.1 mg/L)"
	
* latrine with a functional water seal
gen byte latseal = (q809_9a==1)
	replace latseal = . if (q809_9a==. | q809_9a==888)
	label var latseal "Latrine has functional water seal (q809_9a)"

* no visible feces on the latrine slab or floor
gen byte latfeces = (q809_16!=1)
	replace latfeces = . if (q809_16==. | q809_16==888)
	label var latfeces "No visible feces on slab/floor of latrine"

* no human feces in house or compound
gen byte humfeces = 1
	replace humfeces = 0 if (q4201>0) | (q4203>0) | (q4205>0)
	replace humfeces = . if inlist(q4201,99,.) & inlist(q4203,99,.) & inlist(q4205,99,.)
	label var humfeces "No human feces observed in house/compound (q4201,4203,4205)"
	
* primary handwashing station has water
gen byte hwsw = (q704_1==1)
	replace hwsw = . if q704_1==.
	label var hwsw "Prim handwashing loc has water (q704_1)"

* primary handwashing station has soap
gen byte hwss = (q704_2==1|q704_3==1|q704_4==1|q704_5==1|q704_6==1)
	replace hwss = . if (q704_2==.&q704_3==.&q704_4==.&q704_5==.&q704_6==.)
	label var hwss "Prim handwashing loc has soap (q704_2-6)"
	
* primary handwashing station has water + soap
gen byte hwsws = (hwsw==1&hwss==1)
	replace hwsws = . if (hwsw==. & hwss==.)
	label var hwsws "Prim handwashing loc has water+soap (q704_1-6)"

	
* mean sachets of LNS fed in prior week to index child 6-24 mos
* using 5.95 months because of 1 child who was 5.98 months due to rounding
gen double agem = agedays/30.4167
	label var agem "Index child, age in months"
	
recode n1404 n1405 n1406 n1407 (99=.) (999=.)

  * Sachets consumed per week
  gen lnsn = (n1405+n1406-n1407) / (n1404) * 7
  replace lnsn = . if (agem < 5.95) | (agem>24)
  label var lnsn "LNS sachets consumed per week"
  
  * percent of expected
  gen lnsp = (lnsn/14)
  label var lnsp "Percent of expected LNS sachets consumed"
  
  list n1404 - n1407 lnsn lnsp if lnsn<0

sum lnsn lnsp, d
count if (lnsn<0 | lnsn>28) & (lnsn!=.)
sum lnsn lnsp if lnsn>=0 & lnsn<=28, d
  
* mean sachets of LNS reported fed in prior week to index child 6-24 mos

	* Sachets consumed per week
	gen rlnsn = (n1408*n1409)
	replace rlnsn = . if (agem < 5.95) | (agem>24)
	label var rlnsn "LNS sachets consumed per week (reported)"
 
	* percent of expected
	gen rlnsp = (rlnsn/14)
	label var rlnsp "Percent of expected LNS sachets consumed (reported)"
	
	list n1408 n1409 rlnsn rlnsp if n1409>2 & n1409<.
	
sum rlnsn rlnsp, d
	
	
*--------------------------------------------
* a single Handwashing compound has 
* a detectible free chlorine measurement,
* which is a data-entry error (was not measured
* in non-water arms.  correct it.
*--------------------------------------------
replace freechl = 0 if dataid=="64204" & svy==1


*--------------------------------------------
* Save an uptake analysis dataset
* compound-survey level observations
*--------------------------------------------
label var dataid "Compound ID"

* restrict to household level variables used in the analysis
* and save the data
keep dataid clusterid block tr svy svydate storewat freechl latseal latfeces humfeces hwsw hwss hwsws *lns*

order dataid clusterid block tr svy svydate storewat freechl latseal latfeces humfeces hwsw hwss hwsws *lns*


*** DROP TREATMENT ASSIGNMENTS to keep the data fully blinded (analyses can merge in treatement)
drop tr

compress
sort dataid svy
label data "Bangladesh uptake analysis dataset (compound-svy obs), created by 6-bangladesh-dm-uptake.do"
saveold "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-uptake.dta", replace version(12)
outsheet using "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-uptake.csv", comma replace

* write a codebook for the dataset
log close
log using "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-uptake-codebook.txt", text replace
desc
codebook, c
codebook
log close


exit

