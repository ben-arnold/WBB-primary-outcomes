capture log close
set more off
clear all


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/3-bangladesh-dm-diar.log", text replace

*--------------------------------------------
* 3-bangladesh-dm-diar.do
*
* ben arnold (benarnold@berkeley.edu)
*
* process the Bangladesh trial datasets
* to create a diarrhea analysis dataset
* including data from enrollment, year 1 and year 2
*
* The diarrhea primary analysis includes all measurements
* for children <36 months at enrollment
*
* 
*--------------------------------------------




*--------------------------------------------
* input files:
*
*  Treatment assignments 
*	washb-bang-tr.dta (or washb-bang-blind-tr.dta if blinded)
*
*  ENROLLMENT
*  1. WASHB_Baseline_main_survey.dta
*  2. WASHB_Baseline_census.dta
*  3. WASHB_Baseline_childinfo.dta
*  4. WASHB_Baseline_childhealthinfo.dta
*
*  YEAR 1
*  1. WASHB_Midline_main_survey_cleaned.dta
*  3. WASHB_Midline_childinfo_cleaned.dta
*  4. WASHB_Midline_childhealthinfo_cleaned.dta
*
*  YEAR 2
*  04. WASHB_Endline_main_survey_cleaned.dta
*  07. WASHB_Endline_childinformation_cleaned.dta
*  08. WASHB_Endline_childhealthinfo_cleaned.dta
*
* output files:
*  final/washb-bangladesh-diar.dta / .csv
*
*--------------------------------------------



*--------------------------------------------
* grab the treatment assignments
*--------------------------------------------
* blinded treatment assignments:
use "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bang-blind-tr.dta", clear

* real treatment assignments (not used to keep data blinded)
* use "/Volumes/0-Treatment-assignments/washb-bang-tr.dta", clear

sort clusterid
tempfile trdata
save `trdata'



*--------------------------------------------
* ENROLLMENT
* merge the household and child data
*--------------------------------------------


use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/1_Baseline/1. WASHB_Baseline_main_survey.dta", clear
sort dataid 
tempfile bmain
save `bmain'

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/1_Baseline/3. WASHB_Baseline_childinfo.dta", clear
sort dataid
tempfile bchildinfo
save `bchildinfo'

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/1_Baseline/4. WASHB_Baseline_childhealthinfo.dta", clear
sort dataid childid
tempfile bchildhealth
save `bchildhealth'

* merge child information to the household information
use `bmain', clear
merge 1:m dataid using `bchildinfo'
assert _merge !=2
rename _merge _merge1

* merge child health information to the data
sort dataid childid
merge 1:1 dataid childid using `bchildhealth'
assert _merge!=2
rename _merge _merge2

* list 2 siblings who were enrolled but do not have records in the child health info dataset:
tab _merge1 _merge2
list dataid childid dob if _merge1==3 & _merge2==1
drop _merge1 _merge2


sort dataid childid
gen svy = 0
tempfile enrol
save `enrol'



*--------------------------------------------
* YEAR 1
* merge the household and child data
*--------------------------------------------


use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/2_Midline/1. WASHB_Midline_main_survey_cleaned.dta", clear
sort dataid 
tempfile main
save `main'

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/2_Midline/3. WASHB_Midline_childinfo_cleaned.dta", clear
sort dataid hhid
tempfile childinfo
save `childinfo'

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/2_Midline/4. WASHB_Midline_childhealthinfo_cleaned.dta", clear
sort dataid childid
tempfile childhealth
save `childhealth'


use `main', clear
merge 1:m dataid using `childinfo'
assert _merge!=1
drop if _merge==2
drop _merge

* merge in the child health information
sort dataid childid
merge 1:1 dataid childid using `childhealth'

* list 9 children who are not matching, keep the empty illness records for the missing data tally
list dataid hhid childid dob if _merge==1

* there is one child in the illness records file that has no child info. drop him/her
list dataid hhid childid dob if _merge==2
drop if _merge==2
drop _merge

sort dataid childid
gen svy = 1
tempfile year1
save `year1'



*--------------------------------------------
* YEAR 2
* merge the household and child data
*--------------------------------------------

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/04. WASHB_Endline_main_survey_cleaned.dta", clear
sort dataid 
tempfile main
save `main'

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/07. WASHB_Endline_childinformation_cleaned.dta", clear
sort dataid hhid
tempfile childinfo
save `childinfo'

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/08. WASHB_Endline_childhealthinfo_cleaned.dta", clear
sort dataid childid
tempfile childhealth
save `childhealth'


use `main', clear
merge 1:m dataid using `childinfo'
assert _merge==3
drop _merge

* merge in the child health information
sort dataid childid
merge 1:1 dataid childid using `childhealth'

* list 3 children who are not matching, keep the empty illness records for the missing data tally
list dataid hhid childid dob if _merge==1
drop _merge


*--------------------------------------------
* Year 2 data formats differ for many variables, 
* which creates problems in the append. 
* standardize to earlier year formats
*--------------------------------------------

* convert survey date to date format
gen rq4002 = date(q4002,"DMY")
drop q4002
rename rq4002 q4002

* numeric to string
local rtoslist "q4014 q4015 q4016 q4017 q105years q105months q105days q716mins q716secs q721mins q721secs q817years q817months q905adays q905ahours q905amins q912adays q912ahours q912amins q1004days q1004hours q1004mins q1011mins q1015mins q807cahours q807camins q807cbhours q807cbmins q807cchours q807ccmins q807cdhours q807cdmins q807cehours q807cemins"
foreach var of local rtoslist {
	gen r`var' = string(`var')
	drop `var'
	rename r`var' `var'
}

* string to numeric
local storlist "q102 q1303 q1307 q1308 q1310 q1311"
foreach var of local storlist {
	gen r`var' = real(`var')
	drop `var'
	rename r`var' `var'
}

*--------------------------------------------
* APPEND ENROLLMENT, YEAR 1, YEAR 2
*--------------------------------------------
gen svy = 2
	label var svy "Survey round (0,1,2)"
	
append using `enrol'
append using `year1'


*--------------------------------------------
* merge in the treatment assignment information 
*--------------------------------------------

* drop clusterid and block variables to 
* ensure we get a complete set

* merge in the treatment assignment info (keep only matching obs)
gen clusterid = substr(dataid,1,3)
sort clusterid
capture drop  block tr
merge m:1 clusterid using `trdata'
assert _merge == 3
drop _merge


*--------------------------------------------
* identify target children
*--------------------------------------------
gen byte tchild = strpos(childid,"T")>0
	label var tchild "Target child in birth cohort"
	label define tchild 0 "Sibling" 1 "Target child"
	label values tchild tchild
	tab tchild, m

* cross-tab of child observations by round
tab svy tchild, m


*--------------------------------------------
* Age, using survey date and birth date
*--------------------------------------------
label var dob "Date of birth"
format dob %d

gen svydate = q4002
	format svydate %d
	label var svydate "Survey date"
	
	* replace svydate using the entrydt variable for enrollment and year 1 surveys
	* the entrydt variable is more accurate, but is not in the year 2 dataset
	replace svydate = entrydt if inlist(svy,0,1) & (entrydt!=.)
	codebook svydate
	
* Month of measurement
gen month = month(svydate)
	label var month "Month of measurement"

* Age in days and years
gen agedays = svydate-dob
	label var agedays "Age in days"
gen ageyrs = agedays/365.25
	label var ageyrs "Age in years"
codebook ageyrs
bysort tchild svy: sum ageyrs, d



*--------------------------------------------
* Identify new births using the time between
* survey date and enrollment
*--------------------------------------------

gen svydate0 = svydate if svy==0
gen svydate1 = svydate if svy==1

sort dataid
by dataid: egen _date0 = min(svydate0)
by dataid: egen _date1 = min(svydate1)

gen double enrolage = (_date0 - dob)/365.25
label var enrolage "Age at enrollment, years"

gen byte newbirth = enrolage<0
	label var newbirth "New birth"

* create an indicator for non-target new births for easy identification
gen byte sibnewbirth = (tchild==0) & (newbirth==1)
	label var sibnewbirth "Sibling (non-index child) new birth"


*--------------------------------------------
* Identify children who were >36 months 
* at enrollment
*--------------------------------------------
	
* there was one child where enrollment was on his/her birthday, 
* and was exactly 3 years old, so screen >3.01 to retain this child
gen byte gt36mos = enrolage>3.01 & enrolage<.
	label var gt36mos "Older than 36 mos at enrollment"
tab svy gt36mos
sort dataid childid svy
* list dataid childid svy svydate dob enrolage gt36mos if enrolage>3 & enrolage<.

*--------------------------------------------
* Diarrhea, defined as 3+ loose/watery stools in 24 hours or 1+ stool with blood
* calculate separate 2-day recall window and 7-day recall window
* include the (partial) interview day in the 2-day recall window
* to be consistent with the 7-day recall window
*
* assume for 7d prevalence, that if the child had
* prevalent symptoms in the past 2 days that they had symptoms
* within the past 7 days (there are a few inconsistent responses)
*
* likewise, assume that for 2d prevalence if
* the child a measurement for the 7d recall
* but had missing info for the past 2d
* then impute with the 7d recall value
*--------------------------------------------

gen byte d3plus2d = (q203a==1) | (q203b==1) | (q203c==1)
	replace d3plus2d = . if (q203a==.) & (q203b==.) & (q203c==.)
gen byte d3plus7d = (q203d==1)
	replace d3plus7d = . if (q203d==.)
	label var d3plus2d "3+ stools in 24 hr, 2d recall"
	label var d3plus7d "3+ stools in 24 hr, 7d recall"

gen byte dloose2d = (q205a==1) | (q205b==1) | (q205c==1)
	replace dloose2d = . if (q205a==.) & (q205b==.) & (q205c==.)
gen byte dloose7d =  (q205d==1)
	replace dloose7d = . if (q205d==.) 
	label var dloose2d "loose or watery stool, 2d recall"
	label var dloose7d "loose or watery stool, 7d recall"

gen byte dblood2d = (q206a==1) | (q206b==1) | (q206c==1)
	replace dblood2d = . if (q206a==.) & (q206b==.) & (q206c==.)
gen byte dblood7d =  (q206d==1)
	replace dblood7d = . if (q206d==.) 
	label var dblood2d "blood in stool, 2d recall"
	label var dblood7d "blood in stool, 7d recall"

gen byte diar2d = (d3plus2d==1 & dloose2d==1) | (dblood2d==1)
	label var diar2d "Diarrhea case, 2d recall"
	notes diar2d: Defined as 3+ stools in 24 hrs or 1+ stool w/ blood; includes interview day
	replace diar2d = . if (d3plus2d==.) & (dloose2d==.) & (dblood2d==.)
	replace diar2d = . if (d3plus2d==. | dloose2d==.) & (dblood2d==.)	
	
	
gen byte diar7d = (d3plus7d==1 & dloose7d==1) | (dblood7d==1)
	label var diar7d "Diarrhea case, 7d recall"
	notes diar7d: Defined as 3+ stools in 24 hrs or 1+ stool w/ blood; includes interview day
	replace diar7d = 1 if (diar2d==1) & (diar7d==0 | diar7d==.)
	replace diar7d = . if (d3plus7d==.) & (dloose7d==.) & (dblood7d==.) & (diar2d==.)
	replace diar7d = . if (d3plus7d==. | dloose7d==.) & (dblood7d==.) & (diar2d==.)	
	
	*drop d3plus* dloose* dblood*

		
tab1 diar2d diar7d, m
tab diar2d diar7d, m


bysort tchild: tab svy diar7d, row

*--------------------------------------------
* do some diagnostics on children missing
* diarrhea outcome measurements
*--------------------------------------------

* by sibling status
tab tchild diar7d, m

* by randomization block
tab block diar7d, m

* by tr status (useful after unblinding)
tab tr diar7d, m

*--------------------------------------------
* create negative control symptom variables
*--------------------------------------------
gen byte bruise2d = (q211a==1) | (q211b==1) | (q211c==1)
	replace bruise2d = . if (q211a==.) & (q211b==.) & (q211c==.)
	label var bruise2d "Bruising, 2d recall"
gen byte bruise7d = (q211d==1)
	replace bruise7d = . if (q211d==.)
	replace bruise7d = 1 if (bruise2d==1) & (bruise7d==0 | bruise7d==.)
	label var bruise7d "Bruising, 7d recall"
gen byte tooth2d = (q212a==1) | (q212b==1) | (q212c==1)
	replace tooth2d = . if (q212a==.) & (q212b==.) & (q212c==.)
	label var tooth2d "Toothache, 2d recall"
gen byte tooth7d = (q212d==1)
	replace tooth7d = . if (q212d==.)
	replace tooth7d = 1 if (tooth2d==1) & (tooth7d==0 | tooth7d==.)
	label var tooth7d "Toothache, 7d recall"
	
*--------------------------------------------
* Save a diarrhea dataset
*--------------------------------------------

label var dataid "Unique compound ID"
label var hhid "Household ID"
label var childid "Child ID"
label var clusterid "Cluster ID"
label var block "Randomization block ID"


label define sex 0 "female" 1 "male"
label values sex sex
label var sex "Sex (1=male)"


* restrict to identifying vars and child level variables
keep dataid childid tchild clusterid block tr svy svydate month sex dob agedays ageyrs enrolage newbirth sibnewbirth gt36mos  d3plus* dloose* dblood* diar2d diar7d bruise* tooth*
order  dataid childid tchild clusterid block tr svy svydate month sex dob agedays ageyrs enrolage newbirth sibnewbirth gt36mos d3plus* dloose* dblood* diar2d diar7d bruise* tooth*

* drop households at enrollment with no target children and no siblings
assert svy==0 if childid==""
keep if childid!=""
compress
sort dataid childid 
label data "Bangladesh diarrhea dataset, created by 3-bangladesh-dm-diar.do"
saveold "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-diar.dta", replace version(12)
outsheet using "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-diar.csv", comma replace


* write a codebook for the dataset
log close
log using "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-diar-codebook.txt", text replace
desc
codebook, c
codebook
log close


exit

