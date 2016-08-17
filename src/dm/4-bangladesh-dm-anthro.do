capture log close
set more off
clear all


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/4-bangladesh-dm-anthro.log", text replace

*--------------------------------------------
* 4-bangladesh-dm-anthro.do
*
* ben arnold (benarnold@berkeley.edu)
*
* process the Bangladesh trial year 1 and 
* year 2 anthropometry data for analysis
*
*--------------------------------------------




*--------------------------------------------
* input files:
*
*  Treatment assignments 
*	washb-bang-tr.dta (or washb-bang-blind-tr.dta if blinded)
* 
*
*  ENROLLMENT
*  1. WASHB_Baseline_main_survey.dta
*  3. WASHB_Baseline_childinfo.dta
*
*  YEAR 1
*  1. WASHB_Midline_main_survey_cleaned.dta
*  3. WASHB_Midline_childinfo_cleaned.dta
*  8. WASHB_Midline_anthropometry_cleaned.dta
*
*  YEAR 2
*  04. WASHB_Endline_main_survey_cleaned.dta
*  07. WASHB_Endline_childinformation_cleaned.dta
*  09. WASHB_Endline_anthropometry_cleaned.dta
*
* output files:
*  final/washb-bangladesh-anthro.dta / .csv
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

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/2_Midline/8. WASHB_Midline_anthropometry_cleaned.dta", clear
duplicates list dataid childid
sort dataid childid
tempfile childanthro
save `childanthro'


use `main', clear
merge 1:m dataid using `childinfo'
assert _merge!=1
drop if _merge==2
drop _merge

* merge in the child anthro information
sort dataid childid
merge 1:1 dataid childid using `childanthro'


* list 3 children who are not matching
* dropped, after confirming with Kishor
count if _merge==2
list dataid clusterid bariid motherid childid  if _merge==2
drop if _merge==2

keep if _merge==3
drop _merge

* assign survey
gen svy = 1 
	label var svy "Survey round (0,1,2)"
	
sort dataid childid
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

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/09. WASHB_Endline_anthropometry_cleaned.dta", clear


* for now, drop sibling measurements because there is no clean way to merge them
* to the childinfo data
tab childid
drop if childid=="S1"

* rename variables to be consistent with year 1
* sub the "an" prefix for a "c" prefix
foreach var of varlist an* {
	local vname = "`var'"
	local vstub = substr("`vname'",3,.)
	rename `var' c`vstub'
}


* list duplicates
duplicates list dataid childid
sort dataid childid
tempfile childanthro
save `childanthro'

use `main', clear
merge 1:m dataid using `childinfo'
assert _merge==3
drop _merge

* merge in the child anthro information
sort dataid childid
merge 1:1 dataid childid using `childanthro'


* list 4 children who are not matching
* dropped, after confirming with Kishor
count if _merge==2
list dataid childid  if _merge==2
drop if _merge==2
keep if _merge!=1
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
* APPEND YEAR 1, YEAR 2
*--------------------------------------------
gen svy = 2
	label var svy "Survey round (0,1,2)"
	
append using `year1'



*--------------------------------------------
* merge in the treatment assignment information 
*--------------------------------------------

* drop clusterid and block variables to 
* ensure we get a complete set

* merge in the treatment assignment info (keep only matching obs)
capture drop clusterid
gen clusterid = substr(dataid,1,3)
sort clusterid
merge m:1 clusterid using `trdata'
assert _merge != 1
drop if _merge==2
drop _merge


*--------------------------------------------
* identify target children
*--------------------------------------------
gen byte tchild = strpos(childid,"T")>0
	label var tchild "Target child in birth cohort"
	label define tchild 0 "Sibling" 1 "Target child"
	label values tchild tchild
	tab tchild

	
* check for duplicates (there should be none)
duplicates report dataid childid svy

tab svy tchild


*--------------------------------------------
* for a small number of children without age
* or sex information in one of the measurement
* rounds, fill it in using the other round
*--------------------------------------------


*--------------------------------------------
* Age, using measurement date and birth date
*--------------------------------------------
label var sex "Sex"
label define sex 0 "female" 1 "male"
label values sex sex

label var dob "Date of birth"
format dob %d

* date was saved in a different variable in Y2
replace c4002 = date(c01,"DMY") if svy==2

gen anthrodate = c4002
	format anthrodate %d
	label var anthrodate "Date of anthro measurement"
	
	* for 4 missing anthro dates at midline, impute with
	* either the child info date (only available at midline)
	* or the household survey date
	replace anthrodate=entrydt if anthrodate==.
	assert anthrodate!=.	
	codebook anthrodate

gen aged = anthrodate-dob
	label var aged "Age in days (anthro meas)"
gen double agem = aged/30.4167
	label var agem "Age in months (anthro meas)"
gen double agey = aged/365.25
	label var agey "Age in years (anthro meas)"
codebook agey

* Month of measurement
gen month = month(anthrodate)
	label var month "Month of measurement"


* check children who seem too old
list dataid childid motherid dob aged agem if agem>15 & agem<. & svy==1
list dataid childid motherid dob aged agem if agem>32 & agem<. & svy==2


*--------------------------------------------
* birth order of target children
* collected at year 2, so backfill
*--------------------------------------------
gen int birthord = .
	replace birthord = q4020_1b if childid=="T1"
	replace birthord = q4021_1b if childid=="T2"
	label var birthord "Birth order (target child)"
	
sort dataid childid svy
by dataid childid: egen _x = min(birthord)
replace birthord = _x if birthord==.
drop _x

* some missing values in the measure
* 243 at year 1
tab birthord if svy==1, m
tab birthord if svy==2, m

*--------------------------------------------
* ensure that all the anthro measurements 
* are rounded to the correct level of precision
* no more than 2 decimal places
*--------------------------------------------
for any c414 c415 c416 c404 c405 c406 c408 c409 c410 c411 c412 c413 c418 c419 c420: replace X = round(X,0.01)

*--------------------------------------------
* Calculate median length measurements
*--------------------------------------------
label var c414 "Child length meas 1"
label var c415 "Child length meas 2"
label var c416 "Child length meas 3"

gen len1 = c414
gen len2 = c415
gen len3 = c416
replace len1 = . if len1>999 | len1<=0
replace len2 = . if len2>999 | len2<=0
replace len3 = . if len3>999 | len3<=0


egen float length = rowmedian(len1 len2 len3)
	replace length = . if length>999
	replace length = round(length,0.001)
	label var length "Child length (median), cm"
	notes length: Median of replicate measures

drop len1-len3

*--------------------------------------------
* Calculate median weight measurements
*--------------------------------------------

label var c404 "Maternal weight meas 1"
label var c405 "Maternal weight meas 2"
label var c406 "Maternal weight meas 3"
label var c408 "Child weight w/ mother meas 1"
label var c409 "Child weight w/ mother meas 2"
label var c410 "Child weight w/ mother meas 3"
label var c411 "Child weight meas 1"
label var c412 "Child weight meas 2"
label var c413 "Child weight meas 3"
for any c404 c405 c406 c408 c409 c410 c411 c412 c413 : replace X = . if (X>99 | X<=0)


* create median of mom + child and mom alone
* we have to do it this way because this reflects
* the method used in the field to determine if a 
* 3rd measurement was required
egen float wmc = rowmedian(c408 c409 c410)
egen float wmm = rowmedian(c404 c405 c406)

* create median of the child (if measured alone)
egen float wch = rowmedian(c411 c412 c413)

gen float weight = wch
	replace weight = wmc-wmm if (weight==.) & (wmc!=.) & (wmm!=.)
	replace weight = round(weight,0.001)
	label var weight "Child weight (median), kg"
	notes weight: Median of replicate measures

drop  wmc wmm wch


*--------------------------------------------
* Calculate median head circumference measurements
*--------------------------------------------
label var c418 "Head circumference meas 1"
label var c419 "Head circumference meas 2"
label var c420 "Head circumference meas 3"

gen hc1 = c418
gen hc2 = c419
gen hc3 = c420
replace hc1 = . if hc1>99 | hc1<=0
replace hc2 = . if hc2>99 | hc2<=0
replace hc3 = . if hc3>99 | hc3<=0

egen float headcir = rowmedian(hc1 hc2 hc3)
	replace headcir = . if headcir>99
	replace headcir = round(headcir,0.001)
	label var headcir "Child head circumference (median), cm"
	notes headcir: Median of replicate measures


*--------------------------------------------
* Calculate laz, waz, whz, using the zscore06
* add-on package. do not specify oedema
*--------------------------------------------

zscore06, a(agem) s(sex) h(length) w(weight) female(0) male(1) measure(c417) recum(1) stand(2)

* save a temporary dataset to merge with the WHO igrowup output
save "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-working-anthro.dta", replace


*--------------------------------------------
* Calculate laz, waz, whz, and wcz using the 
* WHO -igrowup- macro. 
*
* this is necessary because the zscore06 
* package does not calculate z-scores for 
* head circumference.  we are using both
* approaches because the zscore06 package is
* ostensibly more accurate (see the package's documentation)
* though as demonstrated below they provide identical results
* (just a good internal validity check)
*
* The WHO -igrowup- macro requires 15 parameters.
* Refer to the documentation for details on this
* finicky macro.
*
* We don't have all of the measurements that
* it processes, so need to create empty
* variables that allow it to run
*---------------------------------------------

/// Indicate to the Stata compiler where the igrowup_standard.ado file is stored ***
adopath + "~/WASHB-Bangladesh-primary-outcomes/src/dm/4-who-igrowup"


gen str reflib="~/WASHB-Bangladesh-primary-outcomes/src/dm/4-who-igrowup"
lab var reflib "Directory of reference tables"

gen str datalib="~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets"
lab var datalib "Directory for datafiles"

gen str datalab="TEMPanthro"
lab var datalab "Working file"

gen str ageunit="days"
gen str measure="L" if c417==1
replace measure="H" if c417==2
replace measure=" " if measure == ""
label var measure "Height measured lying -L- or standing -H-"

* create a temporary sex variable for WHO coding
gen whosex = sex
	replace whosex = 2 if sex==0

* create missing variables so that macro will run
for any uac triskin subskin oedema: gen X = .

* set sampling wgtghts to negative to make the prevalence
* calculations blow up -- impossible to run that piece of 
* code w/o Stata SE b/c requires 10,000 variables
gen sw = -10


*---------------------------------------------
* Save a temporary dataset and run -igrowup-
* (note: igrowup adds variables to the data in
* memory and saves a copy dataset with the 
* suffix "_z_st.dta"). Dataset name must correspond
* to the "datalab" variable (defined in last chunk)
*---------------------------------------------

keep dataid childid svy aged whosex whosex aged ageunit weight length measure headcir uac triskin subskin oedema sw  reflib datalib datalab 

save "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/TEMPanthro", replace

#delimit;
igrowup_standard reflib datalib datalab 
whosex aged ageunit weight length measure headcir uac triskin subskin oedema sw;
#delimit cr



*---------------------------------------------
* Retrieve WHO calculated output
* "_f" variables identify variables outside of
* reasonable bounds
*
* merge back to the main anthro dataset
*---------------------------------------------

use "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/TEMPanthro_z_st", clear

keep dataid childid svy _zwei _zlen _zbmi _zwfl _zhc _fwei _flen _fbmi _fwfl _fhc
sort dataid childid svy
save  "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/TEMPanthro", replace

use "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-working-anthro.dta", clear
sort dataid childid svy
merge 1:1 dataid childid svy using "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/TEMPanthro"
assert _merge==3
drop _merge


* compare measurements from the 2 packages to ensure they are identical
* (correlation = 1)
corr haz06 _zlen 
corr waz06 _zwei
corr whz06 _zwfl
corr bmiz06 _zbmi


* delete tempfiles
erase "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-working-anthro.dta"
erase "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/TEMPanthro.dta"
erase "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/TEMPanthro_z_st.dta"
erase "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/TEMPanthro_z_st.xls"

*--------------------------------------------
* Set extreme values to missing and flag them
* based on the WHO 2006 standards
*--------------------------------------------
rename haz06 laz
rename waz06 waz
rename whz06 whz
rename bmiz06 bmiz
rename _zhc hcz

gen laz_x = (laz < -6 | laz >6)
	replace laz_x = . if laz==.
	label var laz_x "abs(LAZ)>6, set to missing"
	
gen waz_x = (waz < -6 | waz >5)
	replace waz_x = . if waz==.
	label var waz_x "WAZ < -6 or WAZ > 5, set to missing"

gen whz_x = (whz < -5 | whz >5)
	replace whz_x = . if whz==.
	label var whz_x "abs(WHZ)>5, set to missing"
	
gen bmiz_x = (bmiz < -5 | bmiz >5)
	replace bmiz_x = . if bmiz==.
	label var bmiz_x "abs(BMIZ)>5, set to missing"

gen hcz_x = _fhc
	replace hcz_x = . if hcz==.
	label var hcz_x "abs(HCZ)>5, set to missing"
	
* list extreme values before setting them to missing
list dataid aged length weight laz waz whz if laz_x==1
list dataid aged length weight laz waz whz if waz_x==1
list dataid aged length weight laz waz whz if whz_x==1
list dataid aged length weight laz waz whz if bmiz_x==1
list dataid aged headcir hcz if hcz_x==1

replace laz = . if laz_x==1
replace waz = . if waz_x==1
replace whz = . if whz_x==1
replace bmiz = . if bmiz_x==1
replace hcz = . if hcz_x==1

* drop extra Z-score calculation variables
drop _z* _f*

*--------------------------------------------
* Identify children who are stunted, 
* underweight, or wasted based on their Z-scores
*--------------------------------------------

gen byte lazminus2 = laz < -2
	replace lazminus2 =. if laz==. | laz_x==1
	label var lazminus2 "Stunted (LAZ<-2)"
gen byte lazminus3 = laz < -3
	replace lazminus3 =. if laz==. | laz_x==1
	label var lazminus3 "Severely Stunted (LAZ<-3)"

gen byte wazminus2 = waz < -2
	replace wazminus2 =. if waz==. | waz_x==1
	label var wazminus2 "Underweight (WAZ<-2)"
gen byte wazminus3 = waz < -3
	replace wazminus3 =. if waz==. | waz_x==1
	label var wazminus3 "Severely Underweight (WAZ<-3)"

gen byte whzminus2 = whz < -2
	replace whzminus2 =. if whz==. | whz_x==1
	label var whzminus2 "Wasted (WHZ<-2)"
gen byte whzminus3 = whz < -3
	replace whzminus3 =. if whz==. | whz_x==1
	label var whzminus3 "Severely Wasted (WHZ<-3)"
	
	
*--------------------------------------------
* Save an analysis dataset
*--------------------------------------------

label var dataid "Compound ID"
label var childid "Child ID"
label var clusterid "Cluster ID"
label var block "Randomization block ID"

* restrict to variables used in the analysis
keep dataid childid motherid tchild clusterid block tr svy anthrodate month c404-c406 c408-c410 c411-c413 c414-c416 c418-c420 sex birthord aged agem agey weight length headcir laz* waz* whz* bmiz* hcz* 
order dataid childid motherid tchild clusterid block tr svy anthrodate month c404-c406 c408-c410 c411-c413 c414-c416 c418-c420 sex birthord aged agem agey weight length headcir laz* waz* whz* bmiz* hcz*

compress
sort dataid childid svy
label data "Bangladesh anthropometry analysis dataset, created by 4-bangladesh-dm-anthro.do"
save "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-anthro.dta", replace
outsheet using "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-anthro.csv", comma replace

* write a codebook for the dataset
log close
log using "~/dropbox/washb-bangladesh-data/1-primary-outcome-datasets/washb-bangladesh-anthro-codebook.txt", text replace
desc
codebook, c
codebook
log close

exit

