capture log close
set more off
clear all


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/2-bangladesh-dm-enrol.log", text replace

*--------------------------------------------
* 2-bangladesh-dm-enrol.do
*
* ben arnold (benarnold@berkeley.edu)
*
* process the Bangladesh trial enrollment dataset
* 
* this includes enrollment adjustement covariates
* along with uptake variable measures
*
* it creates 2 files: a household-level
* enrollment dataset with analysis covariates
* and a child-level enrollment dataset with
* diarrhea outcomes
*
*--------------------------------------------




*--------------------------------------------
* input files:
*
*  washb-bang-tr.dta (or washb-bang-blind-tr.dta if blinded still)
*
*  1. WASHB_Baseline_main_survey.dta
*  2. WASHB_Baseline_census.dta
*
* output files:
*  final/washb-bangladesh-enrol.dta / .csv
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
* calculate a few compound-level variables
* from the census for the merge
*--------------------------------------------


use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/1_Baseline/2. WASHB_Baseline_census.dta", clear
sort dataid hhid
bysort dataid: egen Ncomp = sum(a6)
	label var Ncomp "Total indivs in compound"
gen Nlt18 = a4+a5
	label var Nlt18 "N indivs in HH <=18y"
replace hhid = subinstr(hhid,"H","",.)
keep dataid hhid Ncomp Nlt18
sort dataid hhid
tempfile bcensus
save `bcensus'

*--------------------------------------------
* merge the household and census data
*--------------------------------------------

use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/1_Baseline/1. WASHB_Baseline_main_survey.dta", clear
sort dataid 
tempfile bmain
save `bmain'

* merge census information to the data
gen hhid = q4016
sort dataid hhid

merge 1:1 dataid hhid using `bcensus'
assert _merge!=1
* drop households in same compounds that were not actually surveyed 
keep if _merge==3
drop _merge


* merge in the treatment assignment info (keep only matching obs)
gen clusterid = substr(dataid,1,3)
sort clusterid
capture drop block tr
merge m:1 clusterid using `trdata'
assert _merge == 3
drop _merge

order dataid clusterid block tr
sort dataid


* format survey dates
gen svydate = q4002
	format svydate %d
	label var svydate "Survey date"
	codebook svydate

*--------------------------------------------
* Adjustment covariates
*--------------------------------------------

* Administrative Union 
rename q4008 union
	label var union "Administrative union (q4008)"
	
* FRA ID code
rename q4001 fracode 
	label var fracode "FRA code (q4001)"
	
* Mother's age
gen momage = q103 if (q102==1)
	label var momage "Mother's age, years"

* Mother's education
gen momedu = .
	replace momedu = 0 if (q009==0)
	replace momedu = 1 if (q009>=1 & q009<=5)
	replace momedu = 2 if (q009>5 & q009<.)
	label define momedu 0 "No education" 1 "Primary (1-5y)" 2 "Secondary (>5y)"
	label values momedu momedu
	label var momedu "Mother's education (category)"

* Household food insecurity
* See definitions in sections 5.3 and 5.4 of
* Household Food Insecurity Access Scale (HFIAS) for Measurement of Food Access: Indicator Guide
* Version 3 (Coates et al. 2007) available at:
* http://www.fao.org/fileadmin/user_upload/eufao-fsi4dm/doc-training/hfias.pdf
egen int hfias = rowtotal(q110*a)
	label var hfias "HFIAs food insecurity score (0-27)"
gen int hfiacat = .
	replace hfiacat = 1 if (q1101==0 | q1101a==1) & (q1102==0 & q1103==0 & q1104==0 & q1105==0 & q1106==0 & q1107==0 & q1108==0 & q1109==0)
	replace hfiacat = 2 if ((q1101a==2 | q1101a==3 | q1102a==1 | q1102a==2 | q1102a==3 | q1103a==1 | q1104a==1) & q1105==0 & q1106==0 & q1107==0 & q1108==0 & q1109==0)
	replace hfiacat = 3 if ((q1103a==2 | q1103a==3 | q1104a==2 | q1104a==3 | q1105a==1 | q1105a==2 | q1106a==1 | q1106a==2) & q1107==0 & q1108==0 & q1109==0)
	replace hfiacat = 4 if (q1105a==3 | q1106a==3 | q1107a==1 | q1107a==2 | q1107a==3 | q1108a==1 | q1108a==2 | q1108a==3 | q1109a==1 | q1109a==2 | q1109a==3)
	label define hfiacat 1 "Food Secure" 2 "Mildly Food Insecure" 3 "Moderately Food Insecure" 4 "Severely Food Insecure"	
	label values hfiacat hfiacat
	label var hfiacat "HFIAS food insecurity category"
	
* Distance to household's primary drinking water source (minutes)
** Check dataid=="26902", with time to water of 10 hours???
gen watmin = real(q1011hours)*60 + real(q1011mins)
	label var watmin "Dist (mins) to primary water source"
	
* Main roof material of house
gen byte roof = inlist(q4101,2,3)
	replace roof = . if q4101==.
	label var roof "Improved roof material (tin, cement)"
	
* Main wall material of house
gen byte walls = inlist(q4102,2,3,4)
	replace walls = . if q4102==.
	label var walls "Improved wall material (wood, brick, tin)"

* Main floor material of house
gen byte floor = inlist(q4103,2,3)
	replace floor = . if q4103==.
	label var floor "Improved floor (wood, concrete)"

* Electricity
rename q4105_a elec
	label var elec "Household has electricity"

* Household Assests
rename q4105_b n_asset_wardrobe
rename q4105_c n_asset_table
rename q4105_d n_asset_chair
rename q4105_e n_asset_clock
rename q4105_f n_asset_khat
rename q4105_g n_asset_chouki
rename q4105_h asset_radio
rename q4105_i asset_tvbw
	replace asset_tvbw = 0 if asset_tvbw==9
rename q4105_j asset_tvcol
rename q4105_k asset_refrig
rename q4105_l asset_bike
rename q4105_m asset_moto
rename q4105_n asset_sewmach
rename q4105_o n_asset_mobile
rename q4105_p asset_phone

* reduce some household assets to binary variables
gen asset_tv = (asset_tvbw==1) | (asset_tvcol==1)
	replace asset_tv=. if (asset_tvbw==.) & (asset_tvcol==.)
	label var asset_tv "Has bw or color TV"

local assets "wardrobe table chair clock khat chouki mobile"
foreach var of local assets {
	gen asset_`var' = n_asset_`var'>0
		replace asset_`var' = . if n_asset_`var' == .
}
label var asset_wardrobe "Has >=1 wardrobe"
label var asset_table "Has >=1 table"
label var asset_chair "Has >=1 chair"
label var asset_clock "Has >=1 clock"
label var asset_khat "Has >=1 khat"
label var asset_chouki "Has >=1 chouki"
label var asset_mobile "Has >=1 mobile"


* person counts
* Number of children < 18 y in the household (Nlt18 - computed above w/ census)
* Number of individuals living in the compound (Ncomp - computed above w/ census)
	
/* REMAINING PRE-SPECIFIED COVARIATES
Child birth order  (in child-level datasets -- only for index children)
Mother’s height (cm) -- added below

*/

*--------------------------------------------
* additional balance table characteristics
*--------------------------------------------

** maternal
*  age (completed above)
*  parity (not collected at enrollment)

*  years of education q009
gen momeduy = q009
	replace momeduy = . if inlist(q009,.,999)
	label var momeduy "Mother, years of education (q009)"
	
** paternal
*  age (not collected at enrollment)
*  years of education q010
gen dadeduy = q010
	replace dadeduy = . if inlist(q010,.,999)
	label var dadeduy "Father, years of education (q010)"
	
*  works in agriculture (%) q011
gen dadagri = inlist(q011,1,3)
	replace dadagri = . if inlist(q011,.,35,999)
	label var dadagri "Father works in agriculture (q011)"

** household 
*  number of persons q012
gen Nhh = q012
	label var Nhh "Number of persons in household (q012)"

*  has electricity (above)

*  has cement floor
gen byte cement = inlist(q4103,3)
	replace cement = . if q4103==.
	label var cement "Has cement floor (q4103)"

*  acres of agricultural land owned
* decimal to acre conversion: 1/100
* decimal to hectares conversion: 0.004046
* include q4112??
gen landacre = q4110/100
	replace landacre = . if inlist(q4110,.,999)
	label var landacre "Acres of land owned (q4110)"

** drinking water
*  primary water source = shallow tubewell
gen byte tubewell = (q1010==1)
	replace tubewell = . if q1010==.
	label var tubewell "Primary water source: shallow tubewell (q1010)"
	
*  stored water observed at home
recode q1003_6 q1003_7 q1003_8 (.=0)
gen byte storewat = (q1003_6==1 | q1003_7==1 | q1003_8==1)
	replace storewat = . if (q1003_6==. & q1003_7==. & q1003_8==.)
	label var storewat "Store drinking water (q1003)"

*  reported treating water yesterday
gen byte treatwat = (q1007==1) & inlist(q1008,1,2)
	replace treatwat = . if (q1007==.) & inlist(q1008,999)
	label var treatwat "Reported treating water today or yesterday (q1008)"

	
** Sanitation
*  daily OD
*  adult men
gen byte odmen = (q801_a==1)
	replace odmen = . if inlist(q801_a,.,888,999)
	label var odmen "Daily OD, adult men (q801a)"
*  adult women
gen byte odwom = (q801_b==1)
	replace odwom = . if inlist(q801_b,.,888,999)
	label var odwom "Daily OD, adult women (q801b)"
*  children 8-<15 years
gen byte odch815 = (q801_e==1)
	replace odch815 = . if inlist(q801_e,.,888,999)
	label var odch815 "Daily OD, children 8-15 (q801e)"
*  children 3-<8 years
gen byte odch38 = (q801_d==1)
	replace odch38 = . if inlist(q801_d,.,888,999)
	label var odch38 "Daily OD, children 3-8 (q801d)"
*  children <3 years
gen byte odchu3 = (q801_c==1)
	replace odchu3 = . if inlist(q801_c,.,888,999)
	label var odchu3 "Daily OD, children <3 (q801c)"

** latrine
*  owned (%)
gen byte latown = (q816==1)
	replace latown = . if (q808==1) & inlist(q816,.,888)
	label var latown "Own their latrine (not shared) (q816)"

*  concrete slab (%)
gen byte latslab = (q809_7==1)
	replace latslab = . if inlist(q809_7,.,888)
	label var latslab "Latrine has concrete slab (q809_7)"

*  functional water seal
gen byte latseal = (q809_9a==1)
	replace latseal = . if  inlist(q809_9a,.,888)
	label var latseal "Latrine has functional water seal (q809_9a)"
	
*  no visible stool on slab or floor
gen byte latfeces = (q809_16!=1)
	replace latfeces = . if inlist(q809_16,.,888)
	label var latfeces "No visible feces on slab/floor of latrine"

* owned a potty
gen byte potty = (q913==1)
	replace potty = . if inlist(q913,.,999)
	label var potty "Has a potty for child defecation"

* human feces observed in the household compound
gen byte humfeces = 0
	replace humfeces = 1 if (q4201>0) | (q4203>0) | (q4205>0)
	replace humfeces = . if inlist(q4201,99,.) & inlist(q4203,99,.) & inlist(q4205,99,.)
	label var humfeces "Human feces observed in house/compound (q4201,4203,4205)"
	
* human feces observed in the location where child spends the most time
gen byte humfecesch = 0
	replace humfecesch = 1 if (q4203>0)
	replace humfecesch = . if inlist(q4203,99,.)
	label var humfecesch "Human feces obs in loc where child spends most time (q4203)"

** handwashing
*  within 6 steps of latrine
gen byte hwlat = (q703==2 & q707<=6) | (q709==2 & q713<=6)
	replace hwlat = . if (q703==.) & (q709==.)
	label var hwlat "Has handwashing loc w/in 6 steps of latrine"
*   has water
gen byte hwlatwat = (q703==2 & q707<=6 & q704_1==1) | (q709==2 & q713<=6 & q710_1==1)
	replace hwlatwat = . if (q704_1==.) & (q710_1==.)
	label var hwlatwat "Has handwashing loc w/in 6 steps of latrine w/ water"
*   has soap
gen byte hwlatsoap = (q703==2 & (q704_2==1 | q704_3==1 | q704_4==1 | q704_5==1 | q704_6==1) ) | (q709==2 & (q710_2==1 | q710_3==1 | q710_4==1 | q710_5==1 | q710_6==1))
	replace hwlatsoap = . if (q704_2==.) & (q704_3==.) & (q704_4==.) & (q704_5==.) & (q704_6==.) & (q710_2==.) & (q710_3==.) & (q710_4==.) & (q710_5==.) & (q710_6==.)
	label var hwlatsoap "Has handwashing loc w/in 6 steps of latrine w/ soap"

*  within 6 steps of kitchen
gen byte hwkit = (q703==3 & q706<=6) | (q709==3 & q712<=6)
	replace hwkit = . if (q703==.) & (q709==.)
	label var hwkit "Has handwashing loc w/in 6 steps of kitchen"
*   has water
gen byte hwkitwat = (q703==3 & q706<=6 & q704_1==1) | (q709==3 & q712<=6 & q710_1==1)
	replace hwkitwat = . if (q704_1==.)  & (q710_1==.)
	label var hwkitwat "Has handwashing loc w/in 6 steps of kitchen w/ water"
*   has soap
gen byte hwkitsoap = (q703==3 & q706<=6 & (q704_2==1 | q704_3==1 | q704_4==1 | q704_5==1 | q704_6==1) ) | (q709==3 & q712<=6 & (q710_2==1 | q710_3==1 | q710_4==1 | q710_5==1 | q710_6==1))
	replace hwkitsoap = . if  (q704_2==.) & (q704_3==.) & (q704_4==.) & (q704_5==.) & (q704_6==.) & (q710_2==.) & (q710_3==.) & (q710_4==.) & (q710_5==.) & (q710_6==.)
	label var hwkitsoap "Has handwashing loc w/in 6 steps of kitchen w/ soap"


*--------------------------------------------
* Updtake indicators
*--------------------------------------------

* has stored drinking water (above)

	
* store water with detectable free chlorine (not at baseline)
*gen byte freechl = q1027level>0 & q1027level<.
*	replace freechl = . if (svy==0) | (q1027==999)
*	label var freechl "Free chlorine detected in stored water (>0 mg/L)"
	
* latrine with a functional water seal (above)

* no visible feces on the latrine slab or floor (above)

* no human feces in house or compound (above)

	
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

* mean saches of LNS fed in prior week to index child 6-24 mos (not at baseline)


*--------------------------------------------
* Save an enrollment analysis dataset
* household level observations
*--------------------------------------------

label var hhid "Household ID"
label var clusterid "Cluster ID"
label var block "Randomization block ID"

* restrict to household level variables used in the analysis
* and save the data
keep dataid clusterid hhid block tr union fracode svydate Nhh Nlt18 Ncomp momage momedu momeduy dadeduy	dadagri landacre hfias hfiacat tubewell watmin storewat treatwat odmen odwom odch815 odch38 odchu3 latown latslab latseal latfeces potty humfeces humfecesch hwlat* hwkit* hwsw hwss hwsws roof walls floor cement elec asset_* n_asset_*

order dataid clusterid hhid block tr union fracode svydate Nhh Nlt18 Ncomp momage momedu momeduy dadeduy	dadagri landacre hfias hfiacat tubewell watmin storewat treatwat odmen odwom odch815 odch38 odchu3 latown latslab latseal latfeces potty humfeces humfecesch hwlat* hwkit* hwsw hwss hwsws roof walls floor cement elec asset_* n_asset_*

** double check to ensure there are no duplicates
duplicates drop dataid clusterid hhid block tr union fracode svydate Nhh Nlt18 Ncomp momage momedu momeduy dadeduy	dadagri landacre hfias hfiacat tubewell watmin storewat treatwat odmen odwom odch815 odch38 odchu3 latown latslab latseal latfeces potty humfeces humfecesch hwlat* hwkit* hwsw hwss hwsws roof walls floor cement elec asset_* n_asset_*, force

compress
sort dataid 
label data "Bangladesh enrollment analysis dataset (HH obs), created by 2-bangladesh-dm-enrol.do"
saveold "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-enrol.dta", replace version(12)

*--------------------------------------------
* There is a single variable measured at
* other timepoints that we need to merge into
* the enrollment data: maternal height
* get maternal height measurements from
* year 1 and year 2 visits, then merge it
* into the enrollment data
*--------------------------------------------

*** Year 1 data
use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/2_Midline/8. WASHB_Midline_anthropometry_cleaned.dta", clear

* calculate median maternal height
for any c422 c423 c424: replace X = . if X >999
egen float momheight1 = rowmedian(c422 c423 c424)
	replace momheight1 = . if momheight >999
	label var momheight1 "Maternal height (median), year 1 visit"
	
* restrict to measurements without missing maternal height
keep if momheight1 != .

* bring along individual measurements
gen y1hgt1 = c422
gen y1hgt2 = c423
gen y1hgt3 = c424

* restrict to single observation per compound
* (confirmed that maternal height is duplicated in the data
*  or is only reported for one child in the 31 cases where 
*  there are multiple target children in the home)
keep dataid momheight1 y1hgt*
bysort dataid: keep if _n==1

sort dataid
tempfile momhgt1
save `momhgt1'

*** Year 2 data
use "~/dropbox/WASHB-Bangladesh-Data/0-Untouched-data/1-Main-survey/3_Endline/09. WASHB_Endline_anthropometry_cleaned.dta", clear


* calculate median maternal height
for any an422 an423 an424: replace X = . if X >999
egen float momheight2 = rowmedian(an422 an423 an424)
	replace momheight2 = . if momheight >999
	label var momheight2 "Maternal height (median), year 2 visit"
	
* restrict to measurements without missing maternal height
keep if momheight2 != .

* bring along individual measurements
gen y2hgt1 = an422
gen y2hgt2 = an423
gen y2hgt3 = an424

* restrict to single observation per compound
* (confirmed that maternal height is duplicated in the data
*  or is only reported for one child in the cases where 
*  there are multiple target children in the home)
keep dataid momheight2 y2hgt*
bysort dataid: keep if _n==1

* merge the year 1 and year 2 data
sort dataid
merge 1:1 dataid using `momhgt1'


* calculate the overall median
egen float momheight = rowmedian(y1hgt1 y1hgt2 y1hgt3 y2hgt1 y2hgt2 y2hgt3)
	label var momheight "Maternal height (median)"

* replace the height measurement with round 1 median if the two medians differ by >3 cm
gen diff = momheight1-momheight2
replace momheight = momheight1 if abs(diff)>3 & (diff!=.)

* merge in the maternal height variable
keep dataid momheight
sort dataid
tempfile momheight
save `momheight'

use "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-enrol.dta", clear
merge 1:1 dataid using `momheight'
assert _merge !=2
drop _merge

order dataid clusterid hhid block tr union fracode svydate Nhh Nlt18 Ncomp momage momheight

saveold "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-enrol.dta", replace version(12)
outsheet using "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-enrol.csv", comma replace
desc
codebook, c

* write a codebook for the dataset
log close
log using "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-enrol-codebook.txt", text replace
desc
codebook, c
codebook
log close


exit

