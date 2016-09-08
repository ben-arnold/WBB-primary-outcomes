capture log close
set more off
clear all


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/7-bangladesh-hbgdki-data-extract.log", text replace

*--------------------------------------------
* 7-bangladesh-hbgdki-data-extract.do
*
* ben arnold (benarnold@berkeley.edu)
*
* extract data from the WASH Benefits Bangladesh
* trial to share with the HBGDki initiative
*
*--------------------------------------------




*--------------------------------------------
* input files:
*  final/washb-bangladesh-enrol.dta
*  final/washb-bangladesh-anthro.dta
*  final/washb-bangladesh-diar.dta
*  final/washb-bangladesh-tr.dta
*
* output files:
*  HBGDki-extract/washb-bangladesh-hbgdki-enrol.dta / .csv
*  HBGDki-extract/washb-bangladesh-hbgdki-anthro.dta / .csv
*  HBGDki-extract/washb-bangladesh-hbgdki-diar.dta / .csv
*  HBGDki-extract/washb-bangladesh-hbgdki-tr.dta / .csv
*--------------------------------------------

*--------------------------------------------
* loss to follow-up dataset
*--------------------------------------------
use "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-track-compound.dta", clear

keep dataid miss1* miss2*
sort dataid
tempfile track
save `track'


*--------------------------------------------
* enrollment dataset
*--------------------------------------------

use "~/dropbox/WASHB-Bangladesh-Data/1-primary-outcome-datasets/washb-bangladesh-enrol.dta", clear

* drop union names
label drop q4008_lbl

* drop randomiation and treatment assignments
drop block tr

* drop hhid
drop hhid

* merge in the compound tracking information
sort dataid
merge 1:1 dataid using `track'
assert _merge==3
drop _merge

* add some additional notes
notes odch815 : Open defecation among children missing if no children in that age range in household
notes odch38 : Open defecation among children missing if no children in that age range in household
notes odchu3 : Open defecation among children missing if no children in that age range in household
notes asset_clock : Clock measurement missing in 2,859 households

saveold "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-enrol.dta", replace version(12)
outsheet using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-enrol.csv", comma replace
desc
codebook, c

* write a codebook for the dataset
log close
log using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-enrol-codebook.txt", text replace
desc
notes
codebook, c
codebook
log close


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/7-bangladesh-hbgdki-data-extract.log", text append
*--------------------------------------------
* anthropometry dataset
*--------------------------------------------

use "~/dropbox/wbb-primary-analysis/data/final/ben/washb-bangladesh-anthro.dta", clear


* drop randomiation and treatment assignments
drop block tr

* add some additional notes
note: Anthropometry measurements collected at approximately 1 and 2 years after intervention (svy=1|2)
note: Anthropometry Z-scores calculated using the zscore06 Stata macro (laz, waz, whz, bmiz) and the WHO igrowup Stata macro (hcz)
note birthord: Birth order measured at year 2 and is missing for children not present at year 2

saveold "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-anthro.dta", replace version(12)
outsheet using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-anthro.csv", comma replace
desc
codebook, c

* write a codebook for the dataset
log close
log using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-anthro-codebook.txt", text replace
desc
notes
codebook, c
codebook
log close


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/7-bangladesh-hbgdki-data-extract.log", text append
*--------------------------------------------
* diarrhea dataset
*--------------------------------------------

use "~/dropbox/wbb-primary-analysis/data/final/ben/washb-bangladesh-diar.dta", clear

* drop randomiation and treatment assignments
drop block tr



saveold "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-diar.dta", replace version(12)
outsheet using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-diar.csv", comma replace
desc
codebook, c

* write a codebook for the dataset
log close
log using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-diar-codebook.txt", text replace
desc
notes
codebook, c
codebook
log close


log using "~/WASHB-Bangladesh-primary-outcomes/src/dm/7-bangladesh-hbgdki-data-extract.log", text append
*--------------------------------------------
* unblinded treatment assignment dataset
* (password protected - must mount disc img)
*--------------------------------------------
use "/Volumes/0-Treatment-assignments/washb-bangladesh-tr.dta", clear

drop tr
order clusterid block
sort clusterid

saveold "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-tr.dta", replace version(12)
outsheet using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-tr.csv", comma replace
desc
codebook, c

* write a codebook for the dataset
log close
log using "~/dropbox/WBB-primary-analysis/data/HBGDki-extract/washb-bangladesh-hbgdki-tr-codebook.txt", text replace
desc
notes
codebook, c
codebook
log close


exit

