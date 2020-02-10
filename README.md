# WBB-primary-outcomes
This is a repository for the WASH Benefits Bangladesh trial ([NCT01590095](https://clinicaltrials.gov/ct2/show/NCT01590095)) primary analysis. The primary outcomes include length for age Z scores and diarrhea. This repository also includes scripts that compile all estimates for the primary outcomes article, including participant enrollment / loss, baseline characteristics, intervention adherence, and select secondary outcomes.

### Associated protocols and datasets

The pre-specified analysis protocol and all data required to run the analyses will be made available in concert with the publication of the article through the Open Science Framework: [https://osf.io/wvyn4](https://osf.io/wvyn4).

The scripts in the `dm` directory process completely raw data and create the processed analysis datasets that are shared publicly through OSF. The completely raw, unprocessed, data are not publicly avaialable at this time, but will be within approximately 2 years time (e.g., by 2019) to allow for the further development of meta-data documentation and an access platform. We will strive to update this page to link to those files when they are available.

### Notes

This was the single effort within the entire WASH Benefits study intention-to-treat analyses that did not use the [`washb`](https://github.com/ben-arnold/washb) package. We developed the package by building off of scripts written for this analysis (mainly in the `basefns` directory).

Due to R package updates since Jan 2018 used for estimates that rely on pseudo-random number generation, there are slight differences in adjusted estimates and permutation test P-values between published estimates and estimates using the current (February 2020) versions of the coin (version 1.3-1) and tmle (version 1.4.0.1) packages. These differences do not affect inferences.

### Steps to set up directory

1. Clone the WBB-primary-outcomes (https://github.com/ben-arnold/WBB-primary-outcomes) to your local computer, and create an R project within the directory. 

2. Create two folders within the WBB-primary-outcomes directory, one named "data" and another named "results".

3. Download the .csv versions of the public data files from OSF (https://osf.io/pqzj5/), and place them in the "data" folder. The "results" folder will populate with figures, tables, etc. as you run the scripts.

4. Within "src", open the config file "0-config.R". This file contains all relevant packages for the scripts, and is sourced and run at the start of each script. Check that you've installed all the packages listed in the config file.

### Directory structure

__`adherence`__ : intervention adherence measures at enrollment, year 1, and year 2

__`balance`__ : enrollment characteristics by arm

__`basefns`__ : base functions used in the analysis

__`consort`__ : CONSORT population flow diagram

__`diar`__ : diarrhea analyses

__`dm`__ : data management scripts (cannot be run)

__`hcz`__ : head circumference analyses

__`laz`__ : length for age analyses

__`lazminus2`__ : stunting analyses

__`lazminus3`__ : severe stunting analyses

__`mortality`__ : all cause mortality analyses

__`negcontrols`__ : negative control (bruising/abrasion) analyses

__`waz`__ : weight for age analyses

__`wazminus2`__ : underweight analyses

__`whz`__ : weight for height analyses

__`whzminus2`__ : wasting analyses
