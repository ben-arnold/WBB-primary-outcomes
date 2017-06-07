# WBB-primary-outcomes
This is a repository for the WASH Benefits Bangladesh trial ([NCT01590095](https://clinicaltrials.gov/ct2/show/NCT01590095)) primary analysis. The primary outcomes include length for age Z scores and diarrhea. This repository also includes scripts that compile all estimates for the primary outcomes article, including participant enrollment / loss, baseline characteristics, intervention adherence, and select secondary outcomes.

### Associated protocols and datasets

The pre-specified analysis protocol and all data required to run the analyses will be made available in concert with the publication of the article through the Open Science Framework: [https://osf.io/wvyn4](https://osf.io/wvyn4).

The scripts in the `dm` directory process completely raw data and create the processed analysis datasets that are shared publicly through OSF. The completely raw, unprocessed, data are not publicly avaialable at this time, but will be within approximately 2 years time (e.g., by 2019) to allow for the further development of meta-data documentation and an access platform. We will strive to update this page to link to those files when they are available.

### Notes

This was the single effort within the entire WASH Benefits study intention-to-treat analyses that did not use the [`washb`](https://github.com/ben-arnold/washb) package. We developed the package by building off of scripts written for this analysis (mainly in the `basefns` directory).

For all analysis scripts, you will need to change directory statements within them to point them to the files on your local directory. Similar directory statement changes will be needed wherever output files are saved down (e.g., raw estimates, figures). 

### Directory structure

__`adherence`__ : intervention adherence measures at enrollment, year 1, and year 2

__`balance`__ : enrollment characteristics by arm

__`basefns`__ : base functions used in the analysis

__`consort`__ : CONSORT population flow diagram

__`diar`__ : diarrhea analyses

__`dm`__ : data management scripts

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
