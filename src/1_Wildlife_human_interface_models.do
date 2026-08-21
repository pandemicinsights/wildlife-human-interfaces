/******************************************************************
This file contain STATA code for replication of analyses from the publication [add details here]. 

Multivariable statistical analyses to evaluate infection prevalence in animals. Analyses were conducted in STATA (version 16.1, StataCorp, College Station, Texas, USA) using data provided in Source Data 1 “data1_viruses_bats_rodents.csv”. See Methods Section, subheading Analytical Procedures, for more information on factor inclusion criteria and model selection.

This code is specific to: 
    Figures 1b and 1c
    Supplementary Tables 2.2-2.6

See the README file for how to download source data: 
https://github.com/pandemicinsights/wildlife-human-interfaces
******************************************************************/


/******************************************************************
a) Multivariable model investigating factors related to detection of coronaviruses paramyxoviruses, flaviviruses, and influenza viruses among bats, rodents, and shrews combined as presented in Figures 1b and 1c and Supplementary Table 2.2.
******************************************************************/

melogit cr_rna i_bats i_wet i_feces i_guano i_rectal_anal_swab i_urine dwellings raiding_crops raiding_markets ecotourism wildlife_management guano_farm hunted consumption private_sale transit_along_valuechain for_sale_in_small_market medium_market large_markets zoo_sanctuary|| _all: R.labname, or


/******************************************************************
b) Multivariable model predicting factors related to detection of coronaviruses among bats only as presented in Figures 1b and 1c and Supplementary Table 2.3.
******************************************************************/

melogit cr_coronaviridae i_wet i_feces i_guano i_rectal_anal_swab i_urine dwellings raiding_crops guano_farm private_sale for_sale_in_small_market large_markets transit_along_valuechain zoo_sanctuary if i_rodents_shrews ==0 & familytest_coronaviridae==1 || _all: R.labname, or


/******************************************************************
c) Multivariable mixed effects logistic regression models predicting related to the detection of paramyxoviruses among bats as presented in Figures 1b and 1c and Supplementary Table 2.4.
******************************************************************/

melogit cr_paramyxoviridae i_wet i_guano i_rectal_anal_swab i_urine i_urogenital_swab dwellings  large_markets zoo_sanctuary if i_rodents_shrews ==0 & familytest_paramyxoviridae==1 || _all: R.labname, or


/******************************************************************
d) Multivariable mixed effects logistic regression models predicting related to the detection of coronaviruses among rodents and shrews as presented in Figures 1b and 1c and Supplementary Table 2.5.
******************************************************************/

melogit cr_coronaviridae i_wet i_feces i_rectal_anal_swab dwellings transit_along_valuechain medium_market large_markets if i_bats==0 & familytest_coronaviridae==1 || _all: R.labname, or


/******************************************************************
e) Multivariable mixed effects logistic regression models predicting related to the detection of paramyxoviruses among rodents and shrews as presented in Supplementary Table 2.6.
******************************************************************/

melogit cr_paramyxoviridae i_feces i_rectal_anal_swab if i_bats==0 & familytest_paramyxoviridae || _all: R.labname, or



