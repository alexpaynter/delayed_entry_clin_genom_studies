#############################################################################################
######### Evaluating Delayed Study Entry & the Effect of Left Truncation Adjustment ######### 
############# in Overall Survival from Diagnosis from a Toy Data Set ########################
#############################################################################################

## Load libraries
library(tidyverse)
library(survival)
library(survminer)
# adding this on for quasi-independence testing:
library(tranSurv)

## Load data
toy_data = read.csv(here::here("toy_data.csv"))

## Variables of interest
## os_status_dx = patient's overall survival status (1 = dead; 0 = censored)
## tt_os_dx_mos = time (months) from date of cancer diagnosis to overall survival, either death or last follow-up
## cpt_status = genomic sequencing indicator (1 = patient underwent genomic sequencing; 0 = patient did not undergo genomic sequencing)
## tt_cpt_report_mos = time (months) from date of cancer diagnosis to date of genomic sequencing

## Calculate median overall survival from diagnosis for total population
survfit(Surv(time = tt_os_dx_mos, 
             event = os_status_dx) ~ 1, 
        data = toy_data)

## Subset data on those who underwent genomic sequencing
toy_data_subset <- toy_data %>% 
  filter(cpt_status == 1)

## Calculate unadjusted median overall survival from diagnosis for subsetted population
survfit(
  Surv(tt_os_dx_mos, 
       os_status_dx) ~ 1, 
  data = toy_data_subset)

## Calculate median overall survival from diagnosis for subsetted population after left truncation adjustment
survfit(
  Surv(time = tt_cpt_report_mos, 
       time2 = tt_os_dx_mos, 
       event = os_status_dx) ~ 1, 
  data = toy_data_subset)





# Quasi-independence testing of time to survival and time to genomic screening.
# Martin and Betensky (2005)
with(toy_data_subset,
     tranSurv::cKendall(trun = tt_cpt_report_mos,
                        obs = tt_os_dx_mos,
                        delta = os_status_dx)
)

# Austin and Betensky (2014) inverse prob weighted version.
with(toy_data_subset,
     tranSurv::cKendall(trun = tt_cpt_report_mos,
                        obs = tt_os_dx_mos,
                        delta = os_status_dx,
                        method = "IPW1")
)

# Austin and Betensky (2014) inverse prob weighted version #2
with(toy_data_subset,
     tranSurv::cKendall(trun = tt_cpt_report_mos,
                        obs = tt_os_dx_mos,
                        delta = os_status_dx,
                        method = "IPW1")
)



