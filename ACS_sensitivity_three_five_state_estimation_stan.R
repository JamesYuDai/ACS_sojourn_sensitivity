###
### This is final code for the ACS manuscript to implement rstan for acs analysis
### Sep 2025
###

rm(list=ls())
library("rstan")
library("tidyverse")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
setwd("/Users/jdai/Documents/Sensitivity_estimation/")

cps_dat <- read.delim("/Users/jdai/Documents/CPS-3/datasets/acs_adsl_primary.tsv",header=T)

summary(cps_dat |> filter(cancer_status == "Cancer") |> select(weight))
table(cps_dat |> select(dx_year))

###
###
### all cancers 
### not limit to age >50
### Fit the three state model
###
###

cdat <- cps_dat |> filter(cancer_status=="Cancer", cstage !="In Situ" ) |> select("ctype","days_to_dx","detection","weight", "cstage")
nrow(cdat)
cdat <- cdat |> mutate (days_to_dx=days_to_dx/365)
summary(cdat|> filter(detection =="Yes") |> select (days_to_dx))

cdat <- cdat |> mutate(detection = 1*(detection=="Yes"))
names(cdat) <- c("ctype","dtime","tpos","weight","cstage")


standata <- list(sens_prior=c(0.4,0.6),  # or use c(0.1,0.9) for not informative prior
                 lamb_prior=c(0.25,0.25),
                 dtime=drop(cdat$dtime),
                 tpos=drop(cdat$tpos),
                 wgt=drop(cdat$weight),
                 J=nrow(cdat))

fit_stan <- stan("ACS_sensitivity.stan",data=standata,chains = 4, iter = 8000)
print(fit_stan)
print(fit_stan, pars = "sens")
print(fit_stan, pars = "lamb")
pairs(fit_stan,pars=c("sens","mst"))
traceplot(fit_stan, pars = c("sens", "lamb"), inc_warmup = TRUE, nrow = 2)









###
###  top 12 lethal cancers 
###  not limited by age >50
###  Fit the three-state model 
###

fatal_cancers <- c("Anus","Bladder","Colon/Rectum","Esophagus","Head and Neck","Liver/Bile-duct","Lung","Lymphoma","Ovary","Pancreas","Plasma Cell Neoplasm", "Stomach")

cdat <- cps_dat |> filter(cancer_status=="Cancer", cstage !="In Situ", ctype %in% fatal_cancers) |> select("ctype","days_to_dx","detection","weight", "cstage")
nrow(cdat)

nrow(cdat)
cdat <- cdat |> mutate (days_to_dx=days_to_dx/365)
summary(cdat|> filter(detection =="Yes") |> select (days_to_dx))

cdat <- cdat |> mutate(detection = 1*(detection=="Yes"))
names(cdat) <- c("ctype","dtime","tpos","weight","cstage")







standata <- list(sens_prior=c(0.4,0.6),   # or use c(0.1,0.9) for not informative prior
                 lamb_prior=c(0.25,0.25),
                 dtime=drop(cdat$dtime),
                 tpos=drop(cdat$tpos),
                 wgt=drop(cdat$weight),
                 J=nrow(cdat))

fit_stan <- stan("ACS_sensitivity.stan",data=standata,chains = 4, iter = 8000)
print(fit_stan)
print(fit_stan, pars = "sens")
print(fit_stan, pars = "lamb")
pairs(fit_stan,pars=c("sens","mst"))
traceplot(fit_stan, pars = c("sens", "lamb"), inc_warmup = TRUE, nrow = 2)








###
### prepare the data for stan
### all cancers; not anchoring age >50
### but still using weight
###


cdat <- cps_dat |> filter(cancer_status=="Cancer") |>
  select("ctype","days_to_dx","detection","weight","cstage")
cdat <- cdat |> filter(cdat$cstage != "In Situ" & cdat$cstage != "Non-informative") |>
  mutate(days_to_dx=days_to_dx/365) |>
  mutate(detection=1*(detection=="Yes"))


###
### use CCGA3 informed prior 
###

names(cdat) <- c("ctype","dtime","tpos","weight","cstage")
standata <- list(sens_early_prior=c(0.3,0.7),
                 sens_late_prior=c(7,3),
                 lamb_sum_prior=c(0.25,0.25),
                 lamb_end_prior=c(0.25,0.25),
                 dtime_early=drop(cdat$dtime[cdat$cstage!="Distant"]),
                 tpos_early=drop(cdat$tpos[cdat$cstage!="Distant"]),
                 wgt_early=drop(cdat$weight[cdat$cstage!="Distant"]),
                 dtime_late=drop(cdat$dtime[cdat$cstage=="Distant"]),
                 tpos_late=drop(cdat$tpos[cdat$cstage=="Distant"]),
                 wgt_late=drop(cdat$weight[cdat$cstage=="Distant"]),
                 J1=nrow(cdat[cdat$cstage!="Distant",]),
                 J2=nrow(cdat[cdat$cstage=="Distant",]))

sum(cdat$cstage!="Distant")
sum(cdat$cstage=="Distant")

fit_stan <- stan("ACS_sensitivity_early_late.stan",data=standata,chains = 4, iter = 8000)
print(fit_stan)
pairs(fit_stan,pars=c("sens_early","sens_late", "mst1", "mst2"))





###
### use uninformative prior 
###

names(cdat) <- c("ctype","dtime","tpos","weight","cstage")
standata <- list(sens_early_prior=c(0.1,0.9),
                 sens_late_prior=c(0.5,0.5),
                 lamb_sum_prior=c(0.25,0.25),
                 lamb_end_prior=c(0.25,0.25),
                 dtime_early=drop(cdat$dtime[cdat$cstage!="Distant"]),
                 tpos_early=drop(cdat$tpos[cdat$cstage!="Distant"]),
                 wgt_early=drop(cdat$weight[cdat$cstage!="Distant"]),
                 dtime_late=drop(cdat$dtime[cdat$cstage=="Distant"]),
                 tpos_late=drop(cdat$tpos[cdat$cstage=="Distant"]),
                 wgt_late=drop(cdat$weight[cdat$cstage=="Distant"]),
                 J1=nrow(cdat[cdat$cstage!="Distant",]),
                 J2=nrow(cdat[cdat$cstage=="Distant",]))

sum(cdat$cstage!="Distant")
sum(cdat$cstage=="Distant")

fit_stan <- stan("ACS_sensitivity_early_late.stan",data=standata,chains = 4, iter = 4000)
print(fit_stan)







###
###
### top 12 lethal cancers 
###
###


fatal_cancers <- c("Anus","Bladder","Colon/Rectum","Esophagus","Head and Neck","Liver/Bile-duct","Lung","Lymphoma","Ovary","Pancreas","Plasma Cell Neoplasm", "Stomach")

cdat <- cps_dat |> filter(cancer_status=="Cancer",ctype %in% fatal_cancers) |>
  select("ctype","days_to_dx","detection","weight","cstage")
cdat <- cdat |> filter(cdat$cstage != "In Situ" & cdat$cstage != "Non-informative") |>
  mutate(days_to_dx=days_to_dx/365) |>
  mutate(detection=1*(detection=="Yes"))



names(cdat) <- c("ctype", "dtime","tpos","weight","cstage")
standata <- list(sens_early_prior=c(0.3,0.7),
                 sens_late_prior=c(7,3),
                 lamb_sum_prior=c(0.25,0.25),
                 lamb_end_prior=c(0.25,0.25),
                 dtime_early=drop(cdat$dtime[cdat$cstage!="Distant"]),
                 tpos_early=drop(cdat$tpos[cdat$cstage!="Distant"]),
                 wgt_early=drop(cdat$weight[cdat$cstage!="Distant"]),
                 dtime_late=drop(cdat$dtime[cdat$cstage=="Distant"]),
                 tpos_late=drop(cdat$tpos[cdat$cstage=="Distant"]),
                 wgt_late=drop(cdat$weight[cdat$cstage=="Distant"]),
                 J1=nrow(cdat[cdat$cstage!="Distant",]),
                 J2=nrow(cdat[cdat$cstage=="Distant",]))

sum(cdat$cstage!="Distant")
sum(cdat$cstage=="Distant")

fit_stan <- stan("ACS_sensitivity_early_late.stan",data=standata,chains = 4, iter = 8000)
print(fit_stan)




###
### Use uninformative prior
### 

standata <- list(sens_early_prior=c(0.1,0.9),
                 sens_late_prior=c(0.5,0.5),
                 lamb_sum_prior=c(0.25,0.25),
                 lamb_end_prior=c(0.25,0.25),
                 dtime_early=drop(cdat$dtime[cdat$cstage!="Distant"]),
                 tpos_early=drop(cdat$tpos[cdat$cstage!="Distant"]),
                 wgt_early=drop(cdat$weight[cdat$cstage!="Distant"]),
                 dtime_late=drop(cdat$dtime[cdat$cstage=="Distant"]),
                 tpos_late=drop(cdat$tpos[cdat$cstage=="Distant"]),
                 wgt_late=drop(cdat$weight[cdat$cstage=="Distant"]),
                 J1=nrow(cdat[cdat$cstage!="Distant",]),
                 J2=nrow(cdat[cdat$cstage=="Distant",]))

sum(cdat$cstage!="Distant")
sum(cdat$cstage=="Distant")

fit_stan <- stan("ACS_sensitivity_early_late.stan",data=standata,chains = 4, iter = 4000)
print(fit_stan)




