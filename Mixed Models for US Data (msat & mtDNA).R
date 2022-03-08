################################################### Script for Mixed Models ########################################################

#analyzes and creates mixed models for mtDNA & msat US datasets 

##########################################################################################################################################

######### Set-up #########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(dplyr)
library(lme4)
library(Hmisc)
library(ggplot2)
library(MuMIn)
library(lmerTest)
library(glmmTMB)
library(tidyr)

#read in data
mtdna_data <- read.csv("new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("new_msat_full_US_data.csv", stringsAsFactors = TRUE) #read in 

#Fixed Variables: Crossspp, repeats, fecundity, body length/maxlength, reproduction mode, fertilization method, primernote
#Random Variables: Species, source/study, marker name
#in brood or similar structure --> internal

################################################### mtDNA data set ################################################### 

#####He##### --> should use separate datasets for Hd & for pi

##Add column for final fertilization
mtdna_data$final_fertilization [mtdna_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization
mtdna_data$final_fertilization [mtdna_data$fertilization =="external"]  <- "external"
mtdna_data$final_fertilization [mtdna_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"

##Create character to factor to integer columns
mtdna_data$fertilizations.or.f <- as.factor(mtdna_data$final_fertilization)
mtdna_data$fertilizations.or.f <- as.numeric(mtdna_data$fertilizations.or.f)

mtdna_data$reproductionmodes.or.f <- as.factor(mtdna_data$reproductionmode)
mtdna_data$reproductionmodes.or.f <- as.numeric(mtdna_data$reproductionmodes.or.f)

##Logtransform what is needed
for (i in 1:nrow(mtdna_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_data$logtransform.fecundity_mean.1 <- log10(mtdna_data$fecundity_mean)
}

for (i in 1:nrow(mtdna_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_data$logtransform.maxlength.1 <- log10(mtdna_data$maxlength)
}

for (i in 1:nrow(mtdna_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_data$logtransform.Pi <- log10(mtdna_data$Pi)
}


##scale bp instead of log-transform
mtdna_data$bp_scale <- scale(as.numeric(mtdna_data$bp))

##Fit model for success/failures##

mtdna_data$success <- round(mtdna_data$He*mtdna_data$n)
mtdna_data$failure <- round((1-mtdna_data$He)*mtdna_data$n)

##Remove NA from variable columns
mtdna_data_He <- mtdna_data %>% drop_na(He, logtransform.maxlength.1, logtransform.fecundity_mean.1,fertilizations.or.f,reproductionmodes.or.f,
                                     bp_scale)

##Create full model for HE
binomial_He_full_model_mtDNA <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                                  fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                  (1|spp) + (1|Source) + (1|MarkerName), na.action = "na.fail", 
                                family=binomial, data = mtdna_data_He,
                                control = glmerControl(optimizer = "bobyqa")) 

mtdna_data_He_dredge <- dredge(binomial_He_full_model_mtDNA) #dredge model
View(mtdna_data_He_dredge) #to get a table that can be copy and pasted to Excel
summary(binomial_He_full_model_mtDNA) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.mtDNAHE <- glmer(formula = cbind(success,failure) ~ fertilizations.or.f + reproductionmodes.or.f +
                                        (1|spp) + (1|Source) + (1|MarkerName), na.action = "na.fail", 
                                      family=binomial, data = mtdna_data_He,
                                      control = glmerControl(optimizer = "bobyqa")) 
mtdna_data_He_dredge_minimal <- dredge(topAIC.mtDNAHE) #dredge model
View(mtdna_data_He_dredge_minimal)
summary(topAIC.mtDNAHE) #get SE, p-value, etc.for top AIC model

#####PI##### --> should use separate datasets for He & for pi

##prep data & remove NA
mtdna_data_nona_fecunditymean <- subset(mtdna_data, mtdna_data$logtransform.fecundity_mean.1 != "NaN")
mtdna_data_nona_fecunditymean_bpscale <- subset(mtdna_data_nona_fecunditymean, mtdna_data_nona_fecunditymean$bp_scale != "NA")
mtdna_pi <- subset(mtdna_data_nona_fecunditymean_bpscale, mtdna_data_nona_fecunditymean_bpscale$logtransform.Pi != "NA") #remove any rows where mtdna pi wasn't calculated

##create full model for pi
Pi_full_model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                        fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                        (1|spp) + (1|Source) + (1|MarkerName), na.action = "na.fail", 
                      data = mtdna_pi, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa")) #have Marial switch to bobyqa to get away from convergence issues AND switch to bp_scale AND make sure REML = FALSE

mtdna_pi_dredge <- dredge(Pi_full_model) #dredge model
View(mtdna_pi_dredge) #to get a table that can be copy and pasted to Excel
summary(Pi_full_model) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.mtDNAPI <- lmer(formula = logtransform.Pi ~ fertilizations.or.f + 
                         (1|spp) + (1|Source) + (1|MarkerName), na.action = "na.fail", 
                       data = mtdna_pi, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa")) #use variables from top AIC model
mtdna_data_Pi_dredge_minimal <- dredge(topAIC.mtDNAPI) #dredge model
View(mtdna_data_Pi_dredge_minimal)
summary(topAIC.mtDNAPI) #get SE, p-value, etc.for top AIC model
summary(mtdna_data_Pi_dredge_minimal)
################################################### msat data set ################################################### 

##Logtransform
for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.fecundity_mean.2 <- log10(msat_data$fecundity_mean)
}

for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.maxlength.2 <- log10(msat_data$maxlength)
}

for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.repeat <- log10(msat_data$Repeat)
}

##Add column for final fertilization

msat_data$final_fertilization [msat_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization
msat_data$final_fertilization [msat_data$fertilization =="external"]  <- "external"
msat_data$final_fertilization [msat_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"

##Create character to factor to integer columns
msat_data$fertilizations.or.f2 <- as.factor(msat_data$final_fertilization)
msat_data$fertilizations.or.f2 <- as.numeric(msat_data$fertilizations.or.f2)

msat_data$reproductionmodes.or.f2 <- as.factor(msat_data$reproductionmode)
msat_data$reproductionmodes.or.f2 <- as.numeric(msat_data$reproductionmodes.or.f2)

##Fit model for success/failures

msat_data$success <- round(msat_data$He*msat_data$n)
msat_data$failure <- round((1-msat_data$He)*msat_data$n)

##Remove NA from variable columns
msat_data <- msat_data %>% drop_na(He, logtransform.maxlength.2, logtransform.fecundity_mean.2,logtransform.repeat,fertilizations.or.f2,reproductionmodes.or.f2)

##Create full model for HE
binomial_He_full_model_msat <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + logtransform.repeat +
                                  fertilizations.or.f2 + reproductionmodes.or.f2 + CrossSpp + PrimerNote +
                                  (1|spp) + (1|Source), na.action = "na.fail", 
                                family=binomial, data = msat_data,
                                control = glmerControl(optimizer = "bobyqa")) #have Marial switch to bobyqa to get away from convergence issues AND switch to bp_scale

msat_dataHe <- dredge(binomial_He_full_model_msat) #dredge model
View(msat_dataHe) #to get a table that can be copy and pasted to Excel
summary(binomial_He_full_model_msat) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.msatHE <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.repeat +
                         fertilizations.or.f2 + CrossSpp +
                         (1|spp) + (1|Source), na.action = "na.fail", 
                       family=binomial, data = msat_data,
                       control = glmerControl(optimizer = "bobyqa"))
summary(topAIC.msatHE) #get SE, p-value, etc.for top AIC model
