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
library(DHARMa)

#read in data
mtdna_data <- read.csv("new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("new_msat_full_US_data.csv", stringsAsFactors = TRUE) #read in
cp_info <- read.csv("spp_combined_info.csv", stringsAsFactors = FALSE)
  mtdna_data <- merge(mtdna_data, cp_info[, c('spp', 'Genus', 'Family')], all.x = TRUE) #add Family & Genus for nested spp RE
  msat_data <- merge(msat_data, cp_info[, c('spp', 'Genus', 'Family')], all.x = TRUE)

#Fixed Variables: Crossspp, repeats, fecundity, body length/maxlength, reproduction mode, fertilization method, primernote
#Random Variables: Species, source/study, marker name
#in brood or similar structure --> internal

################################################### mtDNA data set################################################### 

#Transform He for the gtmmTMB to work
for (i in 1:nrow(mtdna_data)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  mtdna_data$transformed_He[i] <- ((mtdna_data$He[i]*(mtdna_data$n[i] - 1)) + 0.5)/mtdna_data$n[i]
}

#Add column for final fertilization
mtdna_data$final_fertilization [mtdna_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization
mtdna_data$final_fertilization [mtdna_data$fertilization =="external"]  <- "external"
mtdna_data$final_fertilization [mtdna_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"

#Create character to factor to integer columns ########write here what numbers stand for
mtdna_data$fertilizations.or.f <- as.factor(mtdna_data$final_fertilization) #external = 1, internal (oviduct) = 2
mtdna_data$fertilizations.or.f <- as.numeric(mtdna_data$fertilizations.or.f)

mtdna_data$reproductionmodes.or.f <- as.factor(mtdna_data$reproductionmode) #dioecism = 1, protogyny = 2
mtdna_data$reproductionmodes.or.f <- as.numeric(mtdna_data$reproductionmodes.or.f)

#Logtransform what is needed
for (i in 1:nrow(mtdna_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_data$logtransform.fecundity_mean.1 <- log10(mtdna_data$fecundity_mean)
}

mtdna_data <- subset(mtdna_data, mtdna_data$logtransform.fecundity_mean.1 != "NaN")

for (i in 1:nrow(mtdna_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_data$logtransform.maxlength.1 <- log10(mtdna_data$maxlength)
}

#scale bp instead of log-transform
mtdna_data$bp_scale <- scale(as.numeric(mtdna_data$bp))

for (i in 1:nrow(mtdna_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_data$logtransform.Pi <- log10(mtdna_data$Pi)
}

#######################################################################################################

## model for mtdna Hd ## --> should use separate datasets for Hd & for pi
##Fit model for success/failures##

#prep data
mtdna_data_nona_fecunditymean <- subset(mtdna_data, mtdna_data$logtransform.fecundity_mean.1 != "NaN")
mtdna_Hd <- subset(mtdna_data, mtdna_data_nona_fecunditymean$He != "NA") #remove any rows where mtdna Hd wasn't calculated


mtdna_Hd$success <- round(mtdna_Hd$He*mtdna_Hd$n)
mtdna_Hd$failure <- round((1-mtdna_Hd$He)*mtdna_Hd$n)

#null (full) model for Hd --> no site for mtdna bc essentially observation-level RE
summary(mtdna_data)

binomial_Hd_full_model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                 fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                 (1|spp) + (1|Source), na.action = "na.fail", 
               family=binomial, data = mtdna_Hd,
               control = glmerControl(optimizer = "bobyqa")) #have Marial switch to bobyqa to get away from convergence issues AND switch to bp_scale

#dharma binomial
binomial_Hd_full_sim <- simulateResiduals(fittedModel = binomial_Hd_full_model, n = 1000, plot = F)
plotQQunif(binomial_Hd_full_sim)
plotResiduals(binomial_Hd_full_model)
testDispersion(binomial_Hd_full_sim)

mtdna_Hd_nona_fecunditymean <- subset(mtdna_Hd, mtdna_Hd$logtransform.fecundity_mean.1 != "NaN")
plotResiduals(binomial_Hd_full_sim, mtdna_Hd_nona_fecunditymean$logtransform.maxlength.1)
plotResiduals(binomial_Hd_full_sim, mtdna_Hd_nona_fecunditymean$logtransform.fecundity_mean.1)
plotResiduals(binomial_Hd_full_sim, mtdna_Hd_nona_fecunditymean$fertilizations.or.f)
plotResiduals(binomial_Hd_full_sim, mtdna_Hd_nona_fecunditymean$reproductionmodes.or.f)
plotResiduals(binomial_Hd_full_sim, mtdna_Hd_nona_fecunditymean$bp_scale)

beta_Hd_full_model <- glmmTMB(transformed_He ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                                fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                (1|spp) + (1|Source) + (1|MarkerName), 
                              family = beta_family, data = mtdna_Hd)

#dharma beta
beta_Hd_full_sim <- simulateResiduals(fittedModel = beta_Hd_full_model, n = 1000, plot = F)
plotQQunif(beta_Hd_full_sim)
plotResiduals(beta_Hd_full_model)
testDispersion(beta_Hd_full_sim)

plotResiduals(beta_Hd_full_sim, mtdna_Hd_nona_fecunditymean$logtransform.maxlength.1)
plotResiduals(beta_Hd_full_sim, mtdna_Hd_nona_fecunditymean$logtransform.fecundity_mean.1)
plotResiduals(beta_Hd_full_sim, mtdna_Hd_nona_fecunditymean$fertilizations.or.f)
plotResiduals(beta_Hd_full_sim, mtdna_Hd_nona_fecunditymean$reproductionmodes.or.f)
plotResiduals(beta_Hd_full_sim, mtdna_Hd_nona_fecunditymean$bp_scale)

#dredge models
binomial_mtdna_Hd_dredge <- dredge(binomial_Hd_full_model)
beta_mtdna_Hd_dredge <- dredge(beta_Hd_full_model)

#didn't converge --> converged outside of dredge (one of top models)
beta_Hd_check_model <- glmmTMB(transformed_He ~ 
                                fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                (1|spp) + (1|Source) + (1|MarkerName), 
                              family = beta_family, data = mtdna_Hd)

#check fits of top models
beta_Hd_top_model <- glmmTMB(transformed_He ~ 
                                 fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                 (1|spp) + (1|Source) + (1|MarkerName), 
                               family = beta_family, data = mtdna_Hd)

beta_Hd_top_sim <- simulateResiduals(fittedModel = beta_Hd_top_model, n = 1000, plot = F)
plotQQunif(beta_Hd_top_sim)
plotResiduals(beta_Hd_top_model)
testDispersion(beta_Hd_top_sim)

#check fits of top models
binomial_Hd_top_model <- glmer(formula = cbind(success,failure) ~  
                                            fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                            (1|spp) + (1|Source), 
                                          family=binomial, data = mtdna_Hd,  na.action = 'na.fail', 
                                          control = glmerControl(optimizer = "bobyqa"))

binomial_Hd_top_sim <- simulateResiduals(fittedModel = binomial_Hd_top_model, n = 1000, plot = F)
plotQQunif(binomial_Hd_top_sim)
plotResiduals(binomial_Hd_top_sim)
testDispersion(binomial_Hd_top_sim)

#########################################################################################################################

## model for mtdna pi ## --> should use separate datasets for Hd & for pi

#prep data
mtdna_data_nona_fecunditymean <- subset(mtdna_data, mtdna_data$logtransform.fecundity_mean.1 != "NaN")
mtdna_data_nona_fecunditymean_bpscale <- subset(mtdna_data_nona_fecunditymean, mtdna_data_nona_fecunditymean$bp_scale != "NA")
mtdna_pi <- subset(mtdna_data_nona_fecunditymean_bpscale, mtdna_data_nona_fecunditymean_bpscale$logtransform.Pi != "NA") #remove any rows where mtdna pi wasn't calculated
mtdna_pi

#null (full) model for pi
Pi_full_model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                                  fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                  (1|spp) + (1|Source), na.action = "na.fail", 
                                data = mtdna_pi, REML = FALSE,
                                control = lmerControl(optimizer = "bobyqa")) #have Marial switch to bobyqa to get away from convergence issues AND switch to bp_scale AND make sure REML = FALSE

#pull p-values
coefs <- data.frame(coef(summary(Pi_full_model)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#dharma linear pi
Pi_full_sim <- simulateResiduals(fittedModel = Pi_full_model, n = 1000, plot = F)
plotQQunif(Pi_full_sim)
plotResiduals(Pi_full_model)
testDispersion(Pi_full_sim)

plotResiduals(Pi_full_sim, mtdna_pi$logtransform.maxlength.1)
plotResiduals(Pi_full_sim, mtdna_pi$logtransform.fecundity_mean.1)
plotResiduals(Pi_full_sim, mtdna_pi$fertilizations.or.f)
plotResiduals(Pi_full_sim, mtdna_pi$reproductionmodes.or.f)
plotResiduals(Pi_full_sim, mtdna_pi$bp_scale)

#dredge models
mtdna_Pi_dredge <- dredge(Pi_full_model)

#check fits of top models
Pi_top_model <- lmer(formula = logtransform.Pi ~  
                                 fertilizations.or.f + logtransform.maxlength.1 + 
                                 (1|spp) + (1|Source), 
                               data = mtdna_pi,  REML = FALSE, na.action = 'na.fail', 
                               control = lmerControl(optimizer = "bobyqa"))

Pi_top_sim <- simulateResiduals(fittedModel = Pi_top_model, n = 1000, plot = F)
plotQQunif(Pi_top_sim)
plotResiduals(Pi_top_sim)
testDispersion(Pi_top_sim)

################################################### msat data set ################################################### 

#Transform He for the gtmmTMB to work
for (i in 1:nrow(msat_data)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  msat_data$transformed_He[i] <- ((msat_data$He[i]*(msat_data$n[i] - 1)) + 0.5)/msat_data$n[i]
}

#Logtransform
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

#Add column for final fertilization

msat_data$final_fertilization [msat_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization
msat_data$final_fertilization [msat_data$fertilization =="external"]  <- "external"
msat_data$final_fertilization [msat_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"

#Create character to factor to integer columns
msat_data$fertilizations.or.f2 <- as.factor(msat_data$final_fertilization)
msat_data$fertilizations.or.f2 <- as.numeric(msat_data$fertilizations.or.f2)

msat_data$reproductionmodes.or.f2 <- as.factor(msat_data$reproductionmode)
msat_data$reproductionmodes.or.f2 <- as.numeric(msat_data$reproductionmodes.or.f2)

#######################################################################################################

## model for msat He## 

#prep data
msat_data <- msat_data[!is.na(msat_data$Repeat), ]
msat_data <- msat_data[!is.na(msat_data$fertilizations.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$reproductionmodes.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$logtransform.maxlength.2), ]
msat_data <- msat_data[!is.na(msat_data$logtransform.fecundity_mean.2), ]
msat_data <- msat_data[!is.na(msat_data$PrimerNote), ]
msat_data <- msat_data[!is.na(msat_data$CrossSpp), ]
msat_data$ID <- c(1:2163)


msat_data$success <- round(msat_data$He*msat_data$n)
msat_data$failure <- round((1-msat_data$He)*msat_data$n)

#null (full) model for He
binomial_He_full_model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + 
                                  fertilizations.or.f2 + reproductionmodes.or.f2 + CrossSpp +
                                  (1|spp) + (1|Source) + (1|ID), na.action = "na.fail", 
                                family=binomial, data = msat_data,
                                control = glmerControl(optimizer = "bobyqa")) #have Marial switch to bobyqa to get away from convergence issues

#dharma binomial
binomial_He_full_sim <- simulateResiduals(fittedModel = binomial_He_full_model, n = 1000, plot = F)
plotQQunif(binomial_He_full_sim)
plotResiduals(binomial_He_full_model)
testDispersion(binomial_He_full_sim)

plotResiduals(binomial_He_full_sim, msat_data$logtransform.maxlength.2)
plotResiduals(binomial_He_full_sim, msat_data$logtransform.fecundity_mean.2)
plotResiduals(binomial_He_full_sim, msat_data$fertilizations.or.f2)
plotResiduals(binomial_He_full_sim, msat_data$reproductionmodes.or.f2)
plotResiduals(binomial_He_full_sim, msat_data$Repeat)
plotResiduals(binomial_He_full_sim, msat_data$PrimerNote)
plotResiduals(binomial_He_full_sim, msat_data$CrossSpp)

beta_He_full_model <- glmmTMB(transformed_He ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + 
                                fertilizations.or.f2 + reproductionmodes.or.f2 + CrossSpp + 
                                (1|spp) + (1|Source), 
                              family = beta_family, data = msat_data)

#dharma beta
beta_He_full_sim <- simulateResiduals(fittedModel = beta_He_full_model, n = 1000, plot = F)
plotQQunif(beta_He_full_sim)
plotResiduals(beta_He_full_model)
testDispersion(beta_He_full_sim)

plotResiduals(beta_He_full_sim, msat_data$logtransform.maxlength.2)
plotResiduals(beta_He_full_sim, msat_data$logtransform.fecundity_mean.2)
plotResiduals(beta_He_full_sim, msat_data$fertilizations.or.f2)
plotResiduals(beta_He_full_sim, msat_data$reproductionmodes.or.f2)
plotResiduals(beta_He_full_sim, msat_data$Repeat)
plotResiduals(beta_He_full_sim, msat_data$CrossSpp)
plotResiduals(beta_He_full_sim, msat_data$PrimerNote)

#dredge models
binomial_msat_He_dredge <- dredge(binomial_He_full_model)
beta_msat_He_dredge <- dredge(beta_He_full_model)

#check fits of top models
beta_He_top_model <- glmmTMB(transformed_He ~ 
                               fertilizations.or.f2 + logtransform.maxlength.2 + CrossSpp + 
                               (1|spp) + (1|Source), 
                             family = beta_family, data = msat_data)

beta_He_top_sim <- simulateResiduals(fittedModel = beta_He_top_model, n = 1000, plot = F)
plotQQunif(beta_He_top_sim)
plotResiduals(beta_He_top_model)
testDispersion(beta_He_top_sim)

#check fits of top models
binomial_He_top_model <- glmer(formula = cbind(success,failure) ~  
                                 fertilizations.or.f2 + logtransform.maxlength.2 + CrossSpp + 
                                 (1|spp) + (1|Source) + (1|ID), 
                               family=binomial, data = msat_data,  na.action = 'na.fail', 
                               control = glmerControl(optimizer = "bobyqa"))

binomial_He_top_sim <- simulateResiduals(fittedModel = binomial_He_top_model, n = 1000, plot = F)
plotQQunif(binomial_He_top_sim)
plotResiduals(binomial_He_top_sim)
testDispersion(binomial_He_top_sim)



