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

################################################### mtDNA data set################################################### 

#####He##### --> should use separate datasets for Hd & for pi

#Transform He for the gtmmTMB to work
for (i in 1:nrow(mtdna_data)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  mtdna_data$transformed_He[i] <- ((mtdna_data$He[i]*(mtdna_data$n[i] - 1)) + 0.5)/mtdna_data$n[i]
}

#Add column for final fertilization
mtdna_data$final_fertilization [mtdna_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization
mtdna_data$final_fertilization [mtdna_data$fertilization =="external"]  <- "external"
mtdna_data$final_fertilization [mtdna_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"

#Create character to factor to integer columns
mtdna_data$fertilizations.or.f <- as.factor(mtdna_data$final_fertilization)
mtdna_data$fertilizations.or.f <- as.numeric(mtdna_data$fertilizations.or.f)

mtdna_data$reproductionmodes.or.f <- as.factor(mtdna_data$reproductionmode)
mtdna_data$reproductionmodes.or.f <- as.numeric(mtdna_data$reproductionmodes.or.f)

#Logtransform what is needed
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


#scale bp instead of log-transform
mtdna_data$bp_scale <- scale(as.numeric(mtdna_data$bp))

##Fit model for success/failures##

mtdna_data$success <- round(mtdna_data$He*mtdna_data$n)
mtdna_data$failure <- round((1-mtdna_data$He)*mtdna_data$n)

mtdna_data <- mtdna_data %>% drop_na(He, logtransform.maxlength.1, logtransform.fecundity_mean.1,fertilizations.or.f,reproductionmodes.or.f,
                                     bp_scale)

binomial_He_full_model_mtDNA <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                                  fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                  (1|spp) + (1|Source), na.action = "na.fail", 
                                family=binomial, data = mtdna_data,
                                control = glmerControl(optimizer = "bobyqa")) 

mtdna_data <- dredge(binomial_He_full_model_mtDNA)
View(mtdna_data) #to get a table that can be copy and pasted to Excel
summary(binomial_He_full_model_mtDNA) #get SE, p-value, etc.

#####PI##### --> should use separate datasets for He & for pi

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

mtdna_pi <- dredge(Pi_full_model)
View(mtdna_pi) #to get a table that can be copy and pasted to Excel
summary(Pi_full_model)

#dredge for pi
model <- glm(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1, data = mtdna_data, na.action = 'na.fail')
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp) + (1|Source)+ (1|MarkerName), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
# REMOVE: model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Source), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName), data = mtdna_data, na.action = 'na.fail',REML=FALSE)

#combo of random variable for pi 
# REMOVE: model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site) + (1|Source), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp) + (1|Source) + (1|MarkerName), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
# REMOVE: model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site) +(1|Source) + (1|MarkerName), data = mtdna_data, na.action = 'na.fail', REML=FALSE)


mtdna.spp.Pi <- dredge(model)
View(mtdna.spp.Pi) #to get a table that can be copy and pasted to Excel

#Check top 4 models individually to see if there were convergence problems for combo of all random variables
fit.bin <- lmer(logtransform.Pi ~ fertilizations.or.f + (1|spp) + (1|Source)+ (1|MarkerName), data=mtdna_data, REML=FALSE)
fit.bin <- lmer(logtransform.Pi ~ logtransform.maxlength.1 + fertilizations.or.f + (1|spp) + (1|Source)+ (1|MarkerName), data=mtdna_data, REML=FALSE)
fit.bin <- lmer(logtransform.Pi ~ fertilizations.or.f + logtransform.bp.1 + (1|spp) + (1|Source)+ (1|MarkerName), data=mtdna_data, REML=FALSE)
fit.bin <- lmer(logtransform.Pi ~ fertilizations.or.f + reproductionmodes.or.f + (1|spp) + (1|Source)+ (1|MarkerName), data=mtdna_data, REML=FALSE)


#Find SE & P-Value for Fixed, SD for Random
sum <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site) + (1|Source), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
sum <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Source) + (1|MarkerName), data = mtdna_data, na.action = 'na.fail', REML=FALSE)

summary(sum)


################################################### msat data set ################################################### 

#Transform He for the gtmmTMB to work
for (i in 1:nrow(msat_data)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  msat_data$transformed_He2[i] <- ((msat_data$He[i]*(msat_data$n[i] - 1)) + 0.5)/msat_data$n[i]
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

##Fit model for success/failures

msat_data$success <- round(msat_data$He*msat_data$n)
msat_data$failure <- round((1-msat_data$He)*msat_data$n)

msat_data <- msat_data %>% drop_na(He, logtransform.maxlength.2, logtransform.fecundity_mean.2,logtransform.repeat,fertilizations.or.f2,reproductionmodes.or.f2)

binomial_He_full_model_msat <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + logtransform.repeat +
                                  fertilizations.or.f2 + reproductionmodes.or.f2 + CrossSpp + PrimerNote +
                                  (1|spp) + (1|Source), na.action = "na.fail", 
                                family=binomial, data = msat_data,
                                control = glmerControl(optimizer = "bobyqa")) #have Marial switch to bobyqa to get away from convergence issues AND switch to bp_scale


msat_dataHe <- dredge(binomial_He_full_model_msat)
View(msat_dataHe) #to get a table that can be copy and pasted to Excel
summary(binomial_He_full_model_msat)

#####################################################################
#dredge? 
msat_data <- msat_data[!is.na(msat_data$Repeat), ]
msat_data <- msat_data[!is.na(msat_data$fertilizations.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$reproductionmodes.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$logtransform.maxlength.2), ]
msat_data <- msat_data[!is.na(msat_data$logtransform.fecundity_mean.2), ]
msat_data <- msat_data[!is.na(msat_data$PrimerNote), ]
msat_data <- msat_data[!is.na(msat_data$CrossSpp), ]

#summary
model <- glm(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp, data = msat_data, family=binomial, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + PrimerNote + (1|spp) + (1|Source)+ (1|PrimerNote), family=binomial, data = msat_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + PrimerNote + (1|spp), family=binomial, data = msat_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + PrimerNote + (1|Source), family=binomial, data = msat_data, na.action = 'na.fail')

#combo of random variable for He
model <- glmmTMB(transformed_He2 ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + PrimerNote + CrossSpp + Repeat + (1|spp) + (1|Source), family = beta_family, data = msat_data)

msat.spp <- dredge(model)
summary(model)

View(msat.spp)
summary_mtdna.spp <-summary(msat.spp)

#Check top 4 models individually to see if there were convergence problems for combo of all random variables
fit.bin <- glmer(cbind(success,failure) ~ logtransform.maxlength.2 + fertilizations.or.f2 + Repeat + CrossSpp + (1|spp) + (1|Source)+ (1|PrimerNote), family=binomial, data=msat_data)
fit.bin <- glmer(cbind(success,failure) ~ fertilizations.or.f2 + Repeat + CrossSpp + (1|spp) + (1|Source)+ (1|PrimerNote), family=binomial, data=msat_data)
fit.bin <- glmer(cbind(success,failure) ~ logtransform.maxlength.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|spp) + (1|Source)+ (1|PrimerNote), family=binomial, data=msat_data)
fit.bin <- glmer(cbind(success,failure) ~ fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|spp) + (1|Source)+ (1|PrimerNote), family=binomial, data=msat_data)
fit.bin <- glmer(cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + Repeat + CrossSpp +  (1|spp) + (1|Source)+ (1|PrimerNote), family=binomial, data=msat_data)

#Find SE & P-Value for Fixed, SD for Random
sum <- glm(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp, data = msat_data, family=binomial, na.action = 'na.fail')
summary(sum)
