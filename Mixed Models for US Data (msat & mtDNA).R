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

#read in data
mtdna_data <- read.csv("new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("msat_full_US_data.csv", stringsAsFactors = FALSE) #read in 

#Fixed Variables: Crossspp, repeats, fecundity, body length/maxlength, reproduction mode, fertilization method
#Random Variables: Species, site, source/study, marker name
#in brood or similar structure --> internal

############### mtDNA data set ############### 

##Create variable to model success/failures
mtdna_data$mtdnas.or.f <- round(mtdna_data$He*mtdna_data$n) & round((1-mtdna_data$He)*mtdna_data$n)  #create column of successes or failures

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
  mtdna_data$logtransform.bp.1 <- log10(mtdna_data$bp)
}

##Fit model for success/failures##

#dredge for He 
mtdna_data <- mtdna_data[!is.na(mtdna_data$bp), ]
mtdna_data <- mtdna_data[!is.na(mtdna_data$mtdnas.or.f), ]
mtdna_data <- mtdna_data[!grepl('NaN',mtdna_data$logtransform.fecundity_mean.1),]
mtdna_data$fecundity_mean <- NULL
mtdna_data$bp <- NULL
summary(mtdna_data)
model <- glm(formula = mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1, family=binomial, data = mtdna_data,na.action = 'na.fail')
model <- glmer(formula = mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp) + (1|Site) + (1|Source)+ (1|MarkerName), family=binomial, data = mtdna_data)
model <- glmer(formula = mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp), family=binomial, data = mtdna_data)
model <- glmer(formula = mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site), family=binomial, data = mtdna_data, na.action = 'na.fail')
model <- glmer(formula = mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Source), family=binomial, data = mtdna_data, na.action = 'na.fail')
model <- glmer(formula = mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName), family=binomial, data = mtdna_data)

mtdna.spp <- dredge(model)
summary_mtdna.spp <-summary(mtdna.spp)

#dredge for pi
model <- glm(formula = Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1, data = mtdna_data,na.action = 'na.fail')
model <- lmer(formula = Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp) + (1|Site) + (1|Source)+ (1|MarkerName), data = mtdna_data, na.action = 'na.fail')
model <- lmer(formula = Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp), data = mtdna_data, na.action = 'na.fail')
model <- lmer(formula = Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site), data = mtdna_data, na.action = 'na.fail')
model <- lmer(formula = Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Source), data = mtdna_data, na.action = 'na.fail')
model <- lmer(formula = Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName), data = mtdna_data, na.action = 'na.fail')

mtdna.spp.He <- dredge(model)

#spp as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + (1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + (1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ fertilizations.or.f + (1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ reproductionmodes.or.f + (1|spp), family=binomial, data=mtdna_data)

#site as Randoom variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ fertilizations.or.f + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ reproductionmodes.or.f + (1|Site), family=binomial, data=mtdna_data)

#source/study as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ fertilizations.or.f + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ reproductionmodes.or.f + (1|Source), family=binomial, data=mtdna_data)

#MarkerName as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + (1|MarkerName), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + (1|MarkerName), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ fertilizations.or.f + (1|MarkerName), family=binomial, data=mtdna_data) ####CONVERGENCE PROBLEM

fit.bin <- glmer(mtdnas.or.f ~ reproductionmodes.or.f + (1|MarkerName), family=binomial, data=mtdna_data)

##More Mixed Models with multiple fixed variables (max length constant for the following code)##

#spp as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + (1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + fertilizations.or.f + (1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + reproductionmodes.or.f + (1|spp), family=binomial, data=mtdna_data)

#site as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + fertilizations.or.f + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + reproductionmodes.or.f + (1|Site), family=binomial, data=mtdna_data)

#source/study as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + fertilizations.or.f + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + reproductionmodes.or.f + (1|Source), family=binomial, data=mtdna_data)

#MarkerName as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + (1|MarkerName), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + fertilizations.or.f + (1|MarkerName), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + reproductionmodes.or.f + (1|MarkerName), family=binomial, data=mtdna_data)

##More Mixed Models with multiple fixed variables (max length & fecundity constant for the following code)##

#spp as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + (1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|spp), family=binomial, data=mtdna_data)

#site as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|Site), family=binomial, data=mtdna_data)

#source/study as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|Source), family=binomial, data=mtdna_data)

#MarkerName as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + (1|MarkerName), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|MarkerName), family=binomial, data=mtdna_data)

##More Mixed Models with multiple fixed variables (He & fecundity constant for the following code)##

#spp as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 +  logtransform.maxlength.1 +(1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + fertilizations.or.f + (1|spp), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|spp), family=binomial, data=mtdna_data)

#site as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + logtransform.maxlength.1 + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + fertilizations.or.f + (1|Site), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|Site), family=binomial, data=mtdna_data)

#source/study as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + logtransform.maxlength.1 + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + fertilizations.or.f + (1|Source), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|Source), family=binomial, data=mtdna_data)

#MarkerName as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + logtransform.maxlength.1 + (1|MarkerName), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + fertilizations.or.f + (1|MarkerName), family=binomial, data=mtdna_data)

fit.bin <- glmer(mtdnas.or.f ~ logtransform.fecundity_mean.1 + reproductionmodes.or.f + (1|MarkerName), family=binomial, data=mtdna_data)


##ALL VARIABLES MODELED TOGETHER##

#spp as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + (1|spp), family=binomial, data=mtdna_data)

#site as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + (1|Site), family=binomial, data=mtdna_data, na.action = na.fail)

#source/study as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + (1|Source), family=binomial, data=mtdna_data)

#MarkerName as Random variable
fit.bin <- glmer(mtdnas.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + (1|MarkerName), family=binomial, data=mtdna_data)

################################################### msat data set ################################################### 

##Create variable to model success/failures
msat_data$msats.or.f <- round(msat_data$He*msat_data$n) & round((1-msat_data$He)*msat_data$n)  #create column of successes or failures

#Logtransform
for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.fecundity_mean.2 <- log10(msat_data$fecundity_mean)
}

for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.maxlength.2 <- log10(msat_data$maxlength)
}

#Add column for final fertilization

msat_data$final_fertilization [msat_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization
msat_data$final_fertilization [msat_data$fertilization =="external"]  <- "external"
msat_data$final_fertilization [msat_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"


#Create character to factor to integer columns
msat_data$fertilizations.or.f2 <- as.factor(msat_data$final_fertilization)
msat_data$fertilizations.or.f2 <- as.numeric(msat_data$fertilizations.or.f)

msat_data$reproductionmodes.or.f2 <- as.factor(msat_data$reproductionmode)
msat_data$reproductionmodes.or.f2 <- as.numeric(msat_data$reproductionmodes.or.f)

##Fit model for success/failures

#dredge? 
msat_data <- msat_data[!is.na(msat_data$Repeat), ]
msat_data <- msat_data[!is.na(msat_data$fertilizations.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$reproductionmodes.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$Repeat), ]
msat_data <- msat_data[!is.na(msat_data$logtransform.maxlength.2), ]
summary(mtdna_data)
model <- glm(formula = msats.or.f ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|spp) + (1|Site) + (1|Source), data = msat_data, na.action = na.fail)
model <- glmer(formula = msats.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + Repeat + CrossSpp + (1|spp) + (1|Site) + (1|Source) + (1|MarkerName), family=binomial, data = mtdna_data)
model <- glmer(formula = msats.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + Repeat + CrossSpp + (1|spp), family=binomial, data = mtdna_data)
model <- glmer(formula = msats.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + Repeat + CrossSpp + (1|Site), family=binomial, data = mtdna_data)
model <- glmer(formula = msats.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + Repeat + CrossSpp + (1|Source), family=binomial, data = mtdna_data, na.action = 'na.fail')
model <- glmer(formula = msats.or.f ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + Repeat + CrossSpp + (1|MarkerName), family=binomial, data = mtdna_data)


mtdna.spp <- dredge(model)
summary_mtdna.spp <-summary(mtdna.spp)

#spp as Random variable
fit.bin <- glmer(msats.or.f ~ He + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ logtransform.maxlength.2 + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ logtransform.fecundity_mean.2 + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ fertilizations.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ reproductionmodes.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ Repeat + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ CrossSpp + (1|spp), family=binomial, data=msat_data)

#site as Random variable
fit.bin <- glmer(msats.or.f ~ He + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ logtransform.maxlength.2 + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ logtransform.fecundity_mean.2 + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ fertilizations.or.f + (1|Site), family=binomial, data=msat_data) ####CONVERGENCE PROBLEM

fit.bin <- glmer(msats.or.f ~ reproductionmodes.or.f + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ Repeat + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ CrossSpp + (1|Site), family=binomial, data=msat_data)

#source/study as Random variable
fit.bin <- glmer(msats.or.f ~ He + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ logtransform.maxlength.2 + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ logtransform.fecundity_mean.2 + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ fertilizations.or.f + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ reproductionmodes.or.f + (1|Source), family=binomial, data=msat_data) #Downdated VtV is not positive definite

fit.bin <- glmer(msats.or.f ~ Repeat + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ CrossSpp + (1|Source), family=binomial, data=msat_data)

##More Mixed Models with multiple fixed variables##

#spp as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.fecundity_mean.2 + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + fertilizations.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + reproductionmodes.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + Repeat + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + CrossSpp + (1|spp), family=binomial, data=msat_data)

#site as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.fecundity_mean.2 + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + fertilizations.or.f + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + reproductionmodes.or.f + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + Repeat + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + CrossSpp + (1|Site), family=binomial, data=msat_data)

#source/study as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.fecundity_mean.2 + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + fertilizations.or.f + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + reproductionmodes.or.f + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + Repeat + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + CrossSpp + (1|Source), family=binomial, data=msat_data)

##More Mixed Models with multiple fixed variables (He & max length constant for the following code)##

#spp as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + fertilizations.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + reproductionmodes.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + Repeat + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + CrossSpp + (1|spp), family=binomial, data=msat_data)

#site as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + fertilizations.or.f + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + reproductionmodes.or.f + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + Repeat + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + CrossSpp + (1|Site), family=binomial, data=msat_data)

#source/study as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + fertilizations.or.f + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + reproductionmodes.or.f + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + Repeat + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + CrossSpp + (1|Source), family=binomial, data=msat_data)

##More Mixed Models with multiple fixed variables (He, max length, fecundity constant for the following code)##

#spp as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + reproductionmodes.or.f + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + Repeat + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + CrossSpp + (1|spp), family=binomial, data=msat_data)

#site as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + reproductionmodes.or.f + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + Repeat + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + CrossSpp + (1|Site), family=binomial, data=msat_data)

#source/study as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + reproductionmodes.or.f + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + Repeat + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + CrossSpp + (1|Source), family=binomial, data=msat_data)

##More Mixed Models with multiple fixed variables (He, max length, fecundity, fertilization constant for the following code)##

#spp as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f +  (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + Repeat + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + CrossSpp + (1|spp), family=binomial, data=msat_data)

#site as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f +  (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + Repeat + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + CrossSpp + (1|Site), family=binomial, data=msat_data)

#source/study as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f +  (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + Repeat + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + CrossSpp + (1|Source), family=binomial, data=msat_data)

##More Mixed Models with multiple fixed variables (He, max length, fecundity, fertilization, reproduction constant for the following code)##

#spp as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + Repeat + (1|spp), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + CrossSpp + (1|spp), family=binomial, data=msat_data)

#site as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + Repeat + (1|Site), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + CrossSpp + (1|Site), family=binomial, data=msat_data)

#source/study as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + Repeat + (1|Source), family=binomial, data=msat_data)

fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + CrossSpp + (1|Source), family=binomial, data=msat_data)

##All variables modeled together##

#spp as Random variable
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + Repeat + CrossSpp + (1|spp), family=binomial, data=msat_data)

#site as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + Repeat + CrossSpp + (1|Site), family=binomial, data=msat_data)

#source/study as Random variable 
fit.bin <- glmer(msats.or.f ~ He + logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f + reproductionmodes.or.f + Repeat + CrossSpp + (1|Source), family=binomial, data=msat_data)

######################################################################################################
#Examine residual plots to check assumptions
par(mfrow=c(1,2))
plot(residuals(fit.bin)~fitted(fit.bin),main="residuals v.s. Fitted")
qqnorm(residuals(fit.bin))

#Use parametric bootstrap 
nBoot <- 1000
lrStat <- rep(NA,nBoot)
ft.null <- glmer(s.or.f ~ 1 + (1|maxlength) ,family=binomial, data=msat_data) #null model
ft.alt <- glmer(s.or.f ~ spp + (1|maxlength) ,family=binomial, data=msat_data) #alternate model

lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat

for(iBoot in 1:nBoot)
{
  mtdna_data$s.or.fsim <- unlist(simulate(ft.null)) #resampled data
  tryCatch({#sometimes the glmer code doesn't converge
    
    bNull <- glmer(s.or.fsim ~ 1 + (1|maxlength) ,family=binomial, data=Estuaries)#null model
    bAlt <- glmer(s.or.fsim ~ source + (1|Estuary) ,family=binomial, data=Estuaries)#alternate model
    lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull) #resampled test stat
  },warning=function(war) {lrStat[iBoot]=NA},error=function(err){lrStat[iBoot]=NA}) #if code doesn't converge skip sim
}
mean(lrStat>lrObs,na.rm=T) #P-value for test of Estuary effect

#Model the counts of hydroids
fit.pois <- glmer(Hydroid ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries)

#Check the assumptions
par(mfrow=c(1,2))
plot(residuals(fit.pois)~fitted(fit.pois),main="Residuals vs. Fitted")
qqnorm(residuals(fit.pois))

#Parametric bootstrap to test for an effect of Modification
nBoot <- 1000
lrStat <- rep(NA,nBoot)
ft.null <- glmer(Hydroid ~ 1 + (1|Estuary) ,family=poisson, data=Estuaries) #null model
ft.alt <- glmer(Hydroid ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries) #alternate model

lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat
for(iBoot in 1:nBoot)
{
  Estuaries$HydroidSim <- unlist(simulate(ft.null)) #resampled data
  tryCatch({
    bNull <- glmer(HydroidSim ~ 1 + (1|Estuary)  ,family=poisson, data=Estuaries)#null model
    bAlt <- glmer(HydroidSim ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries) #alternate model
    lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull) #resampled test stat
  },warning=function(war) {lrStat[iBoot]=NA},error=function(err){lrStat[iBoot]=NA})  #if code doesn't converge skip sim#   lrStat[iBoot]
}
mean(lrStat>lrObs,na.rm=TRUE) #P-value for test of Estuary effect

#Model counts of bryozoan, Schizoporella errata
fit.pois2 <- glmer(Schizoporella.errata ~ Modification + (1|Estuary), family=poisson,  data=Estuaries)
par(mfrow=c(1,2))
plot(residuals(fit.pois)~fitted(fit.pois),main="residuals vs. Fitted")
qqnorm(residuals(fit.pois))

#Communicate results
fit.pois <- glmer(Hydroid ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries)
means <- fitted(fit.pois) #this will give the estimate at each data point
ModEst <- unique(Estuaries[c("Estuary", "Modification")])#find which Estuaries are modified
cols <- as.numeric(ModEst[order(ModEst[,1]),2])+3 #Assign colour by modification
boxplot(Hydroid~ Estuary,data=Estuaries,col=cols,xlab="Estuary",ylab="Count of hydroids")
legend("topleft", inset=.02,
       c("Modified","Pristine"), fill=unique(cols), horiz=TRUE, cex=0.8)

Est.means <- summarize(means,Estuaries$Estuary, mean)$means #extract means by Estuary
stripchart(Est.means~ sort(unique(Estuary)),data=Estuaries,pch=18,col="red",vertical = TRUE,add=TRUE) #plot means by estuary

#Examine residual plots to check assumptions
par(mfrow=c(1,2))
plot(residuals(fit.bin)~fitted(fit.bin),main="residuals v.s. Fitted")
qqnorm(residuals(fit.bin))

#Use parametric bootstrap 
nBoot <- 1000
lrStat <- rep(NA,nBoot)
ft.null <- glmer(mtdnas.or.f ~ 1 + (1|spp) ,family=binomial, data=mtdna_data) #null model
ft.alt <- glmer(mtdnas.or.f ~ He + (1|spp) ,family=binomial, data=mtdna_data) #alternate model

lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat

for(iBoot in 1:nBoot)
{
  mtdna_data$mtdnas.or.f2 <- unlist(simulate(ft.null)) #resampled data
  tryCatch({#sometimes the glmer code doesn't converge
    
    bNull <- glmer(mtdnas.or.f2 ~ 1 + (1|spp) ,family=binomial, data=mtdna_data)#null model
    bAlt <- glmer(mtdnas.or.f2 ~ He + (1|spp) ,family=binomial, data=mtdna_data)#alternate model
    lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull) #resampled test stat
  },warning=function(war) {lrStat[iBoot]=NA},error=function(err){lrStat[iBoot]=NA}) #if code doesn't converge skip sim
}
mean(lrStat>lrObs,na.rm=T) #P-value for test of Estuary effect

#Model the counts of hydroids
fit.pois <- glmer(mtdnas.or.f ~ He + (1|spp) ,family=poisson, data=mtdna_data)

#Check the assumptions
par(mfrow=c(1,2))
plot(residuals(fit.pois)~fitted(fit.pois),main="Residuals vs. Fitted")
qqnorm(residuals(fit.pois))

#Parametric bootstrap to test for an effect of Modification
nBoot <- 1000
lrStat <- rep(NA,nBoot)
ft.null <- glmer(mtdnas.or.f ~ 1 + (1|maxlength) ,family=poisson, data=mtdna_data) #null model
ft.alt <- glmer(mtdnas.or.f ~ spp + (1|maxlength) ,family=poisson, data=mtdna_data) #alternate model

lrObs <- 2*logLik(ft.alt) - 2*logLik(ft.null) #observed test stat
for(iBoot in 1:nBoot)
{
  mtdna_data$mtdnas.or.fx <- unlist(simulate(ft.null)) #resampled data
  tryCatch({
    bNull <- glmer(mtdnas.or.fx ~ 1 + (1|maxlength)  ,family=poisson, data=mtdna_data)#null model
    bAlt <- glmer(mtdnas.or.fx ~ spp + (1|maxlength) ,family=poisson, data=mtdna_data) #alternate model
    lrStat[iBoot] <- 2*logLik(bAlt) - 2*logLik(bNull) #resampled test stat
  },warning=function(war) {lrStat[iBoot]=NA},error=function(err){lrStat[iBoot]=NA})  #if code doesn't converge skip sim#   lrStat[iBoot]
}
mean(lrStat>lrObs,na.rm=TRUE) #P-value for test of Estuary effect

#Model counts of bryozoan, Schizoporella errata
fit.pois2 <- glmer(Schizoporella.errata ~ Modification + (1|Estuary), family=poisson,  data=Estuaries)
par(mfrow=c(1,2))
plot(residuals(fit.pois)~fitted(fit.pois),main="residuals vs. Fitted")
qqnorm(residuals(fit.pois))

#Communicate results
fit.pois <- glmer(Hydroid ~ Modification + (1|Estuary) ,family=poisson, data=Estuaries)
means <- fitted(fit.pois) #this will give the estimate at each data point
ModEst <- unique(Estuaries[c("Estuary", "Modification")])#find which Estuaries are modified
cols <- as.numeric(ModEst[order(ModEst[,1]),2])+3 #Assign colour by modification
boxplot(Hydroid~ Estuary,data=Estuaries,col=cols,xlab="Estuary",ylab="Count of hydroids")
legend("topleft", inset=.02,
       c("Modified","Pristine"), fill=unique(cols), horiz=TRUE, cex=0.8)

Est.means <- summarize(means,Estuaries$Estuary, mean)$means #extract means by Estuary
stripchart(Est.means~ sort(unique(Estuary)),data=Estuaries,pch=18,col="red",vertical = TRUE,add=TRUE) #plot means by estuary


