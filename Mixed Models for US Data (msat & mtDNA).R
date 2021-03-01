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

#read in data
mtdna_data <- read.csv("new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("new_msat_full_US_data.csv", stringsAsFactors = FALSE) #read in 

#Fixed Variables: Crossspp, repeats, fecundity, body length/maxlength, reproduction mode, fertilization method
#Random Variables: Species, site, source/study, marker name
#in brood or similar structure --> internal

############### mtDNA data set ############### 

##Create variable to model success/failures
# REMOVE: mtdna_data$mtdnas.or.f <- round(mtdna_data$He*mtdna_data$n) & round((1-mtdna_data$He)*mtdna_data$n)  #create column of successes or failures
mtdna_data$success <- round(mtdna_data$He*mtdna_data$n)
mtdna_data$failure <- round((1-mtdna_data$He)*mtdna_data$n)

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

for (i in 1:nrow(mtdna_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_data$logtransform.Pi <- log10(mtdna_data$Pi)
}

##Fit model for success/failures##

#dredge for He 
mtdna_data <- mtdna_data[!is.na(mtdna_data$bp), ]
# REMOVE: mtdna_data <- mtdna_data[!is.na(mtdna_data$mtdnas.or.f), ]
mtdna_data <- mtdna_data[!is.na(mtdna_data$success), ]
mtdna_data <- mtdna_data[!is.na(mtdna_data$failure), ]
mtdna_data <- mtdna_data[!grepl('NaN',mtdna_data$logtransform.fecundity_mean.1),]
mtdna_data <- mtdna_data[!is.na(mtdna_data$logtransform.Pi), ]
mtdna_data$fecundity_mean <- NULL
mtdna_data$bp <- NULL
summary(mtdna_data)
model <- glm(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1, family=binomial, data = mtdna_data,na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp) + (1|Source)+ (1|MarkerName), family=binomial, data = mtdna_data,  na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp), family=binomial, data = mtdna_data, na.action = 'na.fail')
# REMOVE: model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site), family=binomial, data = mtdna_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Source), family=binomial, data = mtdna_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName), family=binomial, data = mtdna_data,na.action = 'na.fail')

#combo of random variable for He 
# REMOVE: model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName) + (1|Site), family=binomial, data = mtdna_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName) + (1|Source), family=binomial, data = mtdna_data, na.action = 'na.fail')
# REMOVE: model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName) + (1|Site) + (1|Source), family=binomial, data = mtdna_data, na.action = 'na.fail')


mtdna.spp <- dredge(model)
View(mtdna.spp) #to get a table that can be copy and pasted to Excel

#Check top 4 models individually to see if there were convergence problems for combo of all random variables
fit.bin <- glmer(cbind(success,failure) ~ fertilizations.or.f + reproductionmodes.or.f + (1|spp) + (1|Source)+ (1|MarkerName), family=binomial, data=mtdna_data)
fit.bin <- glmer(cbind(success,failure) ~ fertilizations.or.f + (1|spp) + (1|Source)+ (1|MarkerName), family=binomial, data=mtdna_data)
fit.bin <- glmer(cbind(success,failure) ~ fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp) + (1|Source)+ (1|MarkerName), family=binomial, data=mtdna_data)
fit.bin <- glmer(cbind(success,failure) ~ logtransform.maxlength.1 + fertilizations.or.f + (1|spp) + (1|Source)+ (1|MarkerName), family=binomial, data=mtdna_data)

#Find SE & P-Value for Fixed, SD for Random
sum <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName) + (1|Site), family=binomial, data = mtdna_data, na.action = 'na.fail')
summary(sum)

# extract coefficients
coefs <- data.frame(coef(summary(sum)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#dredge for pi
model <- glm(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1, data = mtdna_data, na.action = 'na.fail')
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp) + (1|Source)+ (1|MarkerName), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|spp), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
# REMOVE: model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Source), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|MarkerName), data = mtdna_data, na.action = 'na.fail',REML=FALSE)

#combo of random variable for pi 
# REMOVE: model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Site) + (1|Source), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + fertilizations.or.f + reproductionmodes.or.f + logtransform.bp.1 + (1|Source) + (1|MarkerName), data = mtdna_data, na.action = 'na.fail', REML=FALSE)
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

#try to fix convergence
length(getME(model,'theta'))
length(fixef(model))


tt <- getME(model,"theta")
ll <- getME(model,"lower")
min(tt[ll==0])

derivs1 <- model@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))

max(pmin(abs(sc_grad1),abs(derivs1$gradient)))


ss <- getME(model,c("theta","fixef"))
m2 <- update(model,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
m3 <- update(model,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

################################################### msat data set ################################################### 

##Create variable to model success/failures
# REMOVE: msat_data$msats.or.f <- round(msat_data$He*msat_data$n) & round((1-msat_data$He)*msat_data$n)  #create column of successes or failures
msat_data$success <- round(msat_data$He*msat_data$n)
msat_data$failure <- round((1-msat_data$He)*msat_data$n)

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

#dredge? 
msat_data <- msat_data[!is.na(msat_data$Repeat), ]
msat_data <- msat_data[!is.na(msat_data$fertilizations.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$reproductionmodes.or.f2), ]
msat_data <- msat_data[!is.na(msat_data$logtransform.maxlength.2), ]
msat_data <- msat_data[!is.na(msat_data$logtransform.fecundity_mean.2), ]
summary(msat_data)

model <- glm(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp, data = msat_data, family=binomial, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|spp) + (1|Source)+ (1|PrimerNote), family=binomial, data = msat_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|spp), family=binomial, data = msat_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|PrimerNote), family=binomial, data = msat_data, na.action = 'na.fail')
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|Source), family=binomial, data = msat_data, na.action = 'na.fail')

#combo of random variable for He
model <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 + fertilizations.or.f2 + reproductionmodes.or.f2 + Repeat + CrossSpp + (1|Site) + (1|Source) , family=binomial, data = msat_data, na.action = 'na.fail')


msat.spp <- dredge(model)
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
