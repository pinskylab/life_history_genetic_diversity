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
library(effects)
library(sjPlot)

#read in data
mtdna_data <- read.csv("Datasets/new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("Datasets/new_msat_full_US_data.csv", stringsAsFactors = FALSE) #read in 

#Fixed Variables for mtDNA He & Pi: Body length/maxlength, fecundity mean, fertilization method, reproduction mode, Bp scale
#Random Variables for mtDNA He & Pi: Species, source
#Fixed Variables for msat He: Body length/maxlength, fecundity mean, fertilization method, reproduction mode, CrossSpp
#Random Variables for msat He: Species, source, ID
#in brood or similar structure --> internal

################################################### mtDNA data set ################################################### 

##### He #####

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
                                  (1|spp) + (1|Source), na.action = "na.fail", 
                                family=binomial, data = mtdna_data_He,
                                control = glmerControl(optimizer = "bobyqa")) 

mtdna_data_He_dredge <- dredge(binomial_He_full_model_mtDNA) #dredge model
View(mtdna_data_He_dredge) #see table
summary(binomial_He_full_model_mtDNA) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.mtDNAHE <- glmer(formula = cbind(success,failure) ~ fertilizations.or.f + reproductionmodes.or.f +
                                        (1|spp) + (1|Source), na.action = "na.fail", 
                                      family=binomial, data = mtdna_data_He,
                                      control = glmerControl(optimizer = "bobyqa")) 
mtdna_data_He_dredge_minimal <- dredge(topAIC.mtDNAHE) #dredge model
View(mtdna_data_He_dredge_minimal) #see table
summary(topAIC.mtDNAHE) #get SE, p-value, etc.for top AIC model

##### PI ##### 

##prep data & remove NA
mtdna_data_nona_fecunditymean <- subset(mtdna_data, mtdna_data$logtransform.fecundity_mean.1 != "NaN")
mtdna_data_nona_fecunditymean_bpscale <- subset(mtdna_data_nona_fecunditymean, mtdna_data_nona_fecunditymean$bp_scale != "NA")
mtdna_pi <- subset(mtdna_data_nona_fecunditymean_bpscale, mtdna_data_nona_fecunditymean_bpscale$logtransform.Pi != "NA") #remove any rows where mtdna pi wasn't calculated

##create full model for pi
Pi_full_model <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                        fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                        (1|spp) + (1|Source), na.action = "na.fail", 
                      data = mtdna_pi, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa")) #switch to bobyqa to get away from convergence issues AND switch to bp_scale AND make sure REML = FALSE

mtdna_pi_dredge <- dredge(Pi_full_model) #dredge model
View(mtdna_pi_dredge) #see table
summary(Pi_full_model) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.mtDNAPI <- lmer(formula = logtransform.Pi ~ fertilizations.or.f + 
                         (1|spp) + (1|Source), na.action = "na.fail", 
                       data = mtdna_pi, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa")) #use variables from top AIC model
mtdna_data_Pi_dredge_minimal <- dredge(topAIC.mtDNAPI) #dredge model
View(mtdna_data_Pi_dredge_minimal)
summary(topAIC.mtDNAPI) #get SE, p-value, etc.for top AIC model
summary(mtdna_data_Pi_dredge_minimal)

################################################### msat data set ################################################### 

##Create ID column
msat_data$ID <- c(1:3145) 

##Logtransform
for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.fecundity_mean.2 <- log10(msat_data$fecundity_mean)
}

for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.maxlength.2 <- log10(msat_data$maxlength)
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
msat_data <- msat_data %>% drop_na(He, logtransform.maxlength.2, logtransform.fecundity_mean.2,fertilizations.or.f2,reproductionmodes.or.f2)

##Create full model for HE
binomial_He_full_model_msat <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 +
                                  fertilizations.or.f2 + reproductionmodes.or.f2 + CrossSpp +
                                  (1|spp) + (1|Source) + (1|ID), na.action = "na.fail", 
                                family=binomial, data = msat_data,
                                control = glmerControl(optimizer = "bobyqa")) #have Marial switch to bobyqa to get away from convergence issues AND switch to bp_scale

msat_dataHe <- dredge(binomial_He_full_model_msat) #dredge model
View(msat_dataHe) #see table
summary(binomial_He_full_model_msat) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.msatHE <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + fertilizations.or.f2 + CrossSpp +
                         (1|spp) + (1|Source) + (1|ID), na.action = "na.fail", 
                       family=binomial, data = msat_data,
                       control = glmerControl(optimizer = "bobyqa"))
msat_dataHe_minimal <- dredge(topAIC.msatHE) #dredge model
View(topAIC.msatHE)
summary(topAIC.msatHE) #get SE, p-value, etc.for top AIC model

##find RVI 
msatRVImax <- sum(msat_dataHe[complete.cases(msat_dataHe$logtransform.maxlength.2), "weight"])
msatRVIfec <- sum(msat_dataHe[complete.cases(msat_dataHe$logtransform.fecundity_mean.2), "weight"])
msatRVIfert <- sum(msat_dataHe[complete.cases(msat_dataHe$fertilizations.or.f2), "weight"])
msatRVIrepro <- sum(msat_dataHe[complete.cases(msat_dataHe$reproductionmodes.or.f2), "weight"])
msatRVIcross <- sum(msat_dataHe[complete.cases(msat_dataHe$CrossSpp), "weight"])

################################################### Graphing results ################################################### 

######### Mixed model figures: Fecundity #######

#### Fecundity & Pi #### 

Pi_fecund <- plot_model(Pi_full_model, type = "pred", terms = "logtransform.fecundity_mean.1")

#pull out marginal effects dataframe
Pi_fecund_data <- as.data.frame(Pi_fecund$data)

#unlog fecundity
Pi_fecund_data$unlog_fecund <- 10^(Pi_fecund_data$x)

#unlog pi
Pi_fecund_data$unlog_pi <- 10^(Pi_fecund_data$predicted)
Pi_fecund_data$unlog_conf.low <- 10^(Pi_fecund_data$conf.low)
Pi_fecund_data$unlog_conf.high <- 10^(Pi_fecund_data$conf.high)

### Plot fecundity ###

mtdna_pi_plot_both <- ggplot() +
  geom_line(data = Pi_fecund_data, 
            aes(x = unlog_fecund, y = unlog_pi), linewidth = 3) + 
  geom_ribbon(data = Pi_fecund_data, 
              aes(x = unlog_fecund, ymin = unlog_conf.low, ymax = unlog_conf.high), alpha = 0.1) + 
  xlab("Fecundity") + ylab("Mitochondrial Pi") +
  geom_rug(data = mtdna_pi, aes(x = fecundity_mean, y = NULL), inherit.aes = F, 
           alpha = 0.6, length = unit(20,"pt"))
  
mtdna_pi_fecundity_plot_annotated_both <- mtdna_pi_plot_both + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), axis.title = element_text(size = 11), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 25, color = "black"), 
        axis.line = element_line(color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        axis.title.x = element_text( size=25, margin = margin(t = 20)),
        axis.title.y = element_text(size=25, margin = margin(r = 20)),
        text = element_text(size = 28)) +
  xlim(c(0, 1.03e7))
mtdna_pi_fecundity_plot_annotated_both

#### Fecundity & mtDNA He #### 

mtDNAHe_fecund <- plot_model(binomial_He_full_model_mtDNA, type = "pred", terms = "logtransform.fecundity_mean.1 [all]")

#pull out marginal effects dataframe
mtDNAHe_fecund_data <- as.data.frame(mtDNAHe_fecund$data)

#unlog fecundity
mtDNAHe_fecund_data$unlog_fecund <- 10^(mtDNAHe_fecund_data$x)

### Plot fecundity ###

mtDNA_He_plot_both <- ggplot() + 
  geom_line(data = mtDNAHe_fecund_data, 
            aes(x = unlog_fecund, y = predicted), linewidth = 3) + 
  geom_ribbon(data = mtDNAHe_fecund_data, 
              aes(x = unlog_fecund, ymin = conf.low, ymax = conf.high), alpha = 0.1) + 
  xlab("Fecundity") + ylab("Mitochondrial Hd") +
  geom_rug(data = mtdna_data_He, aes(x = fecundity_mean, y = NULL), inherit.aes = F, 
           alpha = 0.6, length = unit(20,"pt"))

mtDNA_He_fecundity_plot_annotated_both <- mtDNA_He_plot_both + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), axis.title = element_text(size = 11), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 25, color = "black"), 
        axis.line = element_line(color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        axis.title.x = element_text(size=25, margin = margin(t = 20)),
        axis.title.y = element_text(size=25, margin = margin(r = 20)),
        text = element_text(size = 28)) +
  ylim(c(0, 1.05))
mtDNA_He_fecundity_plot_annotated_both

#### Fecundity & msat He #### 

msatHe_fecund <- plot_model(binomial_He_full_model_msat, type = "pred", terms = "logtransform.fecundity_mean.2 [all]")

#pull out marginal effects dataframe
msatHe_fecund_data <- as.data.frame(msatHe_fecund$data)

#unlog fecundity
msatHe_fecund_data$unlog_fecund <- 10^(msatHe_fecund_data$x)

### Plot fecundity ###

msat_He_plot_both <- ggplot() + 
  geom_line(data = msatHe_fecund_data, 
            aes(x = unlog_fecund, y = predicted), linewidth = 3) + 
  geom_ribbon(data = msatHe_fecund_data, 
              aes(x = unlog_fecund, ymin = conf.low, ymax = conf.high), alpha = 0.1) + 
  xlab("Fecundity") + ylab("Nuclear He") +
  geom_rug(data = msat_data, aes(x = fecundity_mean, y = NULL), inherit.aes = F, 
           alpha = 0.6, length = unit(20,"pt"))

msat_He_fecundity_plot_annotated_both <- msat_He_plot_both + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), axis.title = element_text(size = 11), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 25, color = "black"), 
        axis.line = element_line(color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        axis.title.x = element_text(size=25, margin = margin(t = 20)),
        axis.title.y = element_text(size=25, margin = margin(r = 20)),
        text = element_text(size = 28))
msat_He_fecundity_plot_annotated_both

######### Mixed model figures: Max Length #######

#### Length & Pi #### 

Pi_length <- plot_model(Pi_full_model, type = "pred", terms = "logtransform.maxlength.1")

#pull out marginal effects dataframe
Pi_length_data <- as.data.frame(Pi_length$data)

#unlog length
Pi_length_data$unlog_length <- 10^(Pi_length_data$x)

#unlog pi
Pi_length_data$unlog_pi <- 10^(Pi_length_data$predicted)
Pi_length_data$unlog_conf.low <- 10^(Pi_length_data$conf.low)
Pi_length_data$unlog_conf.high <- 10^(Pi_length_data$conf.high)

### Plot length ###

mtdna_pi_plot_both_length <- ggplot() + 
  geom_line(data = Pi_length_data, 
            aes(x = unlog_length, y = unlog_pi), linewidth = 3) + 
  geom_ribbon(data = Pi_length_data, 
              aes(x = unlog_length, ymin = unlog_conf.low, ymax = unlog_conf.high), alpha = 0.1) + 
  xlab("Maximum Length") + ylab("Mitochondrial Pi") +
  geom_rug(data = mtdna_pi, aes(x = maxlength, y = NULL), inherit.aes = F, 
           alpha = 0.6, length = unit(20,"pt"))

mtdna_pi_length_plot_annotated_both <- mtdna_pi_plot_both_length + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), axis.title = element_text(size = 11), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 25, color = "black"), 
        axis.line = element_line(color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        axis.title.x = element_text(size=25, margin = margin(t = 20)),
        axis.title.y = element_text(size=25, margin = margin(r = 20)),
        text = element_text(size = 28))
mtdna_pi_length_plot_annotated_both

#### Length & mtDNA He #### 

mtDNAHe_length <- plot_model(binomial_He_full_model_mtDNA, type = "pred", terms = "logtransform.maxlength.1 [all]")

#pull out marginal effects dataframe
mtDNAHe_length_data <- as.data.frame(mtDNAHe_length$data)

#unlog length
mtDNAHe_length_data$unlog_length <- 10^(mtDNAHe_length_data$x)

### Plot length ###

mtDNA_He_plot_both_length <- ggplot() + 
  geom_line(data = mtDNAHe_length_data, 
            aes(x = unlog_length, y = predicted), linewidth = 3) + 
  geom_ribbon(data = mtDNAHe_length_data, 
              aes(x = unlog_length, ymin = conf.low, ymax = conf.high), alpha = 0.1) + 
  xlab("Maximum Length") + ylab("Mitochondrial Hd") +
  geom_rug(data = mtdna_data_He, aes(x = maxlength, y = NULL), inherit.aes = F, 
           alpha = 0.6, length = unit(20,"pt"))

mtDNA_He_length_plot_annotated_both <- mtDNA_He_plot_both_length + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), axis.title = element_text(size = 11), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 25, color = "black"), 
        axis.line = element_line(color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        axis.title.x = element_text(size=25, margin = margin(t = 20)),
        axis.title.y = element_text(size=25, margin = margin(r = 20)),
        text = element_text(size = 28))
mtDNA_He_length_plot_annotated_both

#### Length & msat He #### 

msatHe_length <- plot_model(binomial_He_full_model_msat, type = "pred", terms = "logtransform.maxlength.2 [all]")

#pull out marginal effects dataframe
msatHe_length_data <- as.data.frame(msatHe_length$data)

#unlog length
msatHe_length_data$unlog_length <- 10^(msatHe_length_data$x)

### Plot length ###

msat_He_plot_both_length <- ggplot() + 
  geom_line(data = msatHe_length_data, 
            aes(x = unlog_length, y = predicted), linewidth = 3) + 
  geom_ribbon(data = msatHe_length_data, 
              aes(x = unlog_length, ymin = conf.low, ymax = conf.high), alpha = 0.1) + 
  xlab("Maximum Length") + ylab("Nuclear He") +
  geom_rug(data = msat_data, aes(x = maxlength, y = NULL), inherit.aes = F, 
           alpha = 0.6, length = unit(20,"pt"))

msat_He_length_plot_annotated_both <- msat_He_plot_both_length + theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1), axis.title = element_text(size = 11), 
        axis.ticks = element_line(color = "black", linewidth = 1), 
        axis.text = element_text(size = 25, color = "black"), 
        axis.line = element_line(color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        axis.title.x = element_text(size=25, margin = margin(t = 20)),
        axis.title.y = element_text(size=25, margin = margin(r = 20)),
        text = element_text(size = 28)) +
  xlim(c(0, 500.5))
msat_He_length_plot_annotated_both


################################################### Exclude Red List ################################################### 

################################################### mtDNA data set ################################################### 

##### He #####

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

##Read in and exclude spps that are vulernable, endangered, or critically endangered
IUCN <- read.csv("Datasets/IUCN_status.csv", stringsAsFactors = FALSE)

mtdna_data_IUCN <- merge(mtdna_data, IUCN, by=c('spp'), all =F)

mtdna_data_IUCN <- mtdna_data_IUCN %>% filter(!IUCN_status %in% c("Endangered", "Vulnerable",  "Critically endangered")) 

##Remove NA from variable columns
mtdna_data_IUCN_He <- mtdna_data_IUCN %>% drop_na(He, logtransform.maxlength.1, logtransform.fecundity_mean.1,fertilizations.or.f,reproductionmodes.or.f,
                                        bp_scale)

##Create full model for HE
binomial_He_full_model_mtDNA_IUCN <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                                        fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                                        (1|spp) + (1|Source), na.action = "na.fail", 
                                      family=binomial, data = mtdna_data_IUCN_He,
                                      control = glmerControl(optimizer = "bobyqa")) 

mtdna_data_He_dredge_IUCN <- dredge(binomial_He_full_model_mtDNA_IUCN) #dredge model
View(mtdna_data_He_dredge_IUCN) #see table
summary(binomial_He_full_model_mtDNA_IUCN) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.mtDNAHE_IUCN <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.1 + fertilizations.or.f + 
                          (1|spp) + (1|Source), na.action = "na.fail", 
                        family=binomial, data = mtdna_data_IUCN_He,
                        control = glmerControl(optimizer = "bobyqa")) 
mtdna_data_He_dredge_minima_IUCN <- dredge(topAIC.mtDNAHE_IUCN) #dredge model
View(mtdna_data_He_dredge_minima_IUCN) #see table
summary(topAIC.mtDNAHE_IUCN) #get SE, p-value, etc.for top AIC model

##### PI ##### 

##prep data & remove NA
mtdna_data_nona_fecunditymean_IUCN <- subset(mtdna_data_IUCN, mtdna_data_IUCN$logtransform.fecundity_mean.1 != "NaN")
mtdna_data_nona_fecunditymean_bpscale_IUCN <- subset(mtdna_data_nona_fecunditymean_IUCN, mtdna_data_nona_fecunditymean_IUCN$bp_scale != "NA")
mtdna_pi_IUCN <- subset(mtdna_data_nona_fecunditymean_bpscale_IUCN, mtdna_data_nona_fecunditymean_bpscale_IUCN$logtransform.Pi != "NA") #remove any rows where mtdna pi wasn't calculated

##create full model for pi
Pi_full_model_IUCN <- lmer(formula = logtransform.Pi ~ logtransform.maxlength.1 + logtransform.fecundity_mean.1 + 
                        fertilizations.or.f + reproductionmodes.or.f + bp_scale + 
                        (1|spp) + (1|Source), na.action = "na.fail", 
                      data = mtdna_pi_IUCN, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa")) #switch to bobyqa to get away from convergence issues 

mtdna_pi_dredge_IUCN <- dredge(Pi_full_model_IUCN) #dredge model
View(mtdna_pi_dredge_IUCN) #see table
summary(Pi_full_model_IUCN) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.mtDNAPI_IUCN <- lmer(formula = logtransform.Pi +
                         (1|spp) + (1|Source), na.action = "na.fail", 
                       data = mtdna_pi_IUCN, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa")) #use variables from top AIC model
mtdna_data_Pi_dredge_minimal_IUCN <- dredge(topAIC.mtDNAPI_IUCN) #dredge model
View(mtdna_data_Pi_dredge_minimal_IUCN)
summary(topAIC.mtDNAPI_IUCN) #get SE, p-value, etc.for top AIC model

################################################### msat data set ################################################### 

##Create ID column
msat_data$ID <- c(1:3145) 

##Logtransform
for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.fecundity_mean.2 <- log10(msat_data$fecundity_mean)
}

for (i in 1:nrow(msat_data)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_data$logtransform.maxlength.2 <- log10(msat_data$maxlength)
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

##Read in and exclude spps that are vulernable, endangered, or critically endangered
msat_data_IUCN <- merge(msat_data, IUCN, by=c('spp'), all =F)

msat_data_IUCN <- msat_data_IUCN %>% filter(!IUCN_status %in% c("Endangered", "Vulnerable",  "Critically endangered")) 

##Remove NA from variable columns
msat_data_IUCN <- msat_data_IUCN %>% drop_na(He, logtransform.maxlength.2, logtransform.fecundity_mean.2,fertilizations.or.f2,reproductionmodes.or.f2)

##Create full model for HE
binomial_He_full_model_msat_IUCN <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + logtransform.fecundity_mean.2 +
                                       fertilizations.or.f2 + reproductionmodes.or.f2 + CrossSpp +
                                       (1|spp) + (1|Source) + (1|ID), na.action = "na.fail", 
                                     family=binomial, data = msat_data_IUCN,
                                     control = glmerControl(optimizer = "bobyqa")) 

msat_dataHe_IUCN <- dredge(binomial_He_full_model_msat_IUCN) #dredge model
View(msat_dataHe_IUCN) #see table
summary(binomial_He_full_model_msat_IUCN) #get SE, p-value, etc.

##find minimal model (top AIC model)
topAIC.msatHE_IUCN <- glmer(formula = cbind(success,failure) ~ logtransform.maxlength.2 + fertilizations.or.f2 + CrossSpp +
                         (1|spp) + (1|Source) + (1|ID), na.action = "na.fail", 
                       family=binomial, data = msat_data_IUCN,
                       control = glmerControl(optimizer = "bobyqa"))
msat_dataHe_minimal_IUCN <- dredge(topAIC.msatHE_IUCN ) #dredge model
View(msat_dataHe_minimal_IUCN)
summary(topAIC.msatHE_IUCN) #get SE, p-value, etc.for top AIC model
