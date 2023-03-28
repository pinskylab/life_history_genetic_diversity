################################################### Script for dataset ########################################################

#look at dataset and see different characteristics (ex: how many hermaphrodites in study)

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
mtdna_data <- read.csv("new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("new_msat_full_US_data.csv", stringsAsFactors = FALSE) #read in 

#Fixed Variables for mtDNA He & Pi: Body length/maxlength, fecundity mean, fertilization method, reproduction mode, Bp scale
#Random Variables for mtDNA He & Pi: Species, source
#Fixed Variables for msat He: Body length/maxlength, fecundity mean, fertilization method, reproduction mode, CrossSpp
#Random Variables for msat He: Species, source, ID
#in brood or similar structure --> internal

##########################################################################################################################################

############# MTDNA JUST LOOKING @ HD #############

#Number of Hd in mtDNA
mtDNA_He <- mtdna_data %>% drop_na(He) #138 obs.

#Number of repro in mtDNA Hd
mtDNA_He_repro <- mtDNA_He %>% drop_na(reproductionmode) 
table (mtDNA_He_repro$reproductionmode) #dioecism = 115 | protogyny = 17

#Number of fert in mtDNA Hd
mtDNA_He_fert <- mtDNA_He %>% drop_na(fertilization)
table (mtDNA_He_fert$fertilization) #external = 96 | internal (including in brood) = 36

#Number of sources
unique(mtDNA_He$Source) #34 sources

#Number of spps
unique(mtDNA_He$spp) #35 spps

###Fecundity
#Range
mtDNA_He[mtDNA_He == "NaN"] <- NA
range_fec_mtdnaHe <- max(mtDNA_He$fecundity_mean, na.rm =TRUE) - min(mtDNA_He$fecundity_mean, na.rm =TRUE) #5999999 (max = 6e+06, min = 1)

#Mode
tail(names(sort(table(mtDNA_He$fecundity_mean))), 1) #179

#Mean 
mean(mtDNA_He$fecundity_mean, na.rm =TRUE) #394261.5

###Max length (cm)
#Range
range_maxlength_mtdnaHe <- max(mtDNA_He$maxlength, na.rm =TRUE) - min(mtDNA_He$maxlength, na.rm =TRUE) #476.3 (max = 430, min = 5.7)

#Mode
tail(names(sort(table(mtDNA_He$maxlength))), 1) #91

#Mean
mean(mtDNA_He$maxlength) #94.09924

#Get # of genuses
mtDNA_He_genus <- str_split_fixed(mtDNA_He$spp, " ", 2)
mtDNA_He_genus <- data.frame(mtDNA_He_genus); colnames(mtDNA_He_genus) <- c("Genus", "Species")
unique(mtDNA_He_genus$Genus) #29 genuses

############# MTDNA JUST LOOKING @ PI #############

#Number of repro in mtDNA Pi
mtDNA_Pi_repro <- mtdna_data %>% drop_na(reproductionmode) 
table (mtDNA_Pi_repro$reproductionmode) #dioecism = 134 | protogyny = 17

#Number of fert in mtDNA Pi
mtDNA_Pi_fert <- mtDNA_He %>% drop_na(fertilization)
table (mtDNA_Pi_fert$fertilization) #external = 96 | internal (including in brood) = 36

#Number of sources
unique(mtdna_data$Source) #40 sources

#Number of spps
unique(mtdna_data$spp) #40 spps

###Fecundity
#Range
mtdna_data[mtdna_data == "NaN"] <- NA
range_fec_mtdnaPi <- max(mtdna_data$fecundity_mean, na.rm =TRUE) - min(mtdna_data$fecundity_mean, na.rm =TRUE) #5999999(max = 6e+06, min = 1)

#Mode
tail(names(sort(table(mtdna_data$fecundity_mean))), 1) #179

#Mean 
mean(mtdna_data$fecundity_mean, na.rm =TRUE) #394261.5

###Max length (cm)
#Range
range_maxlength_mtdnaPi <- max(mtdna_data$maxlength, na.rm =TRUE) - min(mtdna_data$maxlength, na.rm =TRUE) #445 (max = 455, min = 10)

#Mode
tail(names(sort(table(mtdna_data$maxlength))), 1) #91

#Mean 
mean(mtdna_data$maxlength) #94.09924

#Get # of genuses
mtDNA_Pi_genus <- str_split_fixed(mtdna_data$spp, " ", 2)
mtDNA_Pi_genus <- data.frame(mtDNA_Pi_genus); colnames(mtDNA_Pi_genus) <- c("Genus", "Species")
unique(mtDNA_Pi_genus$Genus) #32 genuses

############# MSAT #############

#Number of repro in msat
msat_repro <- msat_data %>% drop_na(reproductionmode) 
table (msat_repro$reproductionmode) #dioecism = 2786 | protogyny = 248

#Number of fert in msat
msat_fert <- msat_data %>% drop_na(fertilization)
table (msat_fert$fertilization) #external = 2026 | internal (including in brood) = 1014

#Number of sources
unique(msat_data$Source) #40 sources

#Number of spps
unique(msat_data$spp) #40 spps

#Number of crossspp
table(msat_data$CrossSpp) #0 = 1593 | 0.25 = 80 | 0.6 = 96 | 1 = 1376

###Fecundity
#Range
msat_data[msat_data == "NaN"] <- NA
range_fec_msat <- max(msat_data$fecundity_mean, na.rm =TRUE) - min(msat_data$fecundity_mean, na.rm =TRUE) #4655806 (max = 4655806, min = 1)

#Mode
tail(names(sort(table(msat_data$fecundity_mean))), 1) #4655806.5

#Mean 
mean(msat_data$fecundity_mean, na.rm =TRUE) #1385535

###Max length (cm)
#Range
range_maxlength_msat <- max(msat_data$maxlength, na.rm =TRUE) - min(msat_data$maxlength, na.rm =TRUE) #476.3 (max = 482, min = 5.7)

#Mode
tail(names(sort(table(msat_data$maxlength))), 1) #100

#Mean 
mean(msat_data$maxlength) #94.09924

#Get # of genuses
msat_genus <- str_split_fixed(msat_data$spp, " ", 2)
msat_genus <- data.frame(msat_genus); colnames(msat_genus) <- c("Genus", "Species")
unique(msat_genus$Genus) #43 genuses

############# MSAT & MTDNA (full mtdna dataset) #############

#merge
mtdna_mergeready <- mtdna_data %>% select(spp, Source)
mtdna_mergeready <- mtdna_mergeready[!duplicated(mtdna_mergeready),]

msat_mergeready <- msat_data %>% select(spp, Source)
msat_mergeready <- msat_mergeready[!duplicated(msat_mergeready),]

#Common spps& Sources
msat_mergeready <-msat_mergeready %>% rename(Source_msat = Source) #rename rows to make it easier to merge
mtdna_mergeready <-mtdna_mergeready %>% rename(Source_mtdna = Source) #rename rows to make it easier to merge

compare_both <- merge(msat_mergeready, mtdna_mergeready) #merge

compare_both$Source_mtdna[duplicated(compare_both$Source_mtdna)] <- NA #remove duplicate species to find how many common
compare_both$Source_msat[duplicated(compare_both$Source_msat)] <- NA #remove duplicate species to find how many common

#15 common spps and 8 common sources



