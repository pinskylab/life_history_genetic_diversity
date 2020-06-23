################################################### Script for Pulling Life History Traits  #######################################################

#pulls life history traits for both mtdna & msat species

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(rfishbase)
library(tidyverse)

#read in data
msat_spp <- read.csv("fbdat_msat.csv", stringsAsFactors = FALSE) #read in msat species list
mtdna_spp <- read.csv("fbdat_mtdna.csv", stringsAsFactors = FALSE) #read in mtdna species list

#NOTE: Branchiostoma spps in mtdna NOT in fishbase -- will need to pull out before querying FB for data

##########################################################################################################################################

######## Add life history traits for msat ########

######## Add body size ########

#add body size from FB
msat_spp$maxlength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get length data
  cat(paste(i, " ", sep = ''))
  msat_spp$maxlength[i] <- as.numeric(species(msat_spp$fbsci[i], fields = 'Length')$Length)
  if(is.na(msat_spp$maxlength[i])) { #if male maxlength is NA, use female maxlength
    msat_spp$maxlength[i] <- as.numeric(species(msat_spp$fbsci[i], fields = 'LengthFemale')$LengthFemale)
  }
}

summary(msat_spp)

#Ex of code without if statement
#for(i in 1:nrow(spps2)) { #get code that specifies exact species on FB
#  cat(paste(i, " ", sep = ''))
#  spps2$SpecCode[i] <- as.numeric(species(specieslist = spps2$fbsci[i], fields = 'SpecCode')$SpecCode)
#}