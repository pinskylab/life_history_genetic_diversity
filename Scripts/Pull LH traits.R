################################################### Script for Pulling Life History Traits  #######################################################

#pulls life history traits for both mtdna & msat species

##########################################################################################################################################

####### Set-up #######

remove(list = ls())

#load libraries
library(rfishbase)
library(tidyverse)
library(dplyr)
library(sjmisc)

#read in data
msat_spp <- read.csv("Datasets/msat_papers_USA_new.csv", stringsAsFactors = FALSE) #read in msat species list
mtdna_spp_new <- read.csv("Datasets/mtdna_papers_USA_new.csv", stringsAsFactors = FALSE) #read in mtdna species list

#NOTE: Branchiostoma spps in mtdna NOT in fishbase -- will need to pull out before querying FB for data

##########################################################################################################################################

######## Add life history traits for msat ########

######## Add body size ########

#add body size from FB

msat_spp$maxlength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get length data
  cat(paste(i, " ", sep = ''))
  msat_spp$maxlength[i] <- as.numeric(species(msat_spp$spp[i], fields = 'Length'))
  if(is.na(msat_spp$maxlength[i])) { #if male maxlength is NA, use female maxlength
    msat_spp$maxlength[i] <- as.numeric(species(msat_spp$spp[i], fields = 'LengthFemale'))
   }
}

#add fecundity from FB
msat_spp$fecundity <- NA #create column to fill in for fecundity
msat_spp$fecundity_mean <- NA #calculate mean fecundity (when vector is provided)

for (i in 1:nrow(msat_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  msat_spp$fecundity[i] <- fecundity(species=(msat_spp$spp[i]), field='FecundityMean')
  if(sjmisc::is_empty(unlist(msat_spp$fecundity[i]), all.na.empty = TRUE)) { #if fecundity mean is NA, use fecundity minimum
    msat_spp$fecundity[i] <- fecundity(species=msat_spp$spp[i], field= 'FecundityMin')
    if(sjmisc::is_empty(unlist(msat_spp$fecundity[i]), all.na.empty = TRUE)) { #if fecundity minimum is NA, use fecundity maximum
      msat_spp$fecundity[i] <- fecundity(species=(msat_spp$spp[i]), fields = 'FecundityMax')
    }
  }
  msat_spp$fecundity_mean[i] <- mean(unlist(msat_spp$fecundity[i]), na.rm = TRUE) #calculate mean fecundity
}

  msat_spp$fecundity <- as.integer(as.list(msat_spp$fecundity))
  msat_spp$fecundity_mean <- as.integer(msat_spp$fecundity_mean)
  
#add fertilization from FB
msat_spp$fertilization <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get fertilization data
  cat(paste(i, " ", sep = ''))
  msat_spp$fertilization[i] <- reproduction(species=(msat_spp$spp[i]), field='Fertilization')
}

#add reproduction mode from FB
msat_spp$reproductionmode <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get reproduction mode data
  cat(paste(i, " ", sep = ''))
  msat_spp$reproductionmode[i] <- reproduction(species=(msat_spp$spp[i]), field='ReproMode')
}

summary(msat_spp) #see results

##########################################################################################################################################

######## Add life history traits for mtdna_new data ########

#add body size from FB

mtdna_spp_new$maxlength <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get length data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$maxlength[i] <- as.numeric(species(mtdna_spp_new$spp[i], fields = 'Length')$Length)
  if(is.na(mtdna_spp_new$maxlength[i])) { #if male maxlength is NA, use female maxlength
    mtdna_spp_new$maxlength[i] <- as.numeric(species(mtdna_spp_new$spp[i], fields = 'LengthFemale')$LengthFemale)
  }
}

#add fecundity from FB
mtdna_spp_new$fecundity <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$fecundity[i] <- fecundity(species=(mtdna_spp_new$spp[i]), field='FecundityMean')
  if(sjmisc::is_empty(unlist(mtdna_spp_new$fecundity[i]), all.na.empty = TRUE)) { #if fecundity mean is NA, use fecundity minimum
    mtdna_spp_new$fecundity[i] <- fecundity(species=mtdna_spp_new$spp[i], field= 'FecundityMin')
    if(sjmisc::is_empty(unlist(mtdna_spp_new$fecundity[i]), all.na.empty = TRUE)) { #if fecundity minimum is NA, use fecundity maximum
      mtdna_spp_new$fecundity[i] <- fecundity(species=(mtdna_spp_new$spp[i]), fields = 'FecundityMax')
    }
  }
  mtdna_spp_new$fecundity_mean[i] <- mean(unlist(mtdna_spp_new$fecundity[i]), na.rm = TRUE) #calculate mean fecundity
}

#add fertilization from FB
mtdna_spp_new$fertilization <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get fertilization data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$fertilization[i] <- reproduction(species=(mtdna_spp_new$spp[i]), field='Fertilization')
}

#add reproduction mode from FB
mtdna_spp_new$reproductionmode <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get reproduction mode data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$reproductionmode[i] <- reproduction(species=(mtdna_spp_new$spp[i]), field='ReproMode')
}

summary(mtdna_spp_new) #see results


################################################### WRITE CSV #################################################################

### Write csv for msat new data ###

msat_spp$X.1 <- NULL #removed unnecessary column
msat_spp$X <- NULL #removed unnecessary column

msat_spp <- apply(msat_spp,2,as.character)

write.csv(msat_spp, "new_msat_full_US_data.csv")

### Write csv for mtdna new data ###

mtdna_spp_new$X.1 <- NULL #removed unnecessary column
mtdna_spp_new$X <- NULL #removed unnecessary column

mtdna_spp_newfinal <- apply(mtdna_spp_new,2,as.character)

write.csv(mtdna_spp_newfinal, "new_mtdna_full_US_data.csv")

