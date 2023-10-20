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
library(FishLife)
library(taxize)
library(stringr)

#read in data
msat_spp <- read.csv("Datasets/msat_papers_USA_new.csv", stringsAsFactors = FALSE) #read in msat species list
mtdna_spp_new <- read.csv("Datasets/mtdna_papers_USA_new.csv", stringsAsFactors = FALSE) #read in mtdna species list

msat_spp <- read.csv("Datasets/new_msat_full_US_data.csv", stringsAsFactors = FALSE) #read in msat species list
mtdna_spp_new <- read.csv("Datasets/new_mtdna_full_US_data.csv", stringsAsFactors = FALSE)

#NOTE: Branchiostoma spps in mtdna NOT in fishbase -- will need to pull out before querying FB for data

##########################################################################################################################################

######## Add life history traits for msat ########

#add body size from FB

msat_spp$maxlength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get length data
  cat(paste(i, " ", sep = ''))
  msat_spp$maxlength[i] <- as.numeric(species(msat_spp$spp[i], fields = 'Length')$Length)
  if(is.na(msat_spp$maxlength[i])) { #if male maxlength is NA, use female maxlength
    msat_spp$maxlength[i] <- as.numeric(species(msat_spp$spp[i], fields = 'LengthFemale')$LengthFemale)
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
  gsub('["]',"", msat_spp$fecundity[i])
  msat_spp$fecundity_mean[i] <- mean(as.numeric(unlist(msat_spp$fecundity[i])), na.rm = TRUE) #calculate mean fecundity
}

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

#add age at first maturity from FB
msat_spp$maturityage <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  msat_spp$maturityage[i] <- maturity(species=(msat_spp$spp[i]), field='tm')
  msat_spp$maturityage_mean[i] <- mean(unlist(msat_spp$maturityage[i]), na.rm = TRUE) #calculate mean maturity length
}

#add generation time/length from FB
msat_spp$generationtime <- NA #create column to fill in
msat_spp$K <- NA 

for (i in 1:nrow(msat_spp)) { #get generation time/length data
  cat(paste(i, " ", sep = ''))
  msat_spp$K[i] <- popgrowth(species=(msat_spp$spp[i]), field='K')
  msat_spp$K_mean[i] <- mean(unlist(msat_spp$K[i]), na.rm = TRUE)
  msat_spp$generationtime[i] <- log(3)/msat_spp$K_mean[i]
}  

#add adult lifespan from FB
msat_spp$adultlifespan <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get adult lifespan data
  cat(paste(i, " ", sep = ''))
  msat_spp$adultlifespan[i] <- species(species=(msat_spp$spp[i]), field='LongevityWild')
}

#pull taxonomy
msat_spp$class <- NA 
msat_spp$class <- tax_name(msat_spp$spp, get = 'class', db = 'itis')$class
msat_spp <- transform(msat_spp,class = replace(class, grepl('Clupea pallasii pallasii', spp, fixed=TRUE), 'Actinopterygii')) #Had to manually add for this spp (info found on NatureServe: https://explorer.natureserve.org/Taxon/ELEMENT_GLOBAL.2.105816/Clupea_pallasii)

msat_spp$order <- NA 
msat_spp$order <- tax_name(msat_spp$spp, get = 'order', db = 'itis')$order
msat_spp <- transform(msat_spp,order = replace(order, grepl('Clupea pallasii pallasii', spp, fixed=TRUE), 'Clupeiformes')) #Had to manually add for this spp

msat_spp$family <- NA 
msat_spp$family <- tax_name(msat_spp$spp, get = 'family', db = 'itis')$family
msat_spp <- transform(msat_spp,family = replace(family, grepl('Clupea pallasii pallasii', spp, fixed=TRUE), 'Clupeidae')) #Had to manually add for this spp

msat_spp$genus <- NA
msat_spp$genus <- tax_name(msat_spp$spp, get = 'genus', db = 'itis')$genus
msat_spp <- transform(msat_spp,genus = replace(genus, grepl('Clupea pallasii pallasii', spp, fixed=TRUE), 'Clupea')) #Had to manually add for this spp

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
  mtdna_spp_new$fecundity[i] <- fecundity(species=(mtdna_spp_new$spp[i]), field='FecundityMean', quotes = FALSE)
  if(sjmisc::is_empty(unlist(mtdna_spp_new$fecundity[i]), all.na.empty = TRUE)) { #if fecundity mean is NA, use fecundity minimum
    mtdna_spp_new$fecundity[i] <- fecundity(species=mtdna_spp_new$spp[i], field= 'FecundityMin')
    if(sjmisc::is_empty(unlist(mtdna_spp_new$fecundity[i]), all.na.empty = TRUE)) { #if fecundity minimum is NA, use fecundity maximum
      mtdna_spp_new$fecundity[i] <- fecundity(species=(mtdna_spp_new$spp[i]), fields = 'FecundityMax')
    }
  } 
  gsub('["]',"", mtdna_spp_new$fecundity[i])
  mtdna_spp_new$fecundity_mean[i] <- mean(as.numeric(unlist(mtdna_spp_new$fecundity[i])), na.rm = TRUE) 
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

#add age at first maturity from FB
mtdna_spp_new$maturityage <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$maturityage[i] <- maturity(species=(mtdna_spp_new$spp[i]), field='tm')
  mtdna_spp_new$maturityage_mean[i] <- mean(unlist(mtdna_spp_new$maturityage[i]), na.rm = TRUE) #calculate mean maturity length
}

#add generation time/length from FB
mtdna_spp_new$generationtime <- NA #create column to fill in

mtdna_spp_new$K <- NA 
for (i in 1:nrow(mtdna_spp_new)) { #get generation time/length data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$K[i] <- popgrowth(species=(mtdna_spp_new$spp[i]), field='K')
  mtdna_spp_new$K_mean[i] <- mean(unlist(mtdna_spp_new$K[i]), na.rm = TRUE)
  mtdna_spp_new$generationtime[i] <- log(3)/mtdna_spp_new$K_mean[i]
}  

#add adult lifespan from FB
mtdna_spp_new$adultlifespan <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get adult lifespan data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$adultlifespan[i] <- species(species=(mtdna_spp_new$spp[i]), field='LongevityWild')
}

#pull taxonomy
mtdna_spp_new$class <- NA 
mtdna_spp_new$class <- tax_name(mtdna_spp_new$spp, get = 'class', db = 'ncbi')$class

mtdna_spp_new$order <- NA 
mtdna_spp_new$order <- tax_name(mtdna_spp_new$spp, get = 'order', db = 'ncbi')$order
mtdna_spp_new <- transform(mtdna_spp_new,order = replace(order, grepl('Cynoscion nebulosus', spp, fixed=TRUE), 'Perciformes')) #Had to manually add for this spp (info found on NatureServe: https://explorer.natureserve.org/Taxon/ELEMENT_GLOBAL.2.101483/Cynoscion_nebulosus)


mtdna_spp_new$family <- NA #create column to fill in
mtdna_spp_new$family <- tax_name(mtdna_spp_new$spp, get = 'family', db = 'ncbi')$family

mtdna_spp_new$genus <- NA
mtdna_spp_new$genus <- tax_name(mtdna_spp_new$spp, get = 'genus', db = 'ncbi')$genus

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

