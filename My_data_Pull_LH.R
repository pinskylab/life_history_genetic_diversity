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

#add mortality from FB
msat_spp$mortalitywild <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get mortality data
  cat(paste(i, " ", sep = ''))
  msat_spp$mortalitywild[i] <- as.numeric(species(msat_spp$fbsci[i]), field='LongevityWild')
}

#add fecundity from FB
msat_spp$fecundity <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  msat_spp$fecundity[i] <- fecundity(species=(msat_spp$fbsci[i]), field=('FecundityMean'))
}

#add spawning cycles from FB
msat_spp$spawningcycles <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  msat_spp$spawningcycles[i] <- as.numerica(fecundity(species(msat_spp$fbsci[i]), field='SpawningCycles'))
}

#add maturity length from FB
msat_spp$maturitylength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get maturity length data
  cat(paste(i, " ", sep = ''))
  msat_spp$maturitylength[i] <- as.numerica(maturity(species(msat_spp$fbsci[i]), field='LengthMatMin'))
}

#add maturity age from FB
msat_spp$maturityage <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  msat_spp$maturityage[i] <- as.numerica(maturity(species(msat_spp$fbsci[i]), field='AgeMatMin'))
}

#add fertilization from FB
msat_spp$mortalitywild <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get fertilization data
  cat(paste(i, " ", sep = ''))
  msat_spp$mortalitywild[i] <- as.numerica(species(msat_spp$fbsci[i]), field='LongevityWild')
}

#add parental care from FB
msat_spp$parentalcare <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get parental care data
  cat(paste(i, " ", sep = ''))
  msat_spp$parentalcare[i] <- reproduction(species(msat_spp$fbsci[i]), field='ParentalCareQ')
}

#add reproduction mode from FB
msat_spp$reproductionmode <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get reporudction mode data
  cat(paste(i, " ", sep = ''))
  msat_spp$reproductionmode[i] <- reproduction(species(msat_spp$fbsci[i]), field='ReproMode')
}

#add spawning ground from FB
msat_spp$spawningground <- NA #create column to fill in
for (i in 1:nrow(msat_spp)) { #get spawning ground data
  cat(paste(i, " ", sep = ''))
  msat_spp$spawningground[i] <- spawning(species_list=(msat_spp$fbsci[i]), field='SpawningGround')
}

summary(msat_spp)

#Ex of code without if statement
#for(i in 1:nrow(spps2)) { #get code that specifies exact species on FB
#  cat(paste(i, " ", sep = ''))
#  spps2$SpecCode[i] <- as.numeric(species(specieslist = spps2$fbsci[i], fields = 'SpecCode')$SpecCode)
#}