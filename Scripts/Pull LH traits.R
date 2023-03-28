################################################### Script for Pulling Life History Traits  #######################################################

#pulls life history traits for both mtdna & msat species

##########################################################################################################################################

####### Set-up #######

remove(list = ls())

#load libraries
library(rfishbase)
library(tidyverse)
library(dplyr)

#read in data
msat_spp <- read.csv("msat_papers_USA_new.csv", stringsAsFactors = FALSE) #read in msat species list
mtdna_spp_new <- read.csv("mtdna_papers_USA_new.csv", stringsAsFactors = FALSE) #read in mtdna species list for new 

#NOTE: Branchiostoma spps in mtdna NOT in fishbase -- will need to pull out before querying FB for data

##########################################################################################################################################

######## Add life history traits for msat ########

######## Add body size ########

#add body size from FB

msat_spp$maxlength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get length data
  cat(paste(i, " ", sep = ''))
  msat_spp$maxlength[i] <- as.numeric(species(msat_spp$spp[i], fields = 'Length')$Length)
  if(is.na(msat_spp$maxlength[i])) { #if male maxlength is NA, use female maxlength
    msat_spp$maxlength[i] <- as.numeric(species(msat_spp$spp[i], fields = 'LengthFemale')$LengthFemale)
  }
}

#add mortality from FB
msat_spp$mortalitywild <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get mortality data
  cat(paste(i, " ", sep = ''))
  msat_spp$mortalitywild[i] <- species((msat_spp$spp[i]), field='LongevityWild')
}

#add fecundity from FB
library(sjmisc) #install for correct is_empty() function
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

#add maturity length from FB
msat_spp$maturitylength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get maturity length data
  cat(paste(i, " ", sep = ''))
  msat_spp$maturitylength[i] <- maturity(species=(msat_spp$spp[i]), field='LengthMatMin')
  msat_spp$maturitylength_mean[i] <- mean(unlist(msat_spp$maturitylength[i]), na.rm = TRUE) #calculate mean maturity length
}

#add maturity age from FB
msat_spp$maturityage <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  msat_spp$maturityage[i] <- maturity(species=(msat_spp$spp[i]), field='AgeMatMin')
  msat_spp$maturityage_mean[i] <- mean(unlist(msat_spp$maturityage[i]), na.rm = TRUE) #calculate mean maturity length
}

#add spawning cycles from FB
msat_spp$numberofspawningcycles <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  msat_spp$numberofspawningcycles[i] <- reproduction(species=(msat_spp$spp[i]), field='Spawning')
}


#add reproductive season length from FB
msat_spp$reproductiveseasonlength <- NA #create column to fill in

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

#add spawning ground from FB
msat_spp$spawningground <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get spawning ground data
  cat(paste(i, " ", sep = ''))
  msat_spp$spawningground[i] <- spawning(species_list=(msat_spp$spp[i]), field='SpawningGround')
}

#add Batch Spawner information from FB
msat_spp$batchspawner <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get batch spawner data
  cat(paste(i, " ", sep = ''))
  msat_spp$batchspawner[i] <- reproduction(species_list=(msat_spp$spp[i]), field='BatchSpawner')
}

#add parental care from FB
msat_spp$parentalcare <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get parental care data
  cat(paste(i, " ", sep = ''))
  msat_spp$parentalcare[i] <- reproduction(species=(msat_spp$spp[i]), field='ParentalCareQ')
}

#add larvae size from FB
msat_spp$larvaesize <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get larvae data
  cat(paste(i, " ", sep = ''))
  msat_spp$larvaesize[i] <- larvae(species=(msat_spp$spp[i]), field='LhMid')
  if(sjmisc::is_empty(unlist(msat_spp$larvaesize[i]), all.na.empty = TRUE)) { #if larvae median is NA, use larvae size minimum
    msat_spp$larvaesize[i] <- larvae(species=msat_spp$spp[i], field= 'LhMin')
    if(sjmisc::is_empty(unlist(msat_spp$larvaesize[i]), all.na.empty = TRUE)) { #if larvae minimum is NA, use larvae size maximum
      msat_spp$larvaesize[i] <- larvae(species=(msat_spp$spp[i]), fields = 'LhMax')
    }
  }
  msat_spp$larvae_mean[i] <- mean(unlist(msat_spp$larvaesize[i]), na.rm = TRUE) #calculate mean larvae size
}

summary(msat_spp)

#Ex of code without if statement
#for(i in 1:nrow(spps2)) { #get code that specifies exact species on FB
#  cat(paste(i, " ", sep = ''))
#  spps2$SpecCode[i] <- as.numeric(species(specieslist = spps2$fbsci[i], fields = 'SpecCode')$SpecCode)
#}

##########################################################################################################################################

######## Add life history traits for mtdna_spp ########
#add body size from FB

mtdna_spp$maxlength <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get length data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$maxlength[i] <- as.numeric(species(mtdna_spp$fbsci[i], fields = 'Length')$Length)
  if(is.na(mtdna_spp$maxlength[i])) { #if male maxlength is NA, use female maxlength
    mtdna_spp$maxlength[i] <- as.numeric(species(mtdna_spp$fbsci[i], fields = 'LengthFemale')$LengthFemale)
  }
}

#add mortality from FB
mtdna_spp$mortalitywild <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get mortality data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$mortalitywild[i] <- species((mtdna_spp$fbsci[i]), field='LongevityWild')
}


#add fecundity from FB
mtdna_spp$fecundity <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$fecundity[i] <- fecundity(species=(mtdna_spp$fbsci[i]), field='FecundityMean')
  if(sjmisc::is_empty(unlist(mtdna_spp$fecundity[i]), all.na.empty = TRUE)) { #if fecundity mean is NA, use fecundity minimum
    mtdna_spp$fecundity[i] <- fecundity(species=mtdna_spp$fbsci[i], field= 'FecundityMin')
    if(sjmisc::is_empty(unlist(mtdna_spp$fecundity[i]), all.na.empty = TRUE)) { #if fecundity minimum is NA, use fecundity maximum
      mtdna_spp$fecundity[i] <- fecundity(species=(mtdna_spp$fbsci[i]), fields = 'FecundityMax')
    }
  }
  mtdna_spp$fecundity_mean[i] <- mean(unlist(mtdna_spp$fecundity[i]), na.rm = TRUE) #calculate mean fecundity
}

#add maturity length from FB
mtdna_spp$maturitylength <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get maturity length data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$maturitylength[i] <- maturity(species=(mtdna_spp$fbsci[i]), field='LengthMatMin')
  mtdna_spp$maturitylength_mean[i] <- mean(unlist(mtdna_spp$maturitylength[i]), na.rm = TRUE) #calculate mean maturity length
}

#add maturity age from FB
mtdna_spp$maturityage <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$maturityage[i] <- maturity(species=(mtdna_spp$fbsci[i]), field='AgeMatMin')
  mtdna_spp$maturityage_mean[i] <- mean(unlist(mtdna_spp$maturityage[i]), na.rm = TRUE) #calculate mean maturity length
}

#add spawning cycles from FB
mtdna_spp$numberofspawningcycles <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$numberofspawningcycles[i] <- reproduction(species=(mtdna_spp$fbsci[i]), field='Spawning')
}

#add reproductive season length from FB
mtdna_spp$reproductiveseasonlength <- NA #create column to fill in

#add fertilization from FB
mtdna_spp$fertilization <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get fertilization data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$fertilization[i] <- reproduction(species=(mtdna_spp$fbsci[i]), field='Fertilization')
}

#add reproduction mode from FB
mtdna_spp$reproductionmode <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get reproduction mode data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$reproductionmode[i] <- reproduction(species=(mtdna_spp$fbsci[i]), field='ReproMode')
}

#add spawning ground from FB
mtdna_spp$spawningground <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get spawning ground data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$spawningground[i] <- spawning(species_list=(mtdna_spp$fbsci[i]), field='SpawningGround')
}

#add Batch Spawner information from FB
mtdna_spp$batchspawner <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get batch spawner data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$batchspawner[i] <- reproduction(species_list=(mtdna_spp$fbsci[i]), field='BatchSpawner')
}

#add parental care from FB
mtdna_spp$parentalcare <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get parental care data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$parentalcare[i] <- reproduction(species=(mtdna_spp$fbsci[i]), field='ParentalCareQ')
}

#add larvae size from FB
mtdna_spp$larvaesize <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get larvae data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$larvaesize[i] <- larvae(species=(mtdna_spp$fbsci[i]), field='LhMid')
  if(sjmisc::is_empty(unlist(mtdna_spp$larvaesize[i]), all.na.empty = TRUE)) { #if larvae median is NA, use larvae size minimum
    mtdna_spp$larvaesize[i] <- larvae(species=mtdna_spp$fbsci[i], field= 'LhMin')
    if(sjmisc::is_empty(unlist(mtdna_spp$larvaesize[i]), all.na.empty = TRUE)) { #if larvae minimum is NA, use larvae size maximum
      mtdna_spp$larvaesize[i] <- larvae(species=(mtdna_spp$fbsci[i]), fields = 'LhMax')
    }
  }
  mtdna_spp$larvae_mean[i] <- mean(unlist(mtdna_spp$larvaesize[i]), na.rm = TRUE) #calculate mean larvae size
}

summary(mtdna_spp)

#Ex of code without if statement
#for(i in 1:nrow(spps2)) { #get code that specifies exact species on FB
#  cat(paste(i, " ", sep = ''))
#  spps2$SpecCode[i] <- as.numeric(species(specieslist = spps2$fbsci[i], fields = 'SpecCode')$SpecCode)
#}

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

#add mortality from FB
mtdna_spp_new$mortalitywild <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get mortality data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$mortalitywild[i] <- species((mtdna_spp_new$spp[i]), field='LongevityWild')
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

#add maturity length from FB
mtdna_spp_new$maturitylength <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get maturity length data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$maturitylength[i] <- maturity(species=(mtdna_spp_new$spp[i]), field='LengthMatMin')
  mtdna_spp_new$maturitylength_mean[i] <- mean(unlist(mtdna_spp_new$maturitylength[i]), na.rm = TRUE) #calculate mean maturity length
}

#add maturity age from FB
mtdna_spp_new$maturityage <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$maturityage[i] <- maturity(species=(mtdna_spp_new$spp[i]), field='AgeMatMin')
  mtdna_spp_new$maturityage_mean[i] <- mean(unlist(mtdna_spp_new$maturityage[i]), na.rm = TRUE) #calculate mean maturity length
}

#add spawning cycles from FB
mtdna_spp_new$numberofspawningcycles <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$numberofspawningcycles[i] <- reproduction(species=(mtdna_spp_new$spp[i]), field='Spawning')
}

#add reproductive season length from FB
mtdna_spp_new$reproductiveseasonlength <- NA #create column to fill in

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

#add spawning ground from FB
mtdna_spp_new$spawningground <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get spawning ground data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$spawningground[i] <- spawning(species_list=(mtdna_spp_new$spp[i]), field='SpawningGround')
}

#add Batch Spawner information from FB
mtdna_spp_new$batchspawner <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get batch spawner data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$batchspawner[i] <- reproduction(species_list=(mtdna_spp_new$spp[i]), field='BatchSpawner')
}

#add parental care from FB
mtdna_spp_new$parentalcare <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get parental care data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$parentalcare[i] <- reproduction(species=(mtdna_spp_new$spp[i]), field='ParentalCareQ')
}

#add larvae size from FB
mtdna_spp_new$larvaesize <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp_new)) { #get larvae data
  cat(paste(i, " ", sep = ''))
  mtdna_spp_new$larvaesize[i] <- larvae(species=(mtdna_spp_new$spp[i]), field='LhMid')
  if(sjmisc::is_empty(unlist(mtdna_spp_new$larvaesize[i]), all.na.empty = TRUE)) { #if larvae median is NA, use larvae size minimum
    mtdna_spp_new$larvaesize[i] <- larvae(species=mtdna_spp_new$spp[i], field= 'LhMin')
    if(sjmisc::is_empty(unlist(mtdna_spp_new$larvaesize[i]), all.na.empty = TRUE)) { #if larvae minimum is NA, use larvae size maximum
      mtdna_spp_new$larvaesize[i] <- larvae(species=(mtdna_spp_new$spp[i]), fields = 'LhMax')
    }
  }
  mtdna_spp_new$larvae_mean[i] <- mean(unlist(mtdna_spp_new$larvaesize[i]), na.rm = TRUE) #calculate mean larvae size
}

summary(mtdna_spp_new)

#Ex of code without if statement
#for(i in 1:nrow(spps2)) { #get code that specifies exact species on FB
#  cat(paste(i, " ", sep = ''))
#  spps2$SpecCode[i] <- as.numeric(species(specieslist = spps2$fbsci[i], fields = 'SpecCode')$SpecCode)
#}

##########################################################################################################################################

######## Merge data from mtdna_spp_final.csv and mtdna_assembled.csv togther########

mtdna_spp_final_1 = read.csv("mtdna_spp_final.csv", header=TRUE, sep=",")
mtdna_assembled_2 = read.csv("mtdna_assembled.csv", header=TRUE, sep=",")

mtdna_full_data = merge(mtdna_spp_final_1, mtdna_assembled_2, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE, by="spp")

mtdna_full_data$X.x <- NULL #removed unnecessary column
mtdna_full_data$X.y <- NULL #removed unnecessary column

###write csv for mtdna data###

write.csv(mtdna_full_data, "mtdna_final_data.csv")

###write csv for mtdna new data (no merging required###

mtdna_spp_new$X.1 <- NULL #removed unnecessary column
mtdna_spp_new$X <- NULL #removed unnecessary column

mtdna_spp_newfinal <- apply(mtdna_spp_new,2,as.character)

write.csv(mtdna_spp_newfinal, "new_mtdna_full_US_data.csv")

######## Merge data from msat_spp_final.csv and msat_assembled.csv togther########

msat_spp_final_1 = read.csv("msat_spp_final.csv", header=TRUE, sep=",")
msat_assembled_2 = read.csv("msat_assembled.csv", header=TRUE, sep=",")

msat_full_data = merge(msat_spp_final_1, msat_assembled_2, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE, by="spp")

msat_full_data$X.x <- NULL #removed unnecessary column
msat_full_data$X.y <- NULL #removed unnecessary column

###write csv###

write.csv(msat_full_data, "new_msat_final_data.csv")

###write csv for mtdna new data (no merging required###

msat_spp$X.1 <- NULL #removed unnecessary column
msat_spp$X <- NULL #removed unnecessary column

msat_spp <- apply(msat_spp,2,as.character)

write.csv(msat_spp, "new_msat_full_US_data.csv")
