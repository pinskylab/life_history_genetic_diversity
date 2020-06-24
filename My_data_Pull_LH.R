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
  msat_spp$mortalitywild[i] <- species((msat_spp$fbsci[i]), field='LongevityWild')
}


#add fecundity from FB
msat_spp$fecundity <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  msat_spp$fecundity[i] <- fecundity(species=(msat_spp$fbsci[i]), field=('FecundityMean'))
}

#add fecundity min from FB
msat_spp$fecunditymin <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  msat_spp$fecunditymin[i] <- fecundity(species=(msat_spp$fbsci[i]), field=('FecundityMin'))
}

#add fecundity max from FB
msat_spp$fecunditymax <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  msat_spp$fecunditymax[i] <- fecundity(species=(msat_spp$fbsci[i]), field=('FecundityMax'))
}

#add spawning cycles from FB
msat_spp$spawningcycles <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  msat_spp$spawningcycles[i] <- spawning(species=(msat_spp$fbsci[i]), field='SpawningCycles')
}

#add spawning cycles from FB
msat_spp$numberofspawningcycles <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  msat_spp$numberofspawningcycles[i] <- reproduction(species=(msat_spp$fbsci[i]), field='Spawning')
}

#add reproductive season length from FB
msat_spp$reproductiveseasonlength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  msat_spp$reproductiveseasonlength[i] <- spawning(species=(msat_spp$fbsci[i]), field='Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
}

#add maturity length from FB
msat_spp$maturitylength <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get maturity length data
  cat(paste(i, " ", sep = ''))
  msat_spp$maturitylength[i] <- maturity(species=(msat_spp$fbsci[i]), field='LengthMatMin')
}

#add maturity age from FB
msat_spp$maturityage <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  msat_spp$maturityage[i] <- maturity(species=(msat_spp$fbsci[i]), field='AgeMatMin')
}

#add fertilization from FB
msat_spp$mortalitywild <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get fertilization data
  cat(paste(i, " ", sep = ''))
  msat_spp$mortalitywild[i] <- as.numeric(species=(msat_spp$fbsci[i]), field='LongevityWild')
}

#add parental care from FB
msat_spp$parentalcare <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get parental care data
  cat(paste(i, " ", sep = ''))
  msat_spp$parentalcare[i] <- reproduction(species=(msat_spp$fbsci[i]), field='ParentalCareQ')
}

#add reproduction mode from FB
msat_spp$reproductionmode <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get reproduction mode data
  cat(paste(i, " ", sep = ''))
  msat_spp$reproductionmode[i] <- reproduction(species=(msat_spp$fbsci[i]), field='ReproMode')
}

#add spawning ground from FB
msat_spp$spawningground <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get spawning ground data
  cat(paste(i, " ", sep = ''))
  msat_spp$spawningground[i] <- spawning(species_list=(msat_spp$fbsci[i]), field='SpawningGround')
}

#add Batch Spawner information from FB
msat_spp$batchspawner <- NA #create column to fill in

for (i in 1:nrow(msat_spp)) { #get batch spawner data
  cat(paste(i, " ", sep = ''))
  msat_spp$batchspawner[i] <- reproduction(species_list=(msat_spp$fbsci[i]), field='BatchSpawner')
}


summary(msat_spp)

#Ex of code without if statement
#for(i in 1:nrow(spps2)) { #get code that specifies exact species on FB
#  cat(paste(i, " ", sep = ''))
#  spps2$SpecCode[i] <- as.numeric(species(specieslist = spps2$fbsci[i], fields = 'SpecCode')$SpecCode)
#}

##########################################################################################################################################

######## Add life history traits for mtdna ########
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
  mtdna_spp$fecundity[i] <- fecundity(species=(mtdna_spp$fbsci[i]), field=('FecundityMean'))
}

#add fecundity min from FB
mtdna_spp$fecunditymin <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$fecunditymin[i] <- fecundity(species=(mtdna_spp$fbsci[i]), field=('FecundityMin'))
}

#add fecundity max from FB
mtdna_spp$fecunditymax <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get fecundity data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$fecunditymax[i] <- fecundity(species=(mtdna_spp$fbsci[i]), field=('FecundityMax'))
}

#add spawning cycles from FB
mtdna_spp$spawningcycles <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$spawningcycles[i] <- spawning(species=(mtdna_spp$fbsci[i]), field='SpawningCycles')
}

#add spawning cycles from FB
mtdna_spp$numberofspawningcycles <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$numberofspawningcycles[i] <- reproduction(species=(mtdna_spp$fbsci[i]), field='Spawning')
}

#add spawning ground from FB
mtdna_spp$spawningground <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get spawning ground data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$spawningground[i] <- spawning(species_list=(mtdna_spp$fbsci[i]), field='SpawningGround')
}

#add reproductive season length from FB
mtdna_spp$reproductiveseasonlength <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get spawning cycles data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$reproductiveseasonlength[i] <- spawning(species=(mtdna_spp$fbsci[i]), field='Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
}

#add maturity length from FB
mtdna_spp$maturitylength <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get maturity length data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$maturitylength[i] <- maturity(species=(mtdna_spp$fbsci[i]), field='LengthMatMin')
}

#add maturity age from FB
mtdna_spp$maturityage <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get maturity age data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$maturityage[i] <- maturity(species=(mtdna_spp$fbsci[i]), field='AgeMatMin')
}

#add fertilization from FB
mtdna_spp$mortalitywild <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get fertilization data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$mortalitywild[i] <- as.numeric(species=(mtdna_spp$fbsci[i]), field='LongevityWild')
}

#add parental care from FB
mtdna_spp$parentalcare <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get parental care data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$parentalcare[i] <- reproduction(species=(mtdna_spp$fbsci[i]), field='ParentalCareQ')
}

#add reproduction mode from FB
mtdna_spp$reproductionmode <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get reproduction mode data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$reproductionmode[i] <- reproduction(species=(mtdna_spp$fbsci[i]), field='ReproMode')
}

####add if it is a Hermaphrodite: FIX THIS
mtdna_spp$hermaphrodite <- NA

for(i in 1:nrow(mtdna_spp)){
  cat(paste(i, " ", sep = ''))
within(mtdna_spp, mtdna_spp$hermaphrodite[i] <- ifelse(reproductionmode=="protandry","Yes", 
                          ifelse(reproductionmode=="protogyny", "Yes",
                                 ifelse(reproductionmode=="true hermaphroditism","Yes",
                                        ifelse(reproductionmode=="dioecism","No",
                                               ifelse(reproductionmode== is.na(NA), NA))))) )
}
  
#option 2: for(i in 1:nrow(mtdna_spp))
#cat(paste(i, " ", sep = ''))
#mtdna_spp$hermaphrodite[i]=ifelse (mtdna_spp$reproductionmode&c("protandry" | "protogyny" | "true hermaphroditism"),TRUE,FALSE) {
#mtdna_spp$hermaphrodite[i] <- ('Yes')
#} else 
#hermaphrodite<-print('No')


#add Batch Spawner information from FB
mtdna_spp$batchspawner <- NA #create column to fill in

for (i in 1:nrow(mtdna_spp)) { #get batch spawner data
  cat(paste(i, " ", sep = ''))
  mtdna_spp$batchspawner[i] <- reproduction(species_list=(mtdna_spp$fbsci[i]), field='BatchSpawner')
}


summary(mtdna_spp)

#Ex of code without if statement
#for(i in 1:nrow(spps2)) { #get code that specifies exact species on FB
#  cat(paste(i, " ", sep = ''))
#  spps2$SpecCode[i] <- as.numeric(species(specieslist = spps2$fbsci[i], fields = 'SpecCode')$SpecCode)
#}