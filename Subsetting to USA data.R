################################################### Script for Subsetting to USA Data ########################################################

#subsetting marine species' data to Continental USA 

##########################################################################################################################################

######### Set-up #########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(dplyr)

#read in data
mtdna_papers <- read.csv("mtdna_assembled_2.csv", stringsAsFactors = FALSE) #read in 
msat_papers <- read.csv("msat_assembled_2.csv", stringsAsFactors = FALSE) #read in 

##########################################################################################################################################

############### Isolating sites to only Continental USA ############### 

#### Subset data by latitude and longitude columns in both data sets ####

#Subset mtdna data
colnames(mtdna_papers) #get column names to input for the following line of code

mtdna_papers_latlon <- mtdna_papers %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve = TRUE) #subset data between the given latitudes and longitudes

#Subset msat data
colnames(msat_papers) #get column names to input for the following line of code

msat_papers_latlon <- msat_papers %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve = TRUE) #subset data between the given latitudes and longitudes

### Subset data further by Country name in both sets ###

#Subset mtdna data
unique(mtdna_papers_latlon$Country ) #find unique values to allow for the subsetting to just continental USA

mtdna_papers_USA <- mtdna_papers_latlon  %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(Country %in% c("United States", "USA", "USA, Western North Atlantic","Northwestern Pacific Ocean "), .preserve = TRUE) #subset data for just continental USA

#Subset msat data
unique(msat_papers_latlon$Country ) #find unique values to allow for the subsetting to just continental USA

msat_papers_USA <- msat_papers_latlon %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(Country %in% c("USA", "United States"), .preserve = TRUE) #subset data for just continental USA

### Double check to see if it was subsetted correctly by doing the reverse process ###

### Subset data further by Country name in both sets ###

#Subset mtdna data
unique(mtdna_papers$Country ) #find unique values to allow for the subsetting to just continental USA

mtdna_papers_latlon2 <- mtdna_papers  %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(Country %in% c("United States", "USA", "USA, Western North Atlantic","Northwestern Pacific Ocean ")), .preserve = TRUE) #subset data for just continental USA

#Subset msat data
unique(msat_papers$Country ) #find unique values to allow for the subsetting to just continental USA

msat_papers_latlon2 <- msat_papers %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(Country %in% c("USA", "United States"), .preserve = TRUE) #subset data for just continental USA

### Subset data by latitude and longitude columns in both data sets ###

#Subset mtdna data
colnames(mtdna_papers_latlon2) #get column names to input for the following line of code

mtdna_papers_USA2 <- mtdna_papers_latlon2 %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve=TRUE) #subset data between the given latitudes and longitudes

#Subset msat data
colnames(msat_papers_latlon2) #get column names to input for the following line of code

msat_papers_USA2 <- msat_papers_latlon2 %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve=TRUE) #subset data between the given latitudes and longitudes

#### Write subsetted data as a CSV file ####

write.csv(mtdna_papers_USA,'mtdna_papers_USA.csv') #write csv for mtDNA data

write.csv(msat_papers_USA,'msat_papers_USA.csv') #write csv for msat data

##########################################################################################################################################

############### Merge data from mtdna_spp_final.csv and mtdna_papers_USA.csv togther ############### 

mtdna_spp_final_1 = read.csv("mtdna_spp_final.csv", header=TRUE, sep=",")
mtdna_papers_USA_2 = read.csv("mtdna_papers_USA.csv", header=TRUE, sep=",")

mtdna_full_data_USA = merge( mtdna_papers_USA_2, mtdna_spp_final_1, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE, by="spp")

mtdna_full_data_USA$X.1 <- NULL #removed unnecessary column
mtdna_full_data_USA$X.x <- NULL #removed unnecessary column

mtdna_full_data_USA <- mtdna_full_data_USA[!is.na(mtdna_full_data_USA$CommonName),] #omit any species data from mtdna_spp_final that was not found in mtdna_paper_USA_2 to line up better

###write csv for mtdna data###

write.csv(mtdna_full_data_USA, "mtdna_full_US_data.csv")

######## Merge data from msat_spp_final.csv and msat_assembled.csv togther########

msat_spp_final_1 = read.csv("msat_spp_final.csv", header=TRUE, sep=",")
msat_papers_USA._2 = read.csv("msat_papers_USA.csv", header=TRUE, sep=",")

msat_full_data_USA = merge(msat_papers_USA._2, msat_spp_final_1, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE, by="spp")

msat_full_data_USA$X.1 <- NULL #removed unnecessary column
msat_full_data_USA$X.x <- NULL #removed unnecessary column

msat_full_data_USA <- msat_full_data_USA[!is.na(msat_full_data_USA$CommonName),] #omit any species data from mtdna_spp_final that was not found in mtdna_paper_USA_2 to line up better

###write csv###

write.csv(msat_full_data_USA, "msat_full_US_data.csv")

