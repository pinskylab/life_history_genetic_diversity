################################################### Script for Subsetting to US Data ########################################################

#subsetting marine species' data to Continental USA 

##########################################################################################################################################

######### Set-up #########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(dplyr)

#read in data
mtdna_papers_new <- read.csv("Datasets/mtdna_assembled_new.csv", stringsAsFactors = FALSE) #read in 
msat_papers_new <- read.csv("Datasets/msatloci.csv", stringsAsFactors = FALSE) #read in 

##########################################################################################################################################

############### Isolating sites to only Continental USA ############### 

#### Subset data by latitude and longitude columns in both data sets ####

#Subset mtdna_papers_new data
colnames(mtdna_papers_new) #get column names to input for the following line of code

mtdna_papers_latlon_new <- mtdna_papers_new %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve = TRUE) #subset data between the given latitudes and longitudes

#Subset msat data
colnames(msat_papers_new) #get column names to input for the following line of code

msat_papers_latlon_new <- msat_papers_new %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve = TRUE) #subset data between the given latitudes and longitudes

### Subset data further by Country name in both sets ###

#Subset mtdna_papers_new data
unique(mtdna_papers_latlon_new$Country ) #find unique values to allow for the subsetting to just continental USA

mtdna_papers_USA_new <- mtdna_papers_latlon_new  %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(Country %in% c("United States", "USA", "USA, Western North Atlantic","Northwestern Pacific Ocean "), .preserve = TRUE) #subset data for just continental USA

#Subset msat data
unique(msat_papers_latlon_new$Country ) #find unique values to allow for the subsetting to just continental USA

msat_papers_USA_new <- msat_papers_latlon_new %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(Country %in% c("USA", "United States"), .preserve = TRUE) #subset data for just continental USA

### Double check to see if it was subsetted correctly by doing the reverse process ###

### Subset data further by Country name in both sets ###

#Subset mtdna_papers_new data
unique(mtdna_papers_new$Country ) #find unique values to allow for the subsetting to just continental USA

mtdna_papers_latlon2_new <- mtdna_papers_new  %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(Country %in% c("United States", "USA", "USA, Western North Atlantic","Northwestern Pacific Ocean "), .preserve = TRUE) #subset data for just continental USA

#Subset msat data
unique(msat_papers_new$Country ) #find unique values to allow for the subsetting to just continental USA

msat_papers_latlon2_new <- msat_papers_new %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(Country %in% c("USA", "United States"), .preserve = TRUE) #subset data for just continental USA

### Subset data by latitude and longitude columns in both data sets ###

#Subset mtdna data
colnames(mtdna_papers_latlon2_new) #get column names to input for the following line of code

mtdna_papers_USA2_new <- mtdna_papers_latlon2_new %>% 
  select(X, spp, CommonName, Source, Country, Site, lat, lon, stockid, CollectionYear, MarkerName, n, bp, He, Hese, Pi, Pise, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve=TRUE) #subset data between the given latitudes and longitudes

#Subset msat data
colnames(msat_papers_latlon2_new) #get column names to input for the following line of code

msat_papers_USA2_new <- msat_papers_latlon2_new %>% 
  select(X, spp, CommonName, Source, PrimerNote, Country, Site, lat, lon, stockid, CollectionYear, NumMarkers, MarkerName, CrossSpp, n, Repeat, He, Hese, file) %>%
  filter(between(lat, 23, 50) & between(lon, -128, -65), .preserve=TRUE) #subset data between the given latitudes and longitudes

#### Write subsetted data as a CSV file ####

write.csv(mtdna_papers_USA_new,'mtdna_papers_USA_new.csv') #write csv for mtDNA new data

write.csv(msat_papers_USA_new,'msat_papers_USA_new.csv') #write csv for msat data