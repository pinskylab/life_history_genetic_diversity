# Datasets Directory

## 1. Starting CSVs
These csvs are necessary as the start of the entire pathway and were created during intial data collection (see Materials and Methods section in manuscript). Needed to run [`Subsetting to US.R`](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Subsetting%20to%20US%20data.R)

* [msatloci.csv](https://github.com/pinskylab/marial_diversity/blob/master/Datasets/msatloci.csv) : Contains the collected marine species and important information related.
The columns are: 
> 1) Species name
> 2) Common name
> 3) Source
> 4) Primer note- whether the study described new primer pairs useful for studying the species
> 5) Country- here, it will be around continental USA with some sites in parts of northern Mexico 
> 6) Site
> 7) Latitude and longitude
> 8) Stock ID
> 9) Collection year
> 10) Number of markers 
> 11) Marker name 
> 12) Cross spp- indicates whether the marker used for analysis was made for the specific species or for another one instead
> 13) n- samples
> 14) Repeat
> 15) He- heterozygosity
> 16) Hese- heteroozygosity standard error
> 17) File
> 18) He_orig
> 19) SD- standard deviation

* [mtdna_assembled_new.csv](https://github.com/pinskylab/marial_diversity/blob/master/Datasets/mtdna_assembled_new.csv) : Contains the collected marine species and important information related. Columns are similar to msatloci.csv:
> Species name, common name, source, country, site, latitude and longitude, stock ID, collection year, marker name, n, He, Hese, Pi, Pise, and file as well as the inclusion of base pairs (number of base pairs sequenced)

<br />

## 2. Necessary to run [`Pull LH traits.R`](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Pull%20LH%20traits.R)
These CSVs were created from the above CSVs with Subsetting to US.R script.

* [msat_papers_USA_new.csv](https://github.com/pinskylab/marial_diversity/blob/master/Datasets/msat_papers_USA_new.csv) : Contains the collected marine species and important information related subsetted to the appropriate coordinates. Columns are same as the starting CSVs: 
> Species name, common name, source, primer note, country, site, latitude and longitude, stock ID, collection year, number of markers, marker name, cross spp, n, repeat, He, Hese, and file

* [mtdna_papers_USA_new.csv](https://github.com/pinskylab/marial_diversity/blob/master/Datasets/mtdna_papers_USA_new.csv) : Contains the collected marine species and important information related subsetted to the appropriate coordinates. Columns are same as the starting CSVs: 
> Species name, common name, source, country, site, latitude and longitude, stock ID, collection year, marker name, n, base pairs, He, Hese, Pi, Pise, and file 

<br />

## 3. Used for the remaining scripts

These csvs were created the previous step's CSVs + Pull LH traits.R script and are needed to run [`US Graphs`](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Subsetting%20to%20US%20data.R), [`Mixed Models for US Data (msat & mtDNA)`](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Mixed%20Models%20for%20US%20Data%20(msat%20%26%20mtDNA).R), & [`Data Characteristics`](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Data%20charactersitics.R).

* [new_msat_full_US_data.csv](https://github.com/pinskylab/marial_diversity/blob/master/Datasets/new_msat_full_US_data.csv) : Contains the necessary information for the species to run analyses. Columns are similar to msat_papers_USA_new.csv with the inclusion of the traits used in the study:
> Species name, common name, source, primer note, country, site, latitude and longitude, stock ID, collection year, number of markers, marker name, cross spp, n, repeat, He, Hese, file PLUS the traits shown below &#8595;
> 
> 1) **Maxlength**- size of the longest individual recorded
> 2) **Fecundity (absolute)**- total number of eggs in a female
> 3) **Fecundity mean**- mean of fecundity data
> 4) **Fertilization**- whether species' external or internal
> 5) **Reproduction mode**- whether species' reproduced via dioecism, protandry, protogyny, true hermaphroditism, or parthenogenesis (for this study, only dioecism and protogyny is relevant)

* [new_mtdna_full_US_data.csv](https://github.com/pinskylab/marial_diversity/blob/master/Datasets/new_mtdna_full_US_data.csv) : Contains the necessary information for the species to run analyses. Columns are similar to mtdna_papers_USA_new.csv with the inclusion of the traits used in the study same as new_msat_full_US_data.csv:
> Species name, common name, source, country, site, latitude and longitude, stock ID, collection year, marker name, n, base pairs, He, Hese, Pi, Pise, file, PLUS **maxlength, fecundity, fecundity mean, fertilizatin, and reproductive mode**

<br />

## Misc. CSVs

* [IUCN_status.csv](https://github.com/pinskylab/marial_diversity/blob/master/Datasets/IUCN_status.csv)
Contains the IUCN status of all the species used in this study.
