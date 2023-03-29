# Scripts Directory
Contains description of each script and is placed in the order needed to run the next step. *(Note: steps 3, 4, & 5 are interchangeable as they all start off with the same csv created from step 2).*


### 1. [Subsetting to US data.R](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Subsetting%20to%20US%20data.R)
Used to subset marine species' data from the original microsatellite and mitchondrial DNA to just Continental USA. The range taken for latitude was from 23° to 50° while longitude range was from -128° to -65°. This script starts of the rest of the pathway.


### 2. [Pull LH traits.R](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Pull%20LH%20traits.R)
Used to pull life history traits for both microsatellite and mitchondrial DNA species from Fishbase.


### 3. [US Graphs.R](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Subsetting%20to%20US%20data.R)
Used to analyze and graph data for microsatellite and mitchondrial DNA US data. Includes steps to graph:
* Fertilization method
* Reproduction mode
* Maps of Continental USA with sample locations plotted  

### 4. [Mixed Models for US Data (msat & mtDNA).R](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Mixed%20Models%20for%20US%20Data%20(msat%20%26%20mtDNA).R)
Used to analyze and create mixed models for microsatellite and mitchondrial DNA samples. It also includes steps to graph:
* Maximum length
* Fecundity

### 5. [Data Characteristics](https://github.com/pinskylab/marial_diversity/blob/master/Scripts/Data%20charactersitics.R)
Used to look at the datasets and see the different characteristics. For example: looking at number of sources or number of each reproduction category.
