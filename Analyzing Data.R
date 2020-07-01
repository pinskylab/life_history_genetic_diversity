################################################### Script for Analyzing Data  #######################################################

#analyzes data 

##########################################################################################################################################

######## Set-up ########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(ggplot2)

#read in data
mtdna_data <- read.csv("mtdna_final_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("msat_final_data.csv", stringsAsFactors = FALSE) #read in 

############### mtDNA data set ############### 

#####Box Plots: Character Data#####

##Fertilization##

#convert "in brood pouch or similar structure" to internal fertilization

mtdna_data$final_fertilization <- NA

mtdna_data$final_fertilization [mtdna_data$fertilization =="external"]  <- "external"
mtdna_data$final_fertilization [mtdna_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
mtdna_data$final_fertilization [mtdna_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)"

mtdna_final_fertilization_He_no.na <- mtdna_data[!is.na(mtdna_data$final_fertilization) & !is.na(mtdna_data$He),] #exclu

theme_update(plot.title = element_text(hjust = 0.5))

ggplot(mtdna_final_fertilization_He_no.na) + geom_boxplot(aes(x = final_fertilization, y = He)) +
  ggtitle("mtDNA: Fertilization Method vs. He") +
  xlab("Fertilization Method") + ylab("He") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

##Reproduction Mode##

#combine protandry, protogyny, and true hermaphroditism into "Hermpahrodites"

mtdna_data$final_reproductionmode <- NA

mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="protogyny"] <- "Hermaphrodite"
mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="protandry"] <- "Hermaphrodite"
mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite"

mtdna_final_reproductionmode_He_no.na <- mtdna_data[!is.na(mtdna_data$final_reproductionmode ) & !is.na(mtdna_data$He),] #exclude NA's in hermaphrodite and He columns

ggplot(mtdna_final_reproductionmode_He_no.na ) + geom_boxplot(aes(x = final_reproductionmode , y = He)) + #final reproduction mode & He box plot
  ggtitle("mtDNA: Reproduction Mode vs. He") +
  xlab("Reproduction Mode ") + ylab("He") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Plots: Numerical Data#####

##Maxlength##

mtdna_maxlength_He_no.na <- mtdna_data[!is.na(mtdna_data$maxlength) & !is.na(mtdna_data$He),] 

ggplot(mtdna_maxlength_He_no.na, aes(x= He, y= maxlength)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ggtitle("mtDNA: He vs. Max Length") +
  xlab("He") + ylab("Max Length") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

{
t.test(maxlength, He, mu=0, alt= ("two.sided"),conf.level=0.95,
         paired=FALSE, data="mtdna_maxlength_He_no.na")
  }


#shapiro.test(final_fertilization)
t.test(data = mtdna_maxlength_He_no.na, x= maxlength$He, unequal.var = TRUE)

##Fecundity##
    
mtdna_fecundity_He_no.na <- mtdna_data[!is.na(mtdna_data$fecundity_mean) & !is.na(mtdna_data$He),] 
 
ggplot(mtdna_fecundity_He_no.na, aes(y= fecundity_mean, x= He)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ggtitle("mtDNA: He vs. Fecundity Mean") +
  xlab("He") + ylab("Fecundity Mean") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))


####################################################################################

############### msat data set ############### 

#####Box Plots: Character Data#####

##Fertilization##

#convert "in brood pouch or similar structure" to internal fertilization

msat_data$final_fertilization <- NA

msat_data$final_fertilization [msat_data$fertilization =="external"]  <- "external"
msat_data$final_fertilization [msat_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
msat_data$final_fertilization [msat_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)"

msat_final_fertilization_He_no.na <- msat_data[!is.na(msat_data$final_fertilization) & !is.na(msat_data$He),] #exclu

ggplot(msat_final_fertilization_He_no.na) + geom_boxplot(aes(x = final_fertilization, y = He)) +
ggtitle("msat: Fertilization Method vs. He") +
  xlab("Fertilization Method") + ylab("He") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

##Reproduction Mode##

#combine protandry, protogyny, and true hermaphroditism into "Hermpahrodites"

msat_data$final_reproductionmode <- NA

msat_data$final_reproductionmode  [msat_data$reproductionmode =="dioecism"]  <- "Dioecious"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protogyny"] <- "Hermaphrodite"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protandry"] <- "Hermaphrodite"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite"

msat_final_reproductionmode_He_no.na <- msat_data[!is.na(msat_data$final_reproductionmode ) & !is.na(msat_data$He),] #exclude NA's in hermaphrodite and He columns

ggplot(msat_final_reproductionmode_He_no.na ) + geom_boxplot(aes(x = final_reproductionmode , y = He)) + #final reproduction mode & He box plot
  ggtitle("msat: Reproduction Mode vs. He") +
  xlab("Reproduction Mode ") + ylab("He") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Plots: Numerical Data#####

##Maxlength##

msat_maxlength_He_no.na <- msat_data[!is.na(msat_data$maxlength) & !is.na(msat_data$He),] 

ggplot(msat_maxlength_He_no.na, aes(x= He, y= maxlength)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) + 
  ggtitle("msat: He vs. Max Length") +
   xlab("He") + ylab("Max Length") +
    theme(
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(color="blue", size=14, face="bold"),
      axis.title.y = element_text(color="red", size=14, face="bold"))

t.test("maxlength","He", subset(maxlength_He_no.na), alternative = ("two.sided"), mu=0, paired=TRUE, conf.level=0.95, data("maxlength_He_no.na"))


##Fecundity##
       
msat_fecundity_He_no.na <- msat_data[!is.na(msat_data$fecundity_mean) & !is.na(msat_data$He),] 
       
ggplot(msat_fecundity_He_no.na, aes(y= fecundity_mean, x= He)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
    se=TRUE) +
  ggtitle("msat: He vs. Fecundity") +
  xlab("He") + ylab("Fecundity") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))