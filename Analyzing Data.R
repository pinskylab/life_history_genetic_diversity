################################################### Script for Analyzing Data  #######################################################

#analayzes data 

##########################################################################################################################################

######## Set-up ########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(ggplot2)

#read in data
mtdna_data <- read.csv("mtdna_final_data.csv", stringsAsFactors = FALSE) #read in 

######## Set-up ########

#combine protandry, protogyny, and true hermaphroditism into "Hermpahrodites"

mtdna_data$hermaphrodite <- NA

mtdna_data$hermaphrodite [mtdna_data$reproductionmode =="dioecism"]  <- "No"
mtdna_data$hermaphrodite [mtdna_data$reproductionmode =="protogyny"] <- "Yes"
mtdna_data$hermaphrodite [mtdna_data$reproductionmode =="protandry"] <- "Yes"
mtdna_data$hermaphrodite [mtdna_data$reproductionmode =="true hermaphroditism"] <- "Yes"

#####Box Plots: Character Data#####

##Fertilization##

boxplot(He ~ fertilization,data=mtdna_data, main="He vs. Fertilization",
        xlab="He", ylab="fertilization")

#includes boxplot of Na's
ggplot(subset(mtdna_data, !is.na(He), !is.na(fertilization)), aes(x=He, y=fertilization)) + geom_boxplot(aes(x = He, y = fertilization))

##Reproduction Mode##

hermaphrodite_He_no.na <- mtdna_data[!is.na(mtdna_data$hermaphrodite) & !is.na(mtdna_data$He),] #exclude NA's in hermaphrodite and He columns

ggplot(hermaphrodite_He_no.na ) + geom_boxplot(aes(x = hermaphrodite, y = He)) #hermaphrodite & He box plot

#reproductionmode_He_no.na <- mtdna_data[!is.na(mtdna_data$reproductionmode) & !is.na(mtdna_data$He),]

#ggplot(eproductionmode_He_no.na ) + geom_boxplot(aes(x = reproductionmode, y = He, )) #reproduction column 

#####Plots: Numerical Data#####

##Maxlength##

maxlength_He_no.na <- mtdna_data[!is.na(mtdna_data$maxlength) & !is.na(mtdna_data$He),] 

{
ggplot(maxlength_He_no.na, aes(y = maxlength, x = He) +
  geom_point () 
}

ggplot(maxlength_He_no.na, aes(x= He, y= maxlength)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE)  

#plot(y=maxlength_He_no.na$maxlength, x=maxlength_He_no.na$He, main = "He vs. Max Length",
 #   ylab = "maxlength", xlab = "He",
  #    pch = 10, frame = FALSE)
#abline(lm(maxlength ~ He, data = mtdna_data), col = "red")


##Fecundity##
    
fecundity_He_no.na <- mtdna_data[!is.na(mtdna_data$fecundity) & !is.na(mtdna_data$He),] 
 
ggplot(fecundity_He_no.na, aes(x= fecundity, y= He)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE)  

#plot(y=fecundity_He_no.na$fecundity, x=fecundity_He_no.na$He, main = "He vs. Fecundity",
 #    ylab = "fecundity", xlab = "He",
  #            pch = 10, frame = FALSE)
        # abline(lm(fecundity ~ He, data = mtdna_data), col = "red")
         
