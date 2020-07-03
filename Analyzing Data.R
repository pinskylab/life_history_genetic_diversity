################################################### Script for Analyzing Data  #######################################################

#analyzes data 

##########################################################################################################################################

######## Set-up ########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(dplyr)

#read in data
mtdna_data <- read.csv("mtdna_final_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("msat_final_data.csv", stringsAsFactors = FALSE) #read in 

############### mtDNA data set ############### 

#####Box Plots: Character Data#####

#Fertilization#

mtdna_data$final_fertilization <- NA #create new column to categorize fertilization methods

mtdna_data$final_fertilization [mtdna_data$fertilization =="external"]  <- "external"
mtdna_data$final_fertilization [mtdna_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
mtdna_data$final_fertilization [mtdna_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

mtdna_final_fertilization_He_no.na <- mtdna_data[!is.na(mtdna_data$final_fertilization) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest

theme_update(plot.title = element_text(hjust = 0.5)) #centered plot title

ggplot(mtdna_final_fertilization_He_no.na) + geom_boxplot(aes(x = final_fertilization, y = He)) + #final fertilization & He box plot
  ggtitle("mtDNA: Fertilization Method vs. He") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#Reproduction Mode#

#combine protandry, protogyny, and true hermaphroditism into "Hermpahrodites"

mtdna_data$final_reproductionmode <- NA #create new column to categorize reproduction mode

mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
mtdna_data$final_reproductionmode  [mtdna_data$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

mtdna_final_reproductionmode_He_no.na <- mtdna_data[!is.na(mtdna_data$final_reproductionmode ) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest

ggplot(mtdna_final_reproductionmode_He_no.na ) + geom_boxplot(aes(x = final_reproductionmode , y = He)) + #final reproduction mode & He box plot
  ggtitle("mtDNA: Reproduction Mode vs. He") + #add plot title
  xlab("Reproduction Mode ") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#Hermaphrodites#

mtdna_data$hermaphrodite_type <- NA #create new column to categorize hermaphrodite type

mtdna_data$hermaphrodite_type  [mtdna_data$reproductionmode =="dioecism"]  <- NA
mtdna_data$hermaphrodite_type  [mtdna_data$reproductionmode =="protogyny"] <- "protogyny"
mtdna_data$hermaphrodite_type  [mtdna_data$reproductionmode =="protandry"] <- "protandry"
mtdna_data$hermaphrodite_type  [mtdna_data$reproductionmode =="true hermaphroditism"] <- "true hermaphroditism"

mtdna_hermaphrodite_type_He_no.na <- mtdna_data[!is.na(mtdna_data$hermaphrodite_type ) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest

ggplot(mtdna_hermaphrodite_type_He_no.na ) + geom_boxplot(aes(x = hermaphrodite_type , y = He)) + #hermaphrodite type mode & He box plot
  ggtitle("mtDNA: Hermaphrodite Type vs. He") + #add plot title
  xlab("Hermaphrodite") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Scatter Plots: Numerical Data#####

#Max Length#

mtdna_maxlength_He_no.na <- mtdna_data[!is.na(mtdna_data$maxlength) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest

ggplot(mtdna_maxlength_He_no.na, aes(x= maxlength, y= He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("mtDNA: Max Length vs. He") + #add plot title
  xlab("Max Length") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#Fecundity Mean#
options(scipen = 999) #convert data from scientific notation to numeric

mtdna_fecundity_He_no.na <- mtdna_data[!is.na(mtdna_data$fecundity_mean) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest
 
ggplot(mtdna_fecundity_He_no.na, aes(y= He, x= fecundity_mean)) + #fecundity mean & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("mtDNA: Fecundity Mean vs. He") + #add plot title
  xlab("Fecundity Mean") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Density Plots: Numerical Data#####

#Max Length#
plot(density(mtdna_maxlength_He_no.na$maxlength), main="Max Length Density") #create density plot for max length
polygon(density(mtdna_maxlength_He_no.na$maxlength), main="Max Length Density", col="blue") #specify characteristics of plot

#Fecundity Mean#
plot(density(mtdna_fecundity_He_no.na$fecundity_mean), main="Fecundity Mean Density") #create density plot for fecundity mean
polygon(density(mtdna_fecundity_He_no.na$fecundity_mean), main="Fecundity Mean Density", col="blue") #specify characteristics of plot

#####T-Tests: Character Data#####

#Fertilization#

external <- mtdna_final_fertilization_He_no.na$He[mtdna_final_fertilization_He_no.na$final_fertilization=="external"]
internal <- mtdna_final_fertilization_He_no.na$He[mtdna_final_fertilization_He_no.na$final_fertilization=="internal (oviduct)"]

t.test(external, internal, var.equal=TRUE)

#Fecundity Mean#

fecunditymean.mtdna<-data.matrix(mtdna_fecundity_He_no.na$fecundity_mean) #create fecundity mean vector that contains matrix for fecundity mean converted from dataset 
He.mtdna<- data.matrix(mtdna_fecundity_He_no.na$He) #create He vector that contains matrix for He converted from dataset 

t.test(fecunditymean.mtdna, He.mtdna, var.equal=TRUE) #run a t-test using both newly created vectors 

#OR:: t.test(fecunditymean.mtdna, He., alt = "two.sided", mu = 0, paired = TRUE, var.equal = FALSE, conf.level = 0.95)
#OR:: t.test(fecunditymean.mtdna, He., alt = "two.sided", mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
#OR:: t.test(fecunditymean.mtdna, var.equal = FALSE, conf.level = 0.95)

#####Wilcoxon Tests: Numerical Data#####

#Max Length#
wilcox.test( mtdna_maxlength_He_no.na[ ,'maxlength'] , mtdna_maxlength_He_no.na[ , 'He'], paired=T) #run Wilcoxon test on max length & He


#Fecundity Mean#
wilcox.test( mtdna_fecundity_He_no.na[ ,'fecundity_mean'] , mtdna_fecundity_He_no.na[ , 'He'], paired=T) #run Wilcoxon test on max length & He

#####Shapiro-Wilk Tests: Numerical Data#####

#Max Length#
shapiro.test(mtdna_maxlength_He_no.na$maxlength) #run Shapiro-Wilk test on max length & He

#Fecundity Mean#
shapiro.test(mtdna_fecundity_He_no.na$fecundity_mean) #run Shapiro-Wilk test on max length & He


####################################################################################

############### msat Data Set ############### 

#####Box Plots: Character Data#####

#Fertilization#
msat_data$final_fertilization <- NA #create new column to categorize fertilization methods

msat_data$final_fertilization [msat_data$fertilization =="external"]  <- "external"
msat_data$final_fertilization [msat_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
msat_data$final_fertilization [msat_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

msat_final_fertilization_He_no.na <- msat_data[!is.na(msat_data$final_fertilization) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

ggplot(msat_final_fertilization_He_no.na) + geom_boxplot(aes(x = final_fertilization, y = He)) + #final fertilization & He box plot
ggtitle("msat: Fertilization Method vs. He") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#Reproduction Mode#
msat_data$final_reproductionmode <- NA #create new column to categorize reproduction mode

msat_data$final_reproductionmode  [msat_data$reproductionmode =="dioecism"]  <- "Dioecious"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

msat_final_reproductionmode_He_no.na <- msat_data[!is.na(msat_data$final_reproductionmode ) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

ggplot(msat_final_reproductionmode_He_no.na ) + geom_boxplot(aes(x = final_reproductionmode , y = He)) + #final reproduction mode & He box plot
  ggtitle("msat: Reproduction Mode vs. He") + #add plot title
  xlab("Reproduction Mode ") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#Hermaphrodites#

msat_data$hermaphrodite_type <- NA #create new column to categorize hermaphrodite type

msat_data$hermaphrodite_type  [msat_data$reproductionmode =="dioecism"]  <- NA
msat_data$hermaphrodite_type  [msat_data$reproductionmode =="protogyny"] <- "protogyny"
msat_data$hermaphrodite_type  [msat_data$reproductionmode =="protandry"] <- "protandry"
msat_data$hermaphrodite_type  [msat_data$reproductionmode =="true hermaphroditism"] <- "true hermaphroditism"

msat_hermaphrodite_type_He_no.na <- msat_data[!is.na(msat_data$hermaphrodite_type ) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

ggplot(msat_hermaphrodite_type_He_no.na) + geom_boxplot(aes(x = hermaphrodite_type , y = He)) + #hermaphrodite typemode & He box plot
  ggtitle("mtDNA: Hermaphrodite Type vs. He") + #add plot title
  xlab("Hermaphrodite") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Scatter Plots: Numerical Data#####

#Max Length#
msat_maxlength_He_no.na <- msat_data[!is.na(msat_data$maxlength) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

ggplot(msat_maxlength_He_no.na, aes(x= He, y= maxlength)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) + 
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("msat: Max Length vs. He") + #add plot title
   xlab("Max Length") + ylab("He") + #add axis labels
    theme(                                 #specifying characteristics of the plot 
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(color="blue", size=14, face="bold"),
      axis.title.y = element_text(color="red", size=14, face="bold"))

#Fecundity Mean#
msat_fecundity_He_no.na <- msat_data[!is.na(msat_data$fecundity_mean) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest
       
ggplot(msat_fecundity_He_no.na, aes(y= fecundity_mean, x= He)) + #fecundity mean& He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
    se=TRUE) +
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("msat: Fecundity vs. He") + #add plot title
  xlab("He") + ylab("Fecundity") + #add plot title
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Density Plots: Numerical Data#####

#Max Length#
plot(density(msat_maxlength_He_no.na$maxlength), main="msat: Max Length Density") #create density plot for max length
polygon(density(msat_maxlength_He_no.na$maxlength), main="msat: Max Length Density", col="blue") #specify characteristics of plot

#Fecundity Mean#
plot(density(msat_fecundity_He_no.na$fecundity_mean), main="msat: Fecundity Mean Density") #create density plot for fecundity mean
polygon(density(msat_fecundity_He_no.na$fecundity_mean), main="msat: Fecundity Mean Density", col="blue") #specify characteristics of plot

#####T-Tests: Numerical Data#####

#Max Length#

maxlength.msat<-data.matrix(msat_maxlength_He_no.na$maxlength) #create max length vector that contains matrix for max length converted from dataset 
He.msat<- data.matrix(msat_maxlength_He_no.na$He) #create He vector that contains matrix for He converted from dataset 

t.test(maxlength.msat, He.msat, var.equal=TRUE) #run a t-test using both newly created vectors 

#Fecundity Mean#

fecunditymean.msat<-data.matrix(msat_maxlength_He_no.na$fecundity_mean) #create fecundity mean vector that contains matrix for fecundity meanconverted from dataset 
He.msat<- data.matrix(msat_maxlength_He_no.na$He) #create He vector that contains matrix for He converted from dataset 

t.test(fecunditymean.msat, He.msat, var.equal=TRUE) #run a t-test using both newly created vectors 

#OR:: t.test(fecunditymean.msat, He.msat, alt = "two.sided", mu = 0, paired = TRUE, var.equal = FALSE, conf.level = 0.95)
#OR:: t.test(fecunditymean.msat, He.msat, alt = "two.sided", mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95)
#OR:: t.test(fecunditymean.msat, var.equal = FALSE, conf.level = 0.95)


#####Wilcoxon Tests: Numerical Data#####

#Max Length#
wilcox.test( msat_maxlength_He_no.na[ ,'maxlength'] , msat_maxlength_He_no.na[ , 'He'], paired=T) #run Wilcoxon test on max length & He

#Fecundity Mean#
wilcox.test( msat_fecundity_He_no.na[ ,'fecundity_mean'] , msat_fecundity_He_no.na[ , 'He'], paired=T) #run Wilcoxon test on fecundity & He


#####Shapiro-Wilk Tests: Numerical Data#####

shapiro.test(msat_maxlength_He_no.na$maxlength) #run Shapiro-Wilk test on max length & He

shapiro.test(msat_fecundity_He_no.na$fecundity_mean) #run Shapiro-Wilk test on fecundity mean & He

####################################################################################

############### Comparing mtDNA and msat Data ############### 

#####Box Plots: Character Data#####

#Fertilization#
{
ggarrange(
ggplot(mtdna_final_fertilization_He_no.na&msat_final_fertilization_He_no.na) + geom_boxplot(aes(x = final_fertilization, y = He)) +
  ggtitle("mtDNA: Fertilization Method vs. He") +
  xlab("Fertilization Method") + ylab("He") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold")) +

ggplot(msat_final_fertilization_He_no.na) + geom_boxplot(aes(x = final_fertilization, y = He)) + #create boxplot
  ggtitle("msat: Fertilization Method vs. He") +
  xlab("Fertilization Method") + ylab("He") +
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold")))
}

boxplot(ozone, ozone_norm, temp, temp_norm,
        main = "Multiple boxplots for comparision",
        at = c(1,2,4,5),
        names = c("ozone", "normal", "temp", "normal"),
        las = 2,
        col = c("orange","red"),
        border = "brown",
        horizontal = TRUE,
        notch = TRUE
)