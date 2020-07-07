################################################### Script for Analyzing Data  #######################################################

#analyzes data 

##########################################################################################################################################

######### Set-up #########

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

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

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

#Specific Reproduction Modes#

mtdna_data$specific.repro_mode <- NA #create new column to categorize hermaphrodite type

mtdna_data$specific.repro_mode  [mtdna_data$reproductionmode =="dioecism"]  <- "Dioecism"
mtdna_data$specific.repro_mode  [mtdna_data$reproductionmode =="protogyny"] <- "Protogyny"
mtdna_data$specific.repro_mode  [mtdna_data$reproductionmode =="protandry"] <- "Protandry"
mtdna_data$specific.repro_mode  [mtdna_data$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

mtdna_reproduction_type_He_no.na <- mtdna_data[!is.na(mtdna_data$specific.repro_mode ) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest

ggplot(mtdna_reproduction_type_He_no.na ) + geom_boxplot(aes(x = specific.repro_mode , y = He)) + #hermaphrodite type mode & He box plot
  ggtitle("mtDNA: Specific Reproduction Modes vs. He") + #add plot title
  xlab("Specific Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Scatter Plots: Numerical Data#####

#Max Length#

#w/ Rhincodon typus outlier
mtdna_maxlength_He_no.na <- mtdna_data[!is.na(mtdna_data$maxlength) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest

lm_eqn = function(x, y, df){ #set up formula for regression line equation
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == b %.% italic(x) + a,
                   list(a = format(coef(m)[[1]], digits = 2), 
                        b = format(coef(m)[[2]], digits = 2)))
  as.character(as.expression(eq));                 
}

ggplot(mtdna_maxlength_He_no.na, aes(x=maxlength, y = He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, formula = y ~ x) +
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("mtDNA: Max Length vs. He") + #add plot title
  xlab("Max Length") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
    annotate("text", x = 1200, y = 0.67, label = lm_eqn(mtdna_maxlength_He_no.na$maxlength, mtdna_maxlength_He_no.na$He, mtdna_maxlength_He_no.na), color="black", size = 5, parse=TRUE) #add regression line equation


#w/out Rhincodon typus outlier
mtdna_maxlength_He_no.na_nowhaleshark <- mtdna_maxlength_He_no.na[!(mtdna_maxlength_He_no.na$spp=="Rhincodon typus"),] #exclude Rhincodon typus

ggplot(mtdna_maxlength_He_no.na_nowhaleshark, aes(x= maxlength, y= He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("mtDNA: Max Length vs. He", subtitle= "(no Rhincodon typus)") + #add plot title
  xlab("Max Length") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate("text", x = 450, y = 0.67, label = lm_eqn(mtdna_maxlength_He_no.na_nowhaleshark$maxlength, mtdna_maxlength_He_no.na_nowhaleshark$He, mtdna_maxlength_He_no.na_nowhaleshark), color="black", size = 5, parse=TRUE) #add regression line equation

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
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate("text", x = 5000000, y = 0.67, label = lm_eqn(mtdna_fecundity_He_no.na$fecundity_mean, mtdna_fecundity_He_no.na$He, mtdna_fecundity_He_no.na), color="black", size = 5, parse=TRUE) #add regression line equation


#####Density Plots: Numerical Data#####

#Max Length#
plot(density(mtdna_maxlength_He_no.na$maxlength), main="Max Length Density") #create density plot for max length
polygon(density(mtdna_maxlength_He_no.na$maxlength), main="Max Length Density", col="blue") #specify characteristics of plot

#Fecundity Mean#
plot(density(mtdna_fecundity_He_no.na$fecundity_mean), main="Fecundity Mean Density") #create density plot for fecundity mean
polygon(density(mtdna_fecundity_He_no.na$fecundity_mean), main="Fecundity Mean Density", col="blue") #specify characteristics of plot

#####T-Tests: Character Data#####

#Fertilization#

external.mtdna <- mtdna_final_fertilization_He_no.na$He[mtdna_final_fertilization_He_no.na$final_fertilization=="external"]
internal.mtdna <- mtdna_final_fertilization_He_no.na$He[mtdna_final_fertilization_He_no.na$final_fertilization=="internal (oviduct)"]

fertilization_ttest.mtdna <- t.test(external.mtdna, internal.mtdna, var.equal=TRUE)

#Reproduction Mode#

dioecism.mtdna <- mtdna_final_reproductionmode_He_no.na$He[mtdna_final_reproductionmode_He_no.na$final_reproductionmode=="Dioecious"]
hermaphrodite.mtdna <- mtdna_final_reproductionmode_He_no.na$He[mtdna_final_reproductionmode_He_no.na$final_reproductionmode=="Hermaphrodite"]

reproductionmode_ttest.mtdna <- t.test(dioecism.mtdna, hermaphrodite.mtdna, var.equal=TRUE)

#Hermaphrodites: ANOVA TEST#

hermaphroditesanovamtdna <- aov(He ~ hermaphrodite_type, data = mtdna_hermaphrodite_type_He_no.na)

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

#Specific Reproduction Mode#

msat_data$specific.repro_mode <- NA #create new column to categorize reproduction type

msat_data$specific.repro_mode  [msat_data$reproductionmode =="dioecism"]  <- "Dioecism"
msat_data$specific.repro_mode  [msat_data$reproductionmode =="protogyny"] <- "Protogyny"
msat_data$specific.repro_mode  [msat_data$reproductionmode =="protandry"] <- "Protandry"
msat_data$specific.repro_mode  [msat_data$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

msat_reproduction_type_He_no.na <- msat_data[!is.na(msat_data$specific.repro_mode ) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

ggplot(msat_reproduction_type_He_no.na ) + geom_boxplot(aes(x = specific.repro_mode , y = He)) + #hermaphrodite type mode & He box plot
  ggtitle("msat: Specific Reproduction Modes vs. He") + #add plot title
  xlab("Specific Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#####Scatter Plots: Numerical Data#####

#Max Length#

#w/ Rhincodon typus outlier
msat_maxlength_He_no.na <- msat_data[!is.na(msat_data$maxlength) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

ggplot(msat_maxlength_He_no.na, aes(x= maxlength, y= He)) + #max length & He scatter plot
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
      axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate("text", x = 1200, y = 0.67, label = lm_eqn(msat_maxlength_He_no.na$maxlength, msat_maxlength_He_no.na$He, msat_maxlength_He_no.na), color="black", size = 5, parse=TRUE) #add regression line equation


#w/out Rhincodon typus outlier
msat_maxlength_He_no.na_nowhaleshark <- msat_maxlength_He_no.na[!(msat_maxlength_He_no.na$spp=="Rhincodon typus"),] #exclude Rhincodon typus

ggplot(msat_maxlength_He_no.na_nowhaleshark, aes(x= maxlength, y= He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("msat: Max Length vs. He", subtitle= "(no Rhincodon typus)") + #add plot title
  xlab("Max Length") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate("text", x = 450, y = 1, label = lm_eqn(msat_maxlength_He_no.na_nowhaleshark$maxlength, msat_maxlength_He_no.na_nowhaleshark$He, msat_maxlength_He_no.na_nowhaleshark), color="black", size = 5, parse=TRUE) #add regression line equation


#Fecundity Mean#
msat_fecundity_He_no.na <- msat_data[!is.na(msat_data$fecundity_mean) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest
       
ggplot(msat_fecundity_He_no.na, aes(x= fecundity_mean, y= He)) + #fecundity mean& He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
    se=TRUE) +
  ylim(0,1)+
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("msat: Fecundity Mean vs. He") + #add plot title
  xlab("Fecundity Mean") + ylab("He") + #add plot title
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate("text", x = 25000000, y = 0.956, label = lm_eqn(msat_fecundity_He_no.na$fecundity_mean, msat_fecundity_He_no.na$He, msat_fecundity_He_no.na), color="black", size = 5, parse=TRUE) #add regression line equation

#####Density Plots: Numerical Data#####

#Max Length#
plot(density(msat_maxlength_He_no.na$maxlength), main="msat: Max Length Density") #create density plot for max length
polygon(density(msat_maxlength_He_no.na$maxlength), main="msat: Max Length Density", col="blue") #specify characteristics of plot

#Fecundity Mean#
plot(density(msat_fecundity_He_no.na$fecundity_mean), main="msat: Fecundity Mean Density") #create density plot for fecundity mean
polygon(density(msat_fecundity_He_no.na$fecundity_mean), main="msat: Fecundity Mean Density", col="blue") #specify characteristics of plot

#####T-Tests: Character Data#####

#Fertilization#

external.msat <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="external"]
internal.msat <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="internal (oviduct)"]

fertilization_ttest.msat <- t.test(external.msat, internal.msat, var.equal=TRUE)

#Reproduction Mode#

dioecism.msat <- msat_final_reproductionmode_He_no.na$He[msat_final_reproductionmode_He_no.na$final_reproductionmode=="Dioecious"]
hermaphrodite.msat <- msat_final_reproductionmode_He_no.na$He[msat_final_reproductionmode_He_no.na$final_reproductionmode=="Hermaphrodite"]

reproductionmode_ttest.msat <- t.test(dioecism.msat, hermaphrodite.msat, var.equal=TRUE)

#Hermaphrodites: ANOVA TEST#

hermaphroditesanovamsat <- aov(He ~ hermaphrodite_type, data = msat_hermaphrodite_type_He_no.na)

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

final_fertilization_all = merge(mtdna_final_fertilization_He_no.na, msat_final_fertilization_He_no.na, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_fertilization_all$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_fertilization_all$markertype [final_fertilization_all$file == "mtdna101"]  <- "mtDNA"
final_fertilization_all$markertype [final_fertilization_all$file == "mtdna102"]  <- "mtDNA"
final_fertilization_all$markertype [final_fertilization_all$file == "mtdna103"]  <- "mtDNA"
final_fertilization_all$markertype [final_fertilization_all$file == "msats000"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats001"]  <- "msat"
final_fertilization_all$markertype [final_fertilization_all$file == "msats002"]  <- "msat"
final_fertilization_all$markertype [final_fertilization_all$file == "msats200"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats201"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats100"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats101"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats301"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats302"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats303"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats304"]  <- "msat" 
final_fertilization_all$markertype [final_fertilization_all$file == "msats305"]  <- "msat"
final_fertilization_all$markertype [final_fertilization_all$file ==	"ppdat"]  <- "msat" 

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

ggplot(final_fertilization_all) + geom_boxplot(aes(x = markertype, y = He, fill=final_fertilization)) + #final fertilization & He box plot
  ggtitle("Fertilization Method vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("darkcyan", "turquoise2"))

#Reproduction Mode#
reproductionmode_all = merge(msat_final_reproductionmode_He_no.na, mtdna_final_reproductionmode_He_no.na, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

reproductionmode_all$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
reproductionmode_all$markertype [reproductionmode_all$file == "mtdna101"]  <- "mtDNA"
reproductionmode_all$markertype [reproductionmode_all$file == "mtdna102"]  <- "mtDNA"
reproductionmode_all$markertype [reproductionmode_all$file == "mtdna103"]  <- "mtDNA"
reproductionmode_all$markertype [reproductionmode_all$file == "msats000"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats001"]  <- "msat"
reproductionmode_all$markertype [reproductionmode_all$file == "msats002"]  <- "msat"
reproductionmode_all$markertype [reproductionmode_all$file == "msats200"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats201"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats100"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats101"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats301"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats302"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats303"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats304"]  <- "msat" 
reproductionmode_all$markertype [reproductionmode_all$file == "msats305"]  <- "msat"
reproductionmode_all$markertype [reproductionmode_all$file ==	"ppdat"]  <- "msat" 

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

ggplot(reproductionmode_all) + geom_boxplot(aes(x = markertype, y = He, fill=final_reproductionmode)) + #final fertilization & He box plot
  ggtitle("Reproduction Mode vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("darkcyan", "turquoise2"))


#Specific Reproduction Mode#
specificreproductionmode_all = merge(msat_reproduction_type_He_no.na, mtdna_reproduction_type_He_no.na, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

specificreproductionmode_all$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "mtdna101"]  <- "mtDNA"
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "mtdna102"]  <- "mtDNA"
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "mtdna103"]  <- "mtDNA"
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats000"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats001"]  <- "msat"
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats002"]  <- "msat"
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats200"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats201"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats100"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats101"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats301"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats302"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats303"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats304"]  <- "msat" 
specificreproductionmode_all$markertype [specificreproductionmode_all$file == "msats305"]  <- "msat"
specificreproductionmode_all$markertype [specificreproductionmode_all$file ==	"ppdat"]  <- "msat" 

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

ggplot(specificreproductionmode_all) + geom_boxplot(aes(x = markertype, y = He, fill=specific.repro_mode)) + #final fertilization & He box plot
  ggtitle("Specific Reproduction Modes vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Specific Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_brewer(palette="Blues")

