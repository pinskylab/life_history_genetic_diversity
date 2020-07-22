################################################### Script for Analyzing Data  #######################################################

#analyzes mtDNA & msat datasets

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

mtdna_maxlength_He_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlength_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlength_He_no.na$logtransform.maxlength <- log10(mtdna_maxlength_He_no.na$maxlength)
}

lm_eqn = function(x, y, df){ #set up formula for regression line equation
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == b %.% italic(x) + a,
                   list(a = format(coef(m)[[1]], digits = 2), 
                        b = format(coef(m)[[2]], digits = 2)))
  as.character(as.expression(eq));                 
}

ggplot(mtdna_maxlength_He_no.na, aes(x=logtransform.maxlength, y = He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, formula = y ~ x) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("mtDNA: Max Length vs. He") + #add plot title
  xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 2.4, y = .807, label = lm_eqn(mtdna_maxlength_He_no.na$logtransform.maxlength, mtdna_maxlength_He_no.na$He, mtdna_maxlength_He_no.na), 
           color="black", size = 5, parse=TRUE, alpha=0.80)  #add regression line equation

#w/out Rhincodon typus outlier
mtdna_maxlength_He_no.na_nowhaleshark <- mtdna_maxlength_He_no.na[!(mtdna_maxlength_He_no.na$spp=="Rhincodon typus"),] #exclude Rhincodon typus

mtdna_maxlength_He_no.na_nowhaleshark$logtransform.maxlength_noRT <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlength_He_no.na_nowhaleshark)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlength_He_no.na_nowhaleshark$logtransform.maxlength_noRT <- log10(mtdna_maxlength_He_no.na_nowhaleshark$maxlength)
}

ggplot(mtdna_maxlength_He_no.na_nowhaleshark, aes(x= logtransform.maxlength_noRT, y= He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("mtDNA: Max Length vs. He", subtitle= "(no Rhincodon typus)") + #add plot title
  xlab("Max Length (loc(cm))") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.5, y = 0.65, label = lm_eqn(mtdna_maxlength_He_no.na_nowhaleshark$logtransform.maxlength_noRT, mtdna_maxlength_He_no.na_nowhaleshark$He, mtdna_maxlength_He_no.na_nowhaleshark), 
           color="black", size = 5, parse=TRUE, alpha=0.80) #add regression line equation

#Fecundity Mean#
 mtdna_fecundity_He_no.na <- mtdna_data[!is.na(mtdna_data$fecundity_mean) & !is.na(mtdna_data$He),] #create new table that excludes NA's from columns of interest
 
mtdna_fecundity_He_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_fecundity_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_fecundity_He_no.na$logtransform.fecundity <- log10(mtdna_fecundity_He_no.na$fecundity_mean)
}

ggplot(mtdna_fecundity_He_no.na, aes(x= logtransform.fecundity, y= He)) + #fecundity mean & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("mtDNA: Fecundity Mean vs. He") + #add plot title
  xlab("Fecundity Mean (log)") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 4, y = 0.57, label = lm_eqn(mtdna_fecundity_He_no.na$logtransform.fecundity, mtdna_fecundity_He_no.na$He, mtdna_fecundity_He_no.na), 
           color="black", size = 5, parse=TRUE, alpha=0.80) #add regression line equation


#####Density Plots: Numerical Data#####

#Max Length#
plot(density(mtdna_maxlength_He_no.na$maxlength), main="Max Length Density") #create density plot for max length
polygon(density(mtdna_maxlength_He_no.na$maxlength), main="Max Length Density", col="blue") #specify characteristics of plot

#Fecundity Mean#
plot(density(mtdna_fecundity_He_no.na$fecundity_mean), main="Fecundity Mean Density") #create density plot for fecundity mean
polygon(density(mtdna_fecundity_He_no.na$fecundity_mean), main="Fecundity Mean Density", col="blue") #specify characteristics of plot

#####T-Tests: Character Data#####

#Fertilization#

external.mtdna <- mtdna_final_fertilization_He_no.na$He[mtdna_final_fertilization_He_no.na$final_fertilization=="external"] #create vector for one aspect of t-test
internal.mtdna <- mtdna_final_fertilization_He_no.na$He[mtdna_final_fertilization_He_no.na$final_fertilization=="internal (oviduct)"] #create vector

fertilization_ttest.mtdna <- t.test(external.mtdna, internal.mtdna, var.equal=TRUE) #combine created vectors & perform t-test

#Reproduction Mode#

dioecism.mtdna <- mtdna_final_reproductionmode_He_no.na$He[mtdna_final_reproductionmode_He_no.na$final_reproductionmode=="Dioecious"] #create vector for one aspect of t-test
hermaphrodite.mtdna <- mtdna_final_reproductionmode_He_no.na$He[mtdna_final_reproductionmode_He_no.na$final_reproductionmode=="Hermaphrodite"] #create vector

reproductionmode_ttest.mtdna <- t.test(dioecism.mtdna, hermaphrodite.mtdna, var.equal=TRUE) #combine created vectors & perform t-test

#Specific Reproduction Modes: ANOVA TEST#

specific.repro_modeanovamtdna <- aov(He ~ specific.repro_mode, data = mtdna_reproduction_type_He_no.na) #perform anova test

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

msat_maxlength_He_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(msat_maxlength_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlength_He_no.na$logtransform.maxlength <- log10(msat_maxlength_He_no.na$maxlength)
}

ggplot(msat_maxlength_He_no.na, aes(x= logtransform.maxlength, y= He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) + 
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("msat: Max Length vs. He") + #add plot title
   xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
    theme(                                 #specifying characteristics of the plot 
      plot.title = element_text(size=14, face="bold"),
      axis.title.x = element_text(color="blue", size=14, face="bold"),
      axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom = "label", x = 1.8, y = 0.83, label = lm_eqn(msat_maxlength_He_no.na$logtransform.maxlength, msat_maxlength_He_no.na$He, msat_maxlength_He_no.na), 
           color="black", size = 5, parse=TRUE, alpha = 0.8) #add regression line equation


#w/out Rhincodon typus outlier
msat_maxlength_He_no.na_nowhaleshark <- msat_maxlength_He_no.na[!(msat_maxlength_He_no.na$spp=="Rhincodon typus"),] #exclude Rhincodon typus

msat_maxlength_He_no.na_nowhaleshark$logtransform.maxlength_noRT <- NA #add column to do a log transformation for max length

for (i in 1:nrow(msat_maxlength_He_no.na_nowhaleshark)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlength_He_no.na_nowhaleshark$logtransform.maxlength_noRT <- log10(msat_maxlength_He_no.na_nowhaleshark$maxlength)
}

ggplot(msat_maxlength_He_no.na_nowhaleshark, aes(x= logtransform.maxlength_noRT, y= He)) + #max length & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("msat: Max Length vs. He", subtitle= "(no Rhincodon typus)") + #add plot title
  xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.9, y = 0.83, label = lm_eqn(msat_maxlength_He_no.na_nowhaleshark$logtransform.maxlength_noRT, msat_maxlength_He_no.na_nowhaleshark$He, msat_maxlength_He_no.na_nowhaleshark), 
           color="black", size = 5, parse=TRUE, alpha = 0.8) #add regression line equation


#Fecundity Mean#
msat_fecundity_He_no.na <- msat_data[!is.na(msat_data$fecundity_mean) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

msat_fecundity_He_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(msat_fecundity_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_fecundity_He_no.na$logtransform.fecundity <- log10(msat_fecundity_He_no.na$fecundity_mean)
}
       
ggplot(msat_fecundity_He_no.na, aes(x= logtransform.fecundity, y= He)) + #fecundity mean & He scatter plot
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
    se=TRUE) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("msat: Fecundity Mean vs. He") + #add plot title
  xlab("Fecundity Mean (log)") + ylab("He") + #add plot title
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 3.5, y = 0.83, label = lm_eqn(msat_fecundity_He_no.na$logtransform.fecundity, msat_fecundity_He_no.na$He, msat_fecundity_He_no.na), 
           color="black", size = 5, parse=TRUE, alpha = 0.8) #add regression line equation

#####Density Plots: Numerical Data#####

#Max Length#
plot(density(msat_maxlength_He_no.na$maxlength), main="msat: Max Length Density") #create density plot for max length
polygon(density(msat_maxlength_He_no.na$maxlength), main="msat: Max Length Density", col="blue") #specify characteristics of plot

#Fecundity Mean#
plot(density(msat_fecundity_He_no.na$fecundity_mean), main="msat: Fecundity Mean Density") #create density plot for fecundity mean
polygon(density(msat_fecundity_He_no.na$fecundity_mean), main="msat: Fecundity Mean Density", col="blue") #specify characteristics of plot

#####T-Tests: Character Data#####

#Fertilization#

external.msat <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="external"] #create vector for one aspect of t-test
internal.msat <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="internal (oviduct)"] #create vector

fertilization_ttest.msat <- t.test(external.msat, internal.msat, var.equal=TRUE) #combine created vectors & perform t-test

#Reproduction Mode#

dioecism.msat <- msat_final_reproductionmode_He_no.na$He[msat_final_reproductionmode_He_no.na$final_reproductionmode=="Dioecious"] #create vector for one aspect of t-test
hermaphrodite.msat <- msat_final_reproductionmode_He_no.na$He[msat_final_reproductionmode_He_no.na$final_reproductionmode=="Hermaphrodite"] #create vector

reproductionmode_ttest.msat <- t.test(dioecism.msat, hermaphrodite.msat, var.equal=TRUE) #combine created vectors & perform t-test

#Specific Reproduction Modes: ANOVA TEST#

specific.repro_modeanovamsat <- aov(He ~ specific.repro_mode, data = msat_reproduction_type_He_no.na) #perform anova test

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

ggplot(final_fertilization_all) + geom_boxplot(aes(x = final_fertilization, y = He, fill=markertype)) + #final fertilization & He box plot
  ggtitle("Fertilization Method vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("turquoise2", "darkcyan"))
 
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

ggplot(reproductionmode_all) + geom_boxplot(aes(x = final_reproductionmode, y = He, fill= markertype)) + #final fertilization & He box plot
  ggtitle("Reproduction Mode vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("turquoise2", "darkcyan"))


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

ggplot(specificreproductionmode_all) + geom_boxplot(aes(x = specific.repro_mode, y = He, fill=markertype)) + #final fertilization & He box plot
  ggtitle("Specific Reproduction Modes vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Specific Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("turquoise2", "darkcyan")) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5))


#####Scatter Plots: Numerical Data#####

#Max Length#

#w/ Rhincodon typus outlier
final_maxlength_all = merge(mtdna_maxlength_He_no.na, msat_maxlength_He_no.na, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_maxlength_all$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_maxlength_all$markertype [final_maxlength_all$file == "mtdna101"]  <- "mtDNA"
final_maxlength_all$markertype [final_maxlength_all$file == "mtdna102"]  <- "mtDNA"
final_maxlength_all$markertype [final_maxlength_all$file == "mtdna103"]  <- "mtDNA"
final_maxlength_all$markertype [final_maxlength_all$file == "msats000"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats001"]  <- "msat"
final_maxlength_all$markertype [final_maxlength_all$file == "msats002"]  <- "msat"
final_maxlength_all$markertype [final_maxlength_all$file == "msats200"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats201"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats100"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats101"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats301"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats302"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats303"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats304"]  <- "msat" 
final_maxlength_all$markertype [final_maxlength_all$file == "msats305"]  <- "msat"
final_maxlength_all$markertype [final_maxlength_all$file ==	"ppdat"]  <- "msat" 

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

ggplot(final_maxlength_all, aes(x=logtransform.maxlength, y=He, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
             se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. He", subtitle = "(w/ Rhincodon typus) msat vs. mtDNA") + #add plot title
  xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.5, y = 0.55, label = lm_eqn(msat_maxlength_He_no.na$logtransform.maxlength, msat_maxlength_He_no.na$He, msat_maxlength_He_no.na), 
           color="skyblue3", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 2, y = 0.91, label = lm_eqn(mtdna_maxlength_He_no.na$logtransform.maxlength, mtdna_maxlength_He_no.na$He, mtdna_maxlength_He_no.na), 
           color="blue", size = 5, parse=TRUE, alpha = 0.8) + #add regression line equation
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)

#w/out Rhincodon typus outlier
final_maxlength.nowhaleshark_all = merge(msat_maxlength_He_no.na_nowhaleshark, mtdna_maxlength_He_no.na_nowhaleshark, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_maxlength.nowhaleshark_all$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "mtdna101"]  <- "mtDNA"
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "mtdna102"]  <- "mtDNA"
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "mtdna103"]  <- "mtDNA"
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats000"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats001"]  <- "msat"
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats002"]  <- "msat"
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats200"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats201"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats100"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats101"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats301"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats302"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats303"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats304"]  <- "msat" 
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file == "msats305"]  <- "msat"
final_maxlength.nowhaleshark_all$markertype [final_maxlength.nowhaleshark_all$file ==	"ppdat"]  <- "msat" 

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

ggplot(final_maxlength.nowhaleshark_all, aes(x=logtransform.maxlength_noRT, y=He, col=markertype, shape=markertype)) + #max length w/out Rhincodon typus & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. He", subtitle = "(no Rhincodon typus) msat vs. mtDNA") + #add plot title
  xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.4, y = 0.62, label = lm_eqn(msat_maxlength_He_no.na_nowhaleshark$logtransform.maxlength, msat_maxlength_He_no.na_nowhaleshark$He, msat_maxlength_He_no.na_nowhaleshark), 
           color="skyblue2", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 2, y = 0.90, label = lm_eqn(mtdna_maxlength_He_no.na_nowhaleshark$logtransform.maxlength, mtdna_maxlength_He_no.na_nowhaleshark$He, mtdna_maxlength_He_no.na_nowhaleshark), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_color_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)

#Fecundity Mean#

final_fecunditymean_all = merge(msat_fecundity_He_no.na, mtdna_fecundity_He_no.na, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_fecunditymean_all$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "mtdna101"]  <- "mtDNA"
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "mtdna102"]  <- "mtDNA"
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "mtdna103"]  <- "mtDNA"
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats000"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats001"]  <- "msat"
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats002"]  <- "msat"
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats200"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats201"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats100"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats101"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats301"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats302"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats303"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats304"]  <- "msat" 
final_fecunditymean_all$markertype [final_fecunditymean_all$file == "msats305"]  <- "msat"
final_fecunditymean_all$markertype [final_fecunditymean_all$file ==	"ppdat"]  <- "msat" 

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

ggplot(final_fecunditymean_all, aes(x=logtransform.fecundity, y=He, col=markertype, shape=markertype)) + #fecundity mean & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Fecundity Mean vs. He", subtitle = "msat vs. mtDNA") + #add plot title
  xlab("Fecundity Mean (log)") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 2.8, y = 0.9, label = lm_eqn(msat_fecundity_He_no.na$logtransform.fecundity, msat_fecundity_He_no.na$He, msat_fecundity_He_no.na), 
           color="skyblue2", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 3.4, y = 0.5, label = lm_eqn(mtdna_fecundity_He_no.na$logtransform.fecundity, mtdna_fecundity_He_no.na$He, mtdna_fecundity_He_no.na), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_color_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)


#####Anova Test: Combined Marker Character Data#####

#Fertilization#

fertilizationanova.all <- aov(He ~ final_fertilization * markertype, data = final_fertilization_all) #perform anova test for combined data

#Reproduction Mode#

reproductionmodeanova.all <- aov(He ~ final_reproductionmode * markertype, data = reproductionmode_all) #perform anova test for combined data

#Specific Reproduction Modes: ANOVA#

specific.repro_modeanova.all <- aov(He ~ specific.repro_mode * markertype, data = specificreproductionmode_all) #perform anova test for combined data
TukeyHSD(specific.repro_modeanova.all) #perform TukeyHSD to see full table of results

####################################################################################

############### mtDNA: Statistical Models ############### 

#####Numerical Data#####

#Max Length#

#w/ Rhincodon typus

mtdna_maxlength_He_no.na$success <- round(mtdna_maxlength_He_no.na$He*mtdna_maxlength_He_no.na$n) #create column of successes
mtdna_maxlength_He_no.na$failure <- round((1-mtdna_maxlength_He_no.na$He)*mtdna_maxlength_He_no.na$n) #create column of failures

ggplot(mtdna_maxlength_He_no.na, aes(x=maxlength, y=He)) + #max length vs. He scatter plot with statistical model
  geom_point(shape=1) +
  stat_smooth(method="glm",  method.args = list(family = "binomial"), formula= cbind(mtdna_maxlength_He_no.na$success, mtdna_maxlength_He_no.na$failure)~x) + #add generalized liner model line
  ggtitle("mtDNA:Max Length vs. He") + #add plot title
  xlab("Max Length (cm)") + ylab("He") + #add axis labels
  ylim(0,1)+                              #create limits
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#w/out Rhincodon typus

mtdna_maxlength_He_no.na_nowhaleshark$success <- round(mtdna_maxlength_He_no.na_nowhaleshark$He*mtdna_maxlength_He_no.na_nowhaleshark$n) #create column of successes
mtdna_maxlength_He_no.na_nowhaleshark$failure <- round((1-mtdna_maxlength_He_no.na_nowhaleshark$He)*mtdna_maxlength_He_no.na_nowhaleshark$n) #create column of failures

ggplot(mtdna_maxlength_He_no.na_nowhaleshark, aes(x=maxlength, y=He)) + #max length w/out Rhincodon tyous vs. He scatter plot with statistical model
  geom_point(shape=1) +
  stat_smooth(method="glm",  method.args = list(family = "binomial"), formula= cbind(mtdna_maxlength_He_no.na_nowhaleshark$success, mtdna_maxlength_He_no.na_nowhaleshark$failure)~x) + #add generalized liner model line
  ggtitle("mtDNA: Max Length vs. He", subtitle= "(no Rhincodon typus)") + #add plot  title
  xlab("Max Length (cm)") + ylab("He") + #add axis labels
  ylim(0,1)+                              #create limits
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#Fecundity Mean#

mtdna_fecundity_He_no.na$success <- round(mtdna_fecundity_He_no.na$He*mtdna_fecundity_He_no.na$n) #create column of successes
mtdna_fecundity_He_no.na$failure <- round((1-mtdna_fecundity_He_no.na$He)*mtdna_fecundity_He_no.na$n) #create column of failures

ggplot(mtdna_fecundity_He_no.na, aes(x=fecundity_mean, y=He)) + #fecundity mean vs. He scatter plot with statistical model
  geom_point(shape=1) +
  stat_smooth(method="glm",  method.args = list(family = "binomial"), formula= cbind(mtdna_fecundity_He_no.na$success, mtdna_fecundity_He_no.na$failure)~x) + #add generalized liner model line
  ggtitle("mtDNA: Fecundity Mean vs. He") + #add plot title
  xlab("Fecundity Mean") + ylab("He") + #add axis labels
  ylim(0,1)+                              #create limits
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

############### msat: Statistical Models ############### 

#####Numerical Data#####

#Max Length#

#w/ Rhincodon typus

msat_maxlength_He_no.nan <- msat_maxlength_He_no.na[!is.na(msat_maxlength_He_no.na$n),] #exclude NA
msat_maxlength_He_no.nan <- msat_maxlength_He_no.nan[!(msat_maxlength_He_no.nan$n=="25-32"),]  #exclude n=25-32

msat_maxlength_He_no.nan$success <- round(as.numeric(msat_maxlength_He_no.nan$He)*(as.numeric(msat_maxlength_He_no.nan$n))) #create column of successes
msat_maxlength_He_no.nan$failure <- round(as.numeric((1-msat_maxlength_He_no.nan$He)*(as.numeric(msat_maxlength_He_no.nan$n)))) #create column of failures

ggplot(msat_maxlength_He_no.nan, aes(x=maxlength, y=He)) + #max length vs. He scatter plot with statistical model
  geom_point(shape=1) +
  stat_smooth(method="glm",  method.args = list(family = "binomial"), formula= cbind(msat_maxlength_He_no.nan$success, msat_maxlength_He_no.nan$failure)~x) + #add generalized liner model line
  ggtitle("msat: Max Length vs. He") + #add plot title
  xlab("Max Length (cm)") + ylab("He") + #add axis labels
  ylim(0,1)+                             #create limits
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#w/out Rhincodon typus
msat_maxlength_He_no.na_nowhalesharkn <- msat_maxlength_He_no.na_nowhaleshark[!is.na(msat_maxlength_He_no.na_nowhaleshark$n),] #exclude NA
msat_maxlength_He_no.na_nowhalesharkn <- msat_maxlength_He_no.na_nowhalesharkn[!(msat_maxlength_He_no.na_nowhalesharkn$n=="25-32"),] #exclude n=25-32

msat_maxlength_He_no.na_nowhalesharkn$success <- round(as.numeric(msat_maxlength_He_no.na_nowhalesharkn$He)*(as.numeric(msat_maxlength_He_no.na_nowhalesharkn$n))) #create column of successes
msat_maxlength_He_no.na_nowhalesharkn$failure <- round(as.numeric(1-msat_maxlength_He_no.na_nowhalesharkn$He)*(as.numeric(msat_maxlength_He_no.na_nowhalesharkn$n))) #create column of failures

ggplot(msat_maxlength_He_no.na_nowhalesharkn, aes(x=maxlength, y=He)) + #max length w/out Rhincodon typus vs. He scatter plot with statistical model
  geom_point(shape=1) +
  stat_smooth(method="glm",  method.args = list(family = "binomial"), formula= cbind(msat_maxlength_He_no.na_nowhalesharkn$success, msat_maxlength_He_no.na_nowhalesharkn$failure)~x) + #add generalized liner model line
  ggtitle("msat: Max Length vs. He", subtitle= "(no Rhincodon typus)") + #add plot title
  xlab("Max Length (cm)") + ylab("He") + #add axis labels
  ylim(0,1)+                              #create limits
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

#Fecundity Mean#

msat_fecundity_He_no.nan <- msat_fecundity_He_no.na[!is.na(msat_fecundity_He_no.na$n),] #exclude NA 
msat_fecundity_He_no.nan <- msat_fecundity_He_no.nan[!(msat_fecundity_He_no.nan$n=="25-32"),] #exclude n=25-32

msat_fecundity_He_no.nan$success <- round(as.numeric(msat_fecundity_He_no.nan$He)*(as.numeric(msat_fecundity_He_no.nan$n))) #create column of successes
msat_fecundity_He_no.nan$failure <- round(as.numeric(1-msat_fecundity_He_no.nan$He)*(as.numeric(msat_fecundity_He_no.nan$n))) #create column of failures

ggplot(msat_fecundity_He_no.nan , aes(x=fecundity_mean, y=He)) + #fecundity mean vs. He scatter plot with statistical model
  geom_point(shape=1) +
  stat_smooth(method="glm",  method.args = list(family = "binomial"), formula= cbind(msat_fecundity_He_no.nan$success, msat_fecundity_He_no.nan$failure)~x) + #add generalized liner model line
  ggtitle("msat: Fecundity Mean vs. He") + #add plot title
  xlab("Fecundity Mean") + ylab("He") + #add axis labels
  ylim(0,1)+                              #create limits
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))

############### Combined Statistical Models ############### 

#Max Length#

#w/ Rhincodon typus

final_maxlength_alln <- final_maxlength_all[!is.na(final_maxlength_all$n),] #exclude NA 
final_maxlength_alln <- final_maxlength_alln[!(final_maxlength_alln$n=="25-32"),] #exclude n=25-32

final_maxlength_alln$success <- round(as.numeric(final_maxlength_alln$He)*(as.numeric(final_maxlength_alln$n))) #create column of successes
final_maxlength_alln$failure <- round(as.numeric(1-final_maxlength_alln$He)*(as.numeric(final_maxlength_alln$n))) #create column of failures

ggplot(final_maxlength_alln, aes(x=maxlength, y=He, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(data= msat_maxlength_He_no.nan, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), 
              formula= cbind(msat_maxlength_He_no.nan$success, msat_maxlength_He_no.nan$failure)~x, inherit.aes = FALSE,   #Add linear regression line outline
              se=TRUE, color = "black", size = 2, fill = NA) + 
  stat_smooth(data= msat_maxlength_He_no.nan, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), #Add linear regression line
              formula= cbind(msat_maxlength_He_no.nan$success, msat_maxlength_He_no.nan$failure)~x, inherit.aes = FALSE, color= "skyblue2") +
  geom_smooth(data= mtdna_maxlength_He_no.na, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), 
              formula= cbind(mtdna_maxlength_He_no.na$success, mtdna_maxlength_He_no.na$failure)~x, inherit.aes=FALSE,  #Add linear regression line outline
              se=TRUE, color = "black", size = 2, fill = NA) + 
  stat_smooth(data= mtdna_maxlength_He_no.na, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), #Add linear regression line
              formula= cbind(mtdna_maxlength_He_no.na$success, mtdna_maxlength_He_no.na$failure)~x, inherit.aes=FALSE, color= "blue") +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. He", subtitle = "(w/ Rhincodon typus) msat vs. mtDNA") + #add plot title
  xlab("Max Length (cm)") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE) 

#w/out Rhincodon typus

final_maxlength.nowhaleshark_alln <- final_maxlength.nowhaleshark_all[!is.na(final_maxlength.nowhaleshark_all$n),] #exclude NA 
final_maxlength.nowhaleshark_alln <- final_maxlength.nowhaleshark_alln[!(final_maxlength.nowhaleshark_alln$n=="25-32"),] #exclude n=25-32

final_maxlength.nowhaleshark_alln$success <- round(as.numeric(final_maxlength.nowhaleshark_alln$He)*(as.numeric(final_maxlength.nowhaleshark_alln$n))) #create column of successes
final_maxlength.nowhaleshark_alln$failure <- round(as.numeric(1-final_maxlength.nowhaleshark_alln$He)*(as.numeric(final_maxlength.nowhaleshark_alln$n))) #create column of failures

ggplot(final_maxlength.nowhaleshark_alln, aes(x=maxlength, y=He, col=markertype, shape=markertype)) + #max length w/out Rhincodon typus & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(data= msat_maxlength_He_no.na_nowhalesharkn, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), 
              formula= cbind(msat_maxlength_He_no.na_nowhalesharkn$success, msat_maxlength_He_no.na_nowhalesharkn$failure)~x, inherit.aes = FALSE,   #Add linear regression line outline
              se=TRUE, color = "black", size = 2, fill = NA) + 
  stat_smooth(data= msat_maxlength_He_no.na_nowhalesharkn, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), #Add linear regression line
              formula= cbind(msat_maxlength_He_no.na_nowhalesharkn$success, msat_maxlength_He_no.na_nowhalesharkn$failure)~x, inherit.aes = FALSE, color= "skyblue2") +
  geom_smooth(data= mtdna_maxlength_He_no.na_nowhaleshark, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), 
              formula= cbind(mtdna_maxlength_He_no.na_nowhaleshark$success, mtdna_maxlength_He_no.na_nowhaleshark$failure)~x, inherit.aes=FALSE,  #Add linear regression line outline
              se=TRUE, color = "black", size = 2, fill = NA) + 
  stat_smooth(data= mtdna_maxlength_He_no.na_nowhaleshark, aes(x=maxlength, y=He), method="glm",  method.args = list(family = "binomial"), #Add linear regression line
              formula= cbind(mtdna_maxlength_He_no.na_nowhaleshark$success, mtdna_maxlength_He_no.na_nowhaleshark$failure)~x, inherit.aes=FALSE, color= "blue") +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. He", subtitle = "(w/out Rhincodon typus) msat vs. mtDNA") + #add plot title
  xlab("Max Length (cm)") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE) 


#Fecundity Mean

final_fecunditymean_alln <- final_fecunditymean_all[!is.na(final_fecunditymean_all$n),] #exclude NA
final_fecunditymean_alln <- final_fecunditymean_alln[!(final_fecunditymean_alln$n=="25-32"),] #exclude n=25-32

final_fecunditymean_alln$success <- round(as.numeric(final_fecunditymean_alln$He)*(as.numeric(final_fecunditymean_alln$n))) #create column of successes
final_fecunditymean_alln$failure <- round(as.numeric(1-final_fecunditymean_alln$He)*(as.numeric(final_fecunditymean_alln$n))) #create column of failures

ggplot(final_fecunditymean_alln, aes(x=fecundity_mean, y=He, col=markertype, shape=markertype)) + #fecundity mean & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(data= msat_fecundity_He_no.nan , aes(x=fecundity_mean, y=He), method="glm",  method.args = list(family = "binomial"), 
              formula= cbind(msat_fecundity_He_no.nan $success, msat_fecundity_He_no.nan $failure)~x, inherit.aes = FALSE,   #Add linear regression line outline
              se=TRUE, color = "black", size = 2, fill = NA) + 
  stat_smooth(data= msat_fecundity_He_no.nan , aes(x=fecundity_mean, y=He), method="glm",  method.args = list(family = "binomial"), #Add linear regression line
              formula= cbind(msat_fecundity_He_no.nan $success, msat_fecundity_He_no.nan $failure)~x, inherit.aes = FALSE, color= "skyblue2") +
  geom_smooth(data= mtdna_fecundity_He_no.na, aes(x=fecundity_mean, y=He), method="glm",  method.args = list(family = "binomial"), 
              formula= cbind(mtdna_fecundity_He_no.na$success, mtdna_fecundity_He_no.na$failure)~x, inherit.aes=FALSE,  #Add linear regression line outline
              se=TRUE, color = "black", size = 2, fill = NA) + 
  stat_smooth(data= mtdna_fecundity_He_no.na, aes(x=fecundity_mean, y=He), method="glm",  method.args = list(family = "binomial"), #Add linear regression line
              formula= cbind(mtdna_fecundity_He_no.na$success, mtdna_fecundity_He_no.na$failure)~x, inherit.aes=FALSE, color= "blue") +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Fecundity Mean vs. He", subtitle = "msat vs. mtDNA") + #add plot title
  xlab("Fecundity Mean") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE) 
