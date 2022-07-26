################################################### Script for US Data Graphs vs. Full Data Graphs  ########################################################

#graph US mtDNA & msat datasets 

##########################################################################################################################################

######### Set-up #########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(dplyr)

#read in data
mtdna_data_new <- read.csv("new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("new_msat_full_US_data.csv", stringsAsFactors = FALSE) #read in

mtdna_FULL <- read.csv("new_full_mtdna_LH.csv", stringsAsFactors = FALSE) #read in
msat_FULL <- read.csv("new_full_msat_LH.csv", stringsAsFactors = FALSE) #read in

### Needed Equation ###
lm_eqn = function(x, y, df){ #set up formula for regression line equation
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == b %.% italic(x) + a,
                   list(a = format(coef(m)[[1]], digits = 2), 
                        b = format(coef(m)[[2]], digits = 2)))
  as.character(as.expression(eq));                 
}

### Center Title and Subtitle for Graphs ###
theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

###Combined old full data set to new full data set w/ LH traits: No need to go through code again###

#msat
old <- read.csv("msatloci.csv")
new <- read.csv("msat_spp_final.csv")
new_full_msat_LH <- merge(x = old, y = new[ , c("spp","maxlength","fecundity_mean","fertilization","reproductionmode")], by = "spp", all.x=TRUE)
write.csv(new_full_msat_LH,'new_full_msat_LH.csv')

#mtdna
old_mt <- read.csv("mtdna_assembled_new.csv")
new_mt <- read.csv("mtdna_spp_final.csv")
new_full_mtdna_LH <- merge(x = old_mt, y = new_mt[ , c("spp","maxlength","fecundity_mean","fertilization","reproductionmode")], by = "spp", all.x=TRUE)
write.csv(new_full_mtdna_LH,'new_full_mtdna_LH.csv')

############################################ US DATA ############################################ 

#mtDNA: He#

#Final mtDNA fertilization
mtdna_data_new$final_fertilization <- NA #create new column to categorize fertilization methods

mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="external"]  <- "external"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

mtdna_final_fertilization_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_fertilization) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

#Final mtDNA Reproduction mode

mtdna_data_new$final_reproductionmode <- NA #create new column to categorize reproduction mode

mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

mtdna_final_reproductionmode_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_reproductionmode ) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

#Final mtDNA Specific Reproduction mode

mtdna_data_new$specific.repro_mode <- NA #create new column to categorize hermaphrodite type

mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="dioecism"]  <- "Dioecism"
mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="protogyny"] <- "Protogyny"
mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="protandry"] <- "Protandry"
mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

mtdna_reproduction_type_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$specific.repro_mode ) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

#Final mtDNA Max Length

mtdna_maxlength_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$maxlength) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

mtdna_maxlength_He_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlength_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlength_He_no.na$logtransform.maxlength <- log10(mtdna_maxlength_He_no.na$maxlength)
}

#Final mtDNA Fecundity Mean
mtdna_fecundity_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$fecundity_mean) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

mtdna_fecundity_He_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_fecundity_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_fecundity_He_no.na$logtransform.fecundity <- log10(mtdna_fecundity_He_no.na$fecundity_mean)
}

#mtDNA: Pi#

#Final mtDNA Fertilization 

mtdna_data_new$final_fertilization <- NA #create new column to categorize fertilization methods

mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="external"]  <- "external"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

mtdna_final_fertilization_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_fertilization) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

#Final mtDNA Reproduction mode

mtdna_data_new$final_reproductionmode <- NA #create new column to categorize reproduction mode

mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

mtdna_final_reproductionmode_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_reproductionmode ) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

#Final mtDNA Specific Reproduction mode

mtdna_data_new$specific.repro_mode <- NA #create new column to categorize hermaphrodite type

mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="dioecism"]  <- "Dioecism"
mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="protogyny"] <- "Protogyny"
mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="protandry"] <- "Protandry"
mtdna_data_new$specific.repro_mode  [mtdna_data_new$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

mtdna_reproduction_type_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$specific.repro_mode ) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

#Final mtDNA Max Length

mtdna_maxlength_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$maxlength) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

mtdna_maxlength_Pi_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlength_Pi_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlength_Pi_no.na$logtransform.maxlength <- log10(mtdna_maxlength_Pi_no.na$maxlength)
}

#Final mtDNA Fecundity Mean

mtdna_fecundity_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$fecundity_mean) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

mtdna_fecundity_Pi_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_fecundity_Pi_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_fecundity_Pi_no.na$logtransform.fecundity <- log10(mtdna_fecundity_Pi_no.na$fecundity_mean)
}

#msat#

#Final msat fertilization
msat_data$final_fertilization <- NA #create new column to categorize fertilization methods

msat_data$final_fertilization [msat_data$fertilization =="external"]  <- "external"
msat_data$final_fertilization [msat_data$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
msat_data$final_fertilization [msat_data$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

msat_final_fertilization_He_no.na <- msat_data[!is.na(msat_data$final_fertilization) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

#Final msat Reproduction mode


msat_data$final_reproductionmode <- NA #create new column to categorize reproduction mode

msat_data$final_reproductionmode  [msat_data$reproductionmode =="dioecism"]  <- "Dioecious"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

msat_final_reproductionmode_He_no.na <- msat_data[!is.na(msat_data$final_reproductionmode ) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

#Final msat Specific Reproduction mode

msat_data$specific.repro_mode <- NA #create new column to categorize reproduction type

msat_data$specific.repro_mode  [msat_data$reproductionmode =="dioecism"]  <- "Dioecism"
msat_data$specific.repro_mode  [msat_data$reproductionmode =="protogyny"] <- "Protogyny"
msat_data$specific.repro_mode  [msat_data$reproductionmode =="protandry"] <- "Protandry"
msat_data$specific.repro_mode  [msat_data$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

msat_reproduction_type_He_no.na <- msat_data[!is.na(msat_data$specific.repro_mode ) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

#Final msat Max Length

msat_maxlength_He_no.na <- msat_data[!is.na(msat_data$maxlength) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

msat_maxlength_He_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(msat_maxlength_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlength_He_no.na$logtransform.maxlength <- log10(msat_maxlength_He_no.na$maxlength)
}

#Final msat Fecundity Mean: He##
msat_fecundity_He_no.na <- msat_data[!is.na(msat_data$fecundity_mean) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

msat_fecundity_He_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(msat_fecundity_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_fecundity_He_no.na$logtransform.fecundity <- log10(msat_fecundity_He_no.na$fecundity_mean)
}

################################ Graphing and Comparing Data ################################

############## Box Plots: Character Data ##############

########### Fertilization ###########
### msat vs. mtDNA: He
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

Fertplot1 <- ggplot(final_fertilization_all) + geom_boxplot(aes(x = final_fertilization, y = He, fill=markertype)) + #final fertilization & He box plot
  ggtitle("Fertilization Method vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("#00429d", "#a5d5d8"))

### mtDNA: Pi
Fertplot2 <- ggplot(mtdna_final_fertilization_Pi_no.na) + geom_boxplot(aes(x = final_fertilization, y = Pi, fill = final_fertilization)) + #final fertilization & He box plot
  ggtitle("Fertilization Method vs. Pi", subtitle= "mtDNA") + #add plot title
  xlab("Fertilization Method") + ylab("Pi") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  labs(fill = "Fertilization") +
  scale_fill_manual(values=c("#a5d5d8", "#a5d5d8"))

#Graph msat vs. mtDNA: He & mtDNA: Pi side-by-side
Fertplot1 + Fertplot2

########### Reproduction Mode ###########
### msat vs. mtDNA: He
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

Reproplot1 <- ggplot(reproductionmode_all) + geom_boxplot(aes(x = final_reproductionmode, y = He, fill= markertype)) + #final fertilization & He box plot
  ggtitle("Reproduction Mode vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("#00429d", "#a5d5d8"))

### mtDNA: Pi
Reproplot2 <- ggplot(mtdna_final_reproductionmode_Pi_no.na) + geom_boxplot(aes(x = final_reproductionmode, y = Pi, fill = final_reproductionmode)) + #final fertilization & He box plot
  ggtitle("Reproduction Mode vs. Pi", subtitle= "mtDNA") + #add plot title
  xlab("Reproduction Mode") + ylab("Pi") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  labs(fill = "Reproduction Mode") +
  scale_fill_manual(values=c("#a5d5d8", "#a5d5d8"))

#Graph msat vs. mtDNA: He & mtDNA: Pi side-by-side
Reproplot1 + Reproplot2

#####Scatter Plots: Numerical Data#####

########### Max Length ###########
### msat vs. mtDNA: He
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

Maxplot1 <- ggplot(final_maxlength_all, aes(x=logtransform.maxlength, y=He, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL, col=as.factor(markertype))) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. He", subtitle = "msat vs. mtDNA") + #add plot title
  xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.5, y = 0.55, label = lm_eqn(msat_maxlength_He_no.na$logtransform.maxlength, msat_maxlength_He_no.na$He, msat_maxlength_He_no.na), 
           color="#00429d", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 2, y = 0.91, label = lm_eqn(mtdna_maxlength_He_no.na$logtransform.maxlength, mtdna_maxlength_He_no.na$He, mtdna_maxlength_He_no.na), 
           color="#a5d5d8", size = 5, parse=TRUE, alpha = 0.8) + #add regression line equation
  scale_colour_manual(values=c("#00429d", "#a5d5d8")) +
  scale_shape(solid = FALSE)

### mtDNA: Pi
Maxplot2 <-ggplot(mtdna_maxlength_Pi_no.na, aes(x=logtransform.maxlength, y=Pi, color=as.factor(logtransform.maxlength))) + #max length & He scatter plot
  geom_point(aes(shape=2, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ggtitle("Maxlength vs. Pi", subtitle = "mtDNA") + #add plot title
  xlab("Maxlength (log(cm))") + ylab("Pi") + #add axis labels 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.8, y = 0.02, label = lm_eqn(mtdna_maxlength_Pi_no.na$logtransform.maxlength, mtdna_maxlength_Pi_no.na$Pi, mtdna_maxlength_Pi_no.na), 
           color="#a5d5d8", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("#a5d5d8")) +
  scale_shape(solid = FALSE)

#Graph msat vs. mtDNA: He & mtDNA: Pi side-by-side
Maxplot1 + Maxplot2

########### Fecundity Mean ###########
#msat vs. mtDNA: He
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

############## Misc. Graphs ############## 

### Comparing Max length vs. Fecundity: mtDNA ###

mtdna_maxlengthandfecunditymean_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$maxlength) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest
mtdna_maxlengthandfecunditymean_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$fecundity_mean) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

mtdna_maxlengthandfecunditymean_He_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlengthandfecunditymean_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_He_no.na$logtransform.maxlength <- log10(mtdna_maxlengthandfecunditymean_He_no.na$maxlength)
}

mtdna_maxlengthandfecunditymean_He_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_maxlengthandfecunditymean_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_He_no.na$logtransform.fecundity <- log10(mtdna_maxlengthandfecunditymean_He_no.na$fecundity_mean)
}

### Comparing Max length vs. Fecundity: msat ###

msat_maxlengthandfecunditymean_He_no.na <- msat_data[!is.na(msat_data$maxlength) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest
msat_maxlengthandfecunditymean_He_no.na <- msat_data[!is.na(msat_data$fecundity_mean) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

msat_maxlengthandfecunditymean_He_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(msat_maxlengthandfecunditymean_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlengthandfecunditymean_He_no.na$logtransform.maxlength <- log10(msat_maxlengthandfecunditymean_He_no.na$maxlength)
}

msat_maxlengthandfecunditymean_He_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(msat_maxlengthandfecunditymean_He_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlengthandfecunditymean_He_no.na$logtransform.fecundity <- log10(msat_maxlengthandfecunditymean_He_no.na$fecundity_mean)
}

#Graph

final_maxlength_all_mlvsfm = merge(mtdna_maxlengthandfecunditymean_He_no.na, msat_maxlengthandfecunditymean_He_no.na, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_maxlength_all_mlvsfm$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "mtdna101"]  <- "mtDNA"
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "mtdna102"]  <- "mtDNA"
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "mtdna103"]  <- "mtDNA"
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats000"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats001"]  <- "msat"
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats002"]  <- "msat"
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats200"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats201"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats100"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats101"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats301"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats302"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats303"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats304"]  <- "msat" 
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file == "msats305"]  <- "msat"
final_maxlength_all_mlvsfm$markertype [final_maxlength_all_mlvsfm$file ==	"ppdat"]  <- "msat" 

ggplot(final_maxlength_all_mlvsfm, aes(x=logtransform.maxlength, y=logtransform.fecundity, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  #ylim(0,1)+                              #create limits
  #coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. Fecundity Mean", subtitle = "msat vs. mtDNA") + #add plot title
  xlab("Max Length (log(cm))") + ylab("Fecundity Mean (log)") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 2.2, y = 4, label = lm_eqn(msat_maxlengthandfecunditymean_He_no.na$logtransform.maxlength, msat_maxlengthandfecunditymean_He_no.na$logtransform.fecundity, msat_maxlengthandfecunditymean_He_no.na), 
           color="skyblue3", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 2, y = 7, label = lm_eqn(mtdna_maxlengthandfecunditymean_He_no.na$logtransform.maxlength, mtdna_maxlengthandfecunditymean_He_no.na$logtransform.fecundity, mtdna_maxlengthandfecunditymean_He_no.na), 
           color="blue", size = 5, parse=TRUE, alpha = 0.8) + #add regression line equation
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)

### Comparing He & Pi: mtDNA ###

#add marker type based on file type
mtdna_data_new$markertype <- NA #create new column to categorize marker type

mtdna_data_new$markertype [mtdna_data_new$file == "mtdna101"]  <- "mtDNA"
mtdna_data_new$markertype [mtdna_data_new$file == "mtdna102"]  <- "mtDNA"
mtdna_data_new$markertype [mtdna_data_new$file == "mtdna103"]  <- "mtDNA"
mtdna_data_new$markertype [mtdna_data_new$file == "msats000"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats001"]  <- "msat"
mtdna_data_new$markertype [mtdna_data_new$file == "msats002"]  <- "msat"
mtdna_data_new$markertype [mtdna_data_new$file == "msats200"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file =="msats201"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats100"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats101"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats301"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats302"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats303"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats304"]  <- "msat" 
mtdna_data_new$markertype [mtdna_data_new$file == "msats305"]  <- "msat"
mtdna_data_new$markertype [mtdna_data_new$file =="ppdat"]  <- "msat" 

ggplot(mtdna_data_new, aes(x=He, y=Pi, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  #ylim(0,1)+                              #create limits
  #coord_cartesian(ylim = c(0, 1)) +
  ggtitle("He vs. Pi", subtitle = "mtDNA") + #add plot title
  xlab("He") + ylab("Pi") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 0.5, y = 0.02, label = lm_eqn(mtdna_data_new$He, mtdna_data_new$Pi, mtdna_data_new), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("blue")) +
  scale_shape(solid = FALSE)

### Pi Graphs ###

## Fertilization ##


## Reproduction Mode ##


## Specific Repro Mode##
ggplot(mtdna_reproduction_type_Pi_no.na) + geom_boxplot(aes(x = specific.repro_mode, y = Pi, fill = specific.repro_mode)) + #final fertilization & He box plot
  ggtitle("Specific Reproduction Mode vs. Pi", subtitle= "mtDNA") + #add plot title
  xlab("Specific Repro Mode") + ylab("Pi") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  labs(fill = "Specific Repro Mode") +
  scale_fill_manual(values=c("turquoise2", "darkcyan"))

## Maxlength ##
## Fecundity Mean ##
ggplot(mtdna_fecundity_Pi_no.na, aes(x=logtransform.fecundity, y=Pi)) + #max length & He scatter plot
  geom_point(aes(fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ggtitle("Fecundity Mean vs. Pi", subtitle = "mtDNA") + #add plot title
  xlab("Fecundity Mean") + ylab("Pi") + #add axis labels 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 3.5, y = 0.02, label = lm_eqn(mtdna_fecundity_Pi_no.na$logtransform.fecundity, mtdna_fecundity_Pi_no.na$Pi, mtdna_fecundity_Pi_no.na), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("blue")) +
  scale_shape(solid = FALSE)

### Comparing Max length vs. Fecundity for Pi: mtDNA ###

mtdna_maxlengthandfecunditymean_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$maxlength) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest
mtdna_maxlengthandfecunditymean_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$fecundity_mean) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

mtdna_maxlengthandfecunditymean_Pi_no.na$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlengthandfecunditymean_Pi_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_Pi_no.na$logtransform.maxlength <- log10(mtdna_maxlengthandfecunditymean_Pi_no.na$maxlength)
}

mtdna_maxlengthandfecunditymean_Pi_no.na$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_maxlengthandfecunditymean_Pi_no.na)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_Pi_no.na$logtransform.fecundity <- log10(mtdna_maxlengthandfecunditymean_Pi_no.na$fecundity_mean)
}

#Graph Maxlength vs. Fecundity for Pi
ggplot(mtdna_maxlengthandfecunditymean_Pi_no.na, aes(x=logtransform.maxlength, y=logtransform.fecundity)) + #max length & He scatter plot
  geom_point(aes(fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ggtitle("Maxlength vs. Fecundity Mean", subtitle = "mtDNA: US (log transformed data)") + #add plot title
  xlab("Maxlength") + ylab("Fecundity Mean") + #add axis labels 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 2, y = 2, label = lm_eqn(mtdna_maxlengthandfecunditymean_Pi_no.na$logtransform.maxlength, mtdna_maxlengthandfecunditymean_Pi_no.na$logtransform.fecundity, mtdna_maxlengthandfecunditymean_Pi_no.na), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("blue")) +
  scale_shape(solid = FALSE)

############################################ Full DATA ############################################ 
#mtDNA: He#

#Final mtDNA fertilization
mtdna_FULL$final_fertilization <- NA #create new column to categorize fertilization methods

mtdna_FULL$final_fertilization [mtdna_FULL$fertilization =="external"]  <- "external"
mtdna_FULL$final_fertilization [mtdna_FULL$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
mtdna_FULL$final_fertilization [mtdna_FULL$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

mtdna_final_fertilization_He_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$final_fertilization) & !is.na(mtdna_FULL$He),] #create new table that excludes NA's from columns of interest

#Final mtDNA Reproduction mode

mtdna_FULL$final_reproductionmode <- NA #create new column to categorize reproduction mode

mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

mtdna_final_reproductionmode_He_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$final_reproductionmode ) & !is.na(mtdna_FULL$He),] #create new table that excludes NA's from columns of interest

#Final mtDNA Specific Reproduction mode

mtdna_FULL$specific.repro_mode <- NA #create new column to categorize hermaphrodite type

mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="dioecism"]  <- "Dioecism"
mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="protogyny"] <- "Protogyny"
mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="protandry"] <- "Protandry"
mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

mtdna_reproduction_type_He_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$specific.repro_mode ) & !is.na(mtdna_FULL$He),] #create new table that excludes NA's from columns of interest

#Final mtDNA Max Length

mtdna_maxlength_He_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$maxlength) & !is.na(mtdna_FULL$He),] #create new table that excludes NA's from columns of interest

mtdna_maxlength_He_no.na_FULL$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlength_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlength_He_no.na_FULL$logtransform.maxlength <- log10(mtdna_maxlength_He_no.na_FULL$maxlength)
}

#Final mtDNA Fecundity Mean
mtdna_fecundity_He_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$fecundity_mean) & !is.na(mtdna_FULL$He),] #create new table that excludes NA's from columns of interest

mtdna_fecundity_He_no.na_FULL$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_fecundity_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_fecundity_He_no.na_FULL$logtransform.fecundity <- log10(mtdna_fecundity_He_no.na_FULL$fecundity_mean)
}

#mtDNA: Pi#

#Final mtDNA Fertilization 

mtdna_FULL$final_fertilization <- NA #create new column to categorize fertilization methods

mtdna_FULL$final_fertilization [mtdna_FULL$fertilization =="external"]  <- "external"
mtdna_FULL$final_fertilization [mtdna_FULL$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
mtdna_FULL$final_fertilization [mtdna_FULL$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

mtdna_final_fertilization_Pi_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$final_fertilization) & !is.na(mtdna_FULL$Pi),] #create new table that excludes NA's from columns of interest

#Final mtDNA Reproduction mode

mtdna_FULL$final_reproductionmode <- NA #create new column to categorize reproduction mode

mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
mtdna_FULL$final_reproductionmode  [mtdna_FULL$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

mtdna_final_reproductionmode_Pi_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$final_reproductionmode ) & !is.na(mtdna_FULL$Pi),] #create new table that excludes NA's from columns of interest

#Final mtDNA Specific Reproduction mode

mtdna_FULL$specific.repro_mode <- NA #create new column to categorize hermaphrodite type

mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="dioecism"]  <- "Dioecism"
mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="protogyny"] <- "Protogyny"
mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="protandry"] <- "Protandry"
mtdna_FULL$specific.repro_mode  [mtdna_FULL$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

mtdna_reproduction_type_Pi_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$specific.repro_mode ) & !is.na(mtdna_FULL$Pi),] #create new table that excludes NA's from columns of interest

#Final mtDNA Max Length

mtdna_maxlength_Pi_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$maxlength) & !is.na(mtdna_FULL$Pi),] #create new table that excludes NA's from columns of interest

mtdna_maxlength_Pi_no.na_FULL$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlength_Pi_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlength_Pi_no.na_FULL$logtransform.maxlength <- log10(mtdna_maxlength_Pi_no.na_FULL$maxlength)
}

#Final mtDNA Fecundity Mean

mtdna_fecundity_Pi_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$fecundity_mean) & !is.na(mtdna_FULL$Pi),] #create new table that excludes NA's from columns of interest

mtdna_fecundity_Pi_no.na_FULL$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_fecundity_Pi_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_fecundity_Pi_no.na_FULL$logtransform.fecundity <- log10(mtdna_fecundity_Pi_no.na_FULL$fecundity_mean)
}

#msat#

#Final msat fertilization
msat_FULL$final_fertilization <- NA #create new column to categorize fertilization methods

msat_FULL$final_fertilization [msat_FULL$fertilization =="external"]  <- "external"
msat_FULL$final_fertilization [msat_FULL$fertilization =="internal (oviduct)"] <- "internal (oviduct)"
msat_FULL$final_fertilization [msat_FULL$fertilization =="in brood pouch or similar structure"] <- "internal (oviduct)" #convert "in brood pouch or similar structure" to internal fertilization

msat_final_fertilization_He_no.na_FULL <- msat_FULL[!is.na(msat_FULL$final_fertilization) & !is.na(msat_FULL$He),] #create new table that excludes NA's from columns of interest

#Final msat Reproduction mode


msat_FULL$final_reproductionmode <- NA #create new column to categorize reproduction mode

msat_FULL$final_reproductionmode  [msat_FULL$reproductionmode =="dioecism"]  <- "Dioecious"
msat_FULL$final_reproductionmode  [msat_FULL$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
msat_FULL$final_reproductionmode  [msat_FULL$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
msat_FULL$final_reproductionmode  [msat_FULL$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

msat_final_reproductionmode_He_no.na_FULL <- msat_FULL[!is.na(msat_FULL$final_reproductionmode ) & !is.na(msat_FULL$He),] #create new table that excludes NA's from columns of interest

#Final msat Specific Reproduction mode

msat_FULL$specific.repro_mode <- NA #create new column to categorize reproduction type

msat_FULL$specific.repro_mode  [msat_FULL$reproductionmode =="dioecism"]  <- "Dioecism"
msat_FULL$specific.repro_mode  [msat_FULL$reproductionmode =="protogyny"] <- "Protogyny"
msat_FULL$specific.repro_mode  [msat_FULL$reproductionmode =="protandry"] <- "Protandry"
msat_FULL$specific.repro_mode  [msat_FULL$reproductionmode =="true hermaphroditism"] <- "True Hermaphroditism"

msat_reproduction_type_He_no.na_FULL <- msat_FULL[!is.na(msat_FULL$specific.repro_mode ) & !is.na(msat_FULL$He),] #create new table that excludes NA's from columns of interest

#Final msat Max Length

msat_maxlength_He_no.na_FULL <- msat_FULL[!is.na(msat_FULL$maxlength) & !is.na(msat_FULL$He),] #create new table that excludes NA's from columns of interest

msat_maxlength_He_no.na_FULL$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(msat_maxlength_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlength_He_no.na_FULL$logtransform.maxlength <- log10(msat_maxlength_He_no.na_FULL$maxlength)
}

#Final msat Fecundity Mean: He##
msat_fecundity_He_no.na_FULL <- msat_FULL[!is.na(msat_FULL$fecundity_mean) & !is.na(msat_FULL$He),] #create new table that excludes NA's from columns of interest

msat_fecundity_He_no.na_FULL$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(msat_fecundity_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_fecundity_He_no.na_FULL$logtransform.fecundity <- log10(msat_fecundity_He_no.na_FULL$fecundity_mean)
}

#####Box Plots: Character Data#####

#Fertilization#

final_fertilization_all_FULL = merge(mtdna_final_fertilization_He_no.na_FULL, msat_final_fertilization_He_no.na_FULL, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_fertilization_all_FULL$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "mtdna101"]  <- "mtDNA"
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "mtdna102"]  <- "mtDNA"
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "mtdna103"]  <- "mtDNA"
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats000"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats001"]  <- "msat"
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats002"]  <- "msat"
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats200"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats201"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats100"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats101"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats301"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats302"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats303"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats304"]  <- "msat" 
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file == "msats305"]  <- "msat"
final_fertilization_all_FULL$markertype [final_fertilization_all_FULL$file ==	"ppdat"]  <- "msat" 

ggplot(final_fertilization_all_FULL) + geom_boxplot(aes(x = final_fertilization, y = He, fill=markertype)) + #final fertilization & He box plot
  ggtitle("Fertilization Method vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("turquoise2", "darkcyan"))

#Reproduction Mode#
reproductionmode_all_FULL = merge(msat_final_reproductionmode_He_no.na_FULL, mtdna_final_reproductionmode_He_no.na_FULL, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

reproductionmode_all_FULL$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "mtdna101"]  <- "mtDNA"
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "mtdna102"]  <- "mtDNA"
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "mtdna103"]  <- "mtDNA"
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats000"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats001"]  <- "msat"
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats002"]  <- "msat"
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats200"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats201"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats100"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats101"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats301"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats302"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats303"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats304"]  <- "msat" 
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file == "msats305"]  <- "msat"
reproductionmode_all_FULL$markertype [reproductionmode_all_FULL$file ==	"ppdat"]  <- "msat" 

ggplot(reproductionmode_all_FULL) + geom_boxplot(aes(x = final_reproductionmode, y = He, fill= markertype)) + #final fertilization & He box plot
  ggtitle("Reproduction Mode vs. He", subtitle= "msat vs. mtDNA") + #add plot title
  xlab("Reproduction Mode") + ylab("He") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  scale_fill_manual(values=c("turquoise2", "darkcyan"))


#Specific Reproduction Mode#
specificreproductionmode_all_FULL = merge(msat_reproduction_type_He_no.na_FULL, mtdna_reproduction_type_He_no.na_FULL, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

specificreproductionmode_all_FULL$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "mtdna101"]  <- "mtDNA"
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "mtdna102"]  <- "mtDNA"
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "mtdna103"]  <- "mtDNA"
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats000"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats001"]  <- "msat"
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats002"]  <- "msat"
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats200"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats201"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats100"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats101"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats301"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats302"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats303"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats304"]  <- "msat" 
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file == "msats305"]  <- "msat"
specificreproductionmode_all_FULL$markertype [specificreproductionmode_all_FULL$file ==	"ppdat"]  <- "msat" 

ggplot(specificreproductionmode_all_FULL) + geom_boxplot(aes(x = specific.repro_mode, y = He, fill=markertype)) + #final fertilization & He box plot
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

#w/ outliers
final_maxlength_all_FULL = merge(mtdna_maxlength_He_no.na_FULL, msat_maxlength_He_no.na_FULL, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_maxlength_all_FULL$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "mtdna101"]  <- "mtDNA"
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "mtdna102"]  <- "mtDNA"
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "mtdna103"]  <- "mtDNA"
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats000"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats001"]  <- "msat"
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats002"]  <- "msat"
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats200"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats201"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats100"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats101"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats301"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats302"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats303"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats304"]  <- "msat" 
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file == "msats305"]  <- "msat"
final_maxlength_all_FULL$markertype [final_maxlength_all_FULL$file ==	"ppdat"]  <- "msat" 

ggplot(final_maxlength_all_FULL, aes(x=logtransform.maxlength, y=He, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ylim(0,1)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. He", subtitle = "msat vs. mtDNA") + #add plot title
  xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.5, y = 0.55, label = lm_eqn(msat_maxlength_He_no.na_FULL$logtransform.maxlength, msat_maxlength_He_no.na_FULL$He, msat_maxlength_He_no.na_FULL), 
           color="skyblue3", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 2, y = 0.91, label = lm_eqn(mtdna_maxlength_He_no.na_FULL$logtransform.maxlength, mtdna_maxlength_He_no.na_FULL$He, mtdna_maxlength_He_no.na_FULL), 
           color="blue", size = 5, parse=TRUE, alpha = 0.8) + #add regression line equation
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)

#w/out outliers
final_maxlength_all_FULLo = merge(mtdna_maxlength_He_no.na, msat_maxlength_He_no.na, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_maxlength_all_FULLo$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "mtdna101"]  <- "mtDNA"
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "mtdna102"]  <- "mtDNA"
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "mtdna103"]  <- "mtDNA"
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats000"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats001"]  <- "msat"
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats002"]  <- "msat"
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats200"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats201"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats100"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats101"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats301"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats302"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file =="msats303"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats304"]  <- "msat" 
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file == "msats305"]  <- "msat"
final_maxlength_all_FULLo$markertype [final_maxlength_all_FULLo$file =="ppdat"]  <- "msat" 

ggplot(final_maxlength_all_FULLo, aes(x=logtransform.maxlength, y=He, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  xlim(0.95,2.65)+                              #create limits
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. He (w/out smallest & largest spp)", subtitle = "msat vs. mtDNA") + #add plot title
  xlab("Max Length (log(cm))") + ylab("He") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 2, y = 0.60, label = lm_eqn(msat_maxlength_He_no.na$logtransform.maxlength, msat_maxlength_He_no.na$He, msat_maxlength_He_no.na), 
           color="skyblue3", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 2, y = 0.85, label = lm_eqn(mtdna_maxlength_He_no.na$logtransform.maxlength, mtdna_maxlength_He_no.na$He, mtdna_maxlength_He_no.na), 
           color="blue", size = 5, parse=TRUE, alpha = 0.8) + #add regression line equation
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)

#Fecundity Mean#

final_fecunditymean_all_FULL = merge(msat_fecundity_He_no.na_FULL, mtdna_fecundity_He_no.na_FULL, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_fecunditymean_all_FULL$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "mtdna101"]  <- "mtDNA"
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "mtdna102"]  <- "mtDNA"
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "mtdna103"]  <- "mtDNA"
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats000"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats001"]  <- "msat"
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats002"]  <- "msat"
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats200"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats201"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats100"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats101"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats301"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats302"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats303"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats304"]  <- "msat" 
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file == "msats305"]  <- "msat"
final_fecunditymean_all_FULL$markertype [final_fecunditymean_all_FULL$file ==	"ppdat"]  <- "msat" 

ggplot(final_fecunditymean_all_FULL, aes(x=logtransform.fecundity, y=He, col=markertype, shape=markertype)) + #fecundity mean & He scatter plot
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
  annotate(geom="label", x = 2.8, y = 0.9, label = lm_eqn(msat_fecundity_He_no.na_FULL$logtransform.fecundity, msat_fecundity_He_no.na_FULL$He, msat_fecundity_He_no.na_FULL), 
           color="skyblue2", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 3.4, y = 0.5, label = lm_eqn(mtdna_fecundity_He_no.na_FULL$logtransform.fecundity, mtdna_fecundity_He_no.na_FULL$He, mtdna_fecundity_He_no.na_FULL), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_color_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)

############## Misc. Graphs ############## 

### Comparing Max length vs. Fecundity: mtDNA ###

mtdna_maxlengthandfecunditymean_He_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$maxlength) & !is.na(mtdna_FULL$He),] #create new table that excludes NA's from columns of interest
mtdna_maxlengthandfecunditymean_He_no.na_FULL <- mtdna_FULL[!is.na(mtdna_FULL$fecundity_mean) & !is.na(mtdna_FULL$He),] #create new table that excludes NA's from columns of interest

mtdna_maxlengthandfecunditymean_He_no.na_FULL$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(mtdna_maxlengthandfecunditymean_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_He_no.na_FULL$logtransform.maxlength <- log10(mtdna_maxlengthandfecunditymean_He_no.na_FULL$maxlength)
}

mtdna_maxlengthandfecunditymean_He_no.na_FULL$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_maxlengthandfecunditymean_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_He_no.na_FULL$logtransform.fecundity <- log10(mtdna_maxlengthandfecunditymean_He_no.na_FULL$fecundity_mean)
}

### Comparing Max length vs. Fecundity: msat ###

msat_maxlengthandfecunditymean_He_no.na_FULL <- msat_data[!is.na(msat_data$maxlength) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest
msat_maxlengthandfecunditymean_He_no.na_FULL <- msat_data[!is.na(msat_data$fecundity_mean) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

msat_maxlengthandfecunditymean_He_no.na_FULL$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1:nrow(msat_maxlengthandfecunditymean_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlengthandfecunditymean_He_no.na_FULL$logtransform.maxlength <- log10(msat_maxlengthandfecunditymean_He_no.na_FULL$maxlength)
}

v$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(msat_maxlengthandfecunditymean_He_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  msat_maxlengthandfecunditymean_He_no.na_FULL$logtransform.fecundity <- log10(msat_maxlengthandfecunditymean_He_no.na_FULL$fecundity_mean)
}

#Graph

final_maxlength_all_mlvsfmFULL = merge(mtdna_maxlengthandfecunditymean_He_no.na_FULL, msat_maxlengthandfecunditymean_He_no.na_FULL, all=TRUE, no.dups= TRUE, all.x=TRUE, all.y=TRUE) #merge final fertilization data form mtdna and msat together

final_maxlength_all_mlvsfmFULL$markertype <- NA #create new column to categorize marker type

#add marker type based on file type
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "mtdna101"]  <- "mtDNA"
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "mtdna102"]  <- "mtDNA"
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "mtdna103"]  <- "mtDNA"
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats000"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats001"]  <- "msat"
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats002"]  <- "msat"
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats200"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats201"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats100"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats101"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats301"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats302"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats303"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats304"]  <- "msat" 
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file == "msats305"]  <- "msat"
final_maxlength_all_mlvsfmFULL$markertype [final_maxlength_all_mlvsfmFULL$file ==	"ppdat"]  <- "msat" 

ggplot(final_maxlength_all_mlvsfmFULL, aes(x=logtransform.maxlength, y=logtransform.fecundity, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  #ylim(0,1)+                              #create limits
  #coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Max Length vs. Fecundity Mean", subtitle = "msat vs. mtDNA") + #add plot title
  xlab("Max Length (log(cm))") + ylab("Fecundity Mean (log)") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 2.4, y = 6, label = lm_eqn(msat_maxlengthandfecunditymean_He_no.na_FULL$logtransform.maxlength, msat_maxlengthandfecunditymean_He_no.na_FULL$logtransform.fecundity, msat_maxlengthandfecunditymean_He_no.na_FULL), 
           color="skyblue3", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  annotate(geom="label", x = 2, y = 3, label = lm_eqn(mtdna_maxlengthandfecunditymean_He_no.na_FULL$logtransform.maxlength, mtdna_maxlengthandfecunditymean_He_no.na_FULL$logtransform.fecundity, mtdna_maxlengthandfecunditymean_He_no.na_FULL), 
           color="blue", size = 5, parse=TRUE, alpha = 0.8) + #add regression line equation
  scale_colour_manual(values=c("skyblue2","blue")) +
  scale_shape(solid = FALSE)

### Comparing He & Pi: mtDNA ###

#add marker type based on file type
mtdna_FULL$markertype <- NA #create new column to categorize marker type

mtdna_FULL$markertype [mtdna_FULL$file == "mtdna101"]  <- "mtDNA"
mtdna_FULL$markertype [mtdna_FULL$file == "mtdna102"]  <- "mtDNA"
mtdna_FULL$markertype [mtdna_FULL$file == "mtdna103"]  <- "mtDNA"
mtdna_FULL$markertype [mtdna_FULL$file == "msats000"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats001"]  <- "msat"
mtdna_FULL$markertype [mtdna_FULL$file == "msats002"]  <- "msat"
mtdna_FULL$markertype [mtdna_FULL$file == "msats200"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file =="msats201"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats100"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats101"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats301"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats302"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats303"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats304"]  <- "msat" 
mtdna_FULL$markertype [mtdna_FULL$file == "msats305"]  <- "msat"
mtdna_FULL$markertype [mtdna_FULL$file =="ppdat"]  <- "msat" 

ggplot(mtdna_FULL, aes(x=He, y=Pi, col=markertype, shape=markertype)) + #max length & He scatter plot
  geom_point(aes(shape=markertype, fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  #ylim(0,1)+                              #create limits
  #coord_cartesian(ylim = c(0, 1)) +
  ggtitle("He vs. Pi", subtitle = "mtDNA") + #add plot title
  xlab("He") + ylab("Pi") + #add axis labels
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 0.5, y = 0.02, label = lm_eqn(mtdna_FULL$He, mtdna_FULL$Pi, mtdna_FULL), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("blue")) +
  scale_shape(solid = FALSE)

### Pi Graphs ###

## Fertilization ##
ggplot(mtdna_final_fertilization_Pi_no.na_FULL) + geom_boxplot(aes(x = final_fertilization, y = Pi, fill = final_fertilization)) + #final fertilization & He box plot
  ggtitle("Fertilization Method vs. Pi", subtitle= "mtDNA: Full Data set") + #add plot title
  xlab("Fertilization Method") + ylab("Pi") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  labs(fill = "Fertilization") +
  scale_fill_manual(values=c("turquoise2", "darkcyan"))

## Reproduction Mode ##
ggplot(mtdna_final_reproductionmode_Pi_no.na_FULL) + geom_boxplot(aes(x = final_reproductionmode, y = Pi, fill = final_reproductionmode)) + #final fertilization & He box plot
  ggtitle("Reproduction Mode vs. Pi", subtitle= "mtDNA: Full Data set") + #add plot title
  xlab("Reproduction Mode") + ylab("Pi") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold"), 
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  labs(fill = "Reproduction Mode") +
  scale_fill_manual(values=c("turquoise2", "darkcyan"))

## Specific Repro Mode##
ggplot(mtdna_reproduction_type_Pi_no.na_FULL) + geom_boxplot(aes(x = specific.repro_mode, y = Pi, fill = specific.repro_mode)) + #final fertilization & He box plot
  ggtitle("Specific Reproduction Mode vs. Pi", subtitle= "mtDNA: Full Data set") + #add plot title
  xlab("Specific Repro Mode") + ylab("Pi") + #add axis labels
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=14, face="bold", hjust = 0.3), 
    plot.subtitle = element_text(hjust = 0.3),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  labs(fill = "Specific Repro Mode") +
  scale_fill_manual(values=c("turquoise2", "darkcyan","deepskyblue1"))

theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

## Maxlength ##
ggplot(mtdna_maxlength_Pi_no.na_FULL, aes(x=logtransform.maxlength, y=Pi)) + #max length & He scatter plot
  geom_point(aes(fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ggtitle("Maxlength vs. Pi", subtitle = "mtDNA: Full Data set") + #add plot title
  xlab("Maxlength (log(cm))") + ylab("Pi") + #add axis labels 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 1.8, y = 0.02, label = lm_eqn(mtdna_maxlength_Pi_no.na_FULL$logtransform.maxlength, mtdna_maxlength_Pi_no.na_FULL$Pi, mtdna_maxlength_Pi_no.na_FULL), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("blue")) +
  scale_shape(solid = FALSE)

## Fecundity Mean ##
ggplot(mtdna_fecundity_Pi_no.na_FULL, aes(x=logtransform.fecundity, y=Pi)) + #max length & He scatter plot
  geom_point(aes(fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ggtitle("Fecundity Mean vs. Pi", subtitle = "mtDNA: Full Data set") + #add plot title
  xlab("Fecundity Mean") + ylab("Pi") + #add axis labels 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 3.5, y = 0.02, label = lm_eqn(mtdna_fecundity_Pi_no.na_FULL$logtransform.fecundity, mtdna_fecundity_Pi_no.na_FULL$Pi, mtdna_fecundity_Pi_no.na_FULL), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("blue")) +
  scale_shape(solid = FALSE)

### Comparing Max length vs. Fecundity for Pi: mtDNA ###

mtdna_maxlengthandfecunditymean_Pi_no.na_FULL <- mtdna_data_new[!is.na(mtdna_data_new$maxlength) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest
mtdna_maxlengthandfecunditymean_Pi_no.na_FULL <- mtdna_data_new[!is.na(mtdna_data_new$fecundity_mean) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$logtransform.maxlength <- NA #add column to do a log transformation for max length

for (i in 1::nrow(mtdna_maxlengthandfecunditymean_Pi_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$logtransform.maxlength <- log10(mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$maxlength)
}

mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$logtransform.fecundity <- NA #add column to do a log transformation for fecundity mean

for (i in 1:nrow(mtdna_maxlengthandfecunditymean_Pi_no.na_FULL)) { #get log transformation data
  cat(paste(i, " ", sep = ''))
  mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$logtransform.fecundity <- log10(mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$fecundity_mean)
}

#Graph Maxlength vs. Fecundity for Pi
ggplot(mtdna_maxlengthandfecunditymean_Pi_no.na_FULL, aes(x=logtransform.maxlength, y=logtransform.fecundity)) + #max length & He scatter plot
  geom_point(aes(fill=NULL)) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 1.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, size = 1.1, fill = NA) +
  ggtitle("Maxlength vs. Fecundity Mean", subtitle = "mtDNA: Full (log transformed data)") + #add plot title
  xlab("Maxlength") + ylab("Fecundity Mean") + #add axis labels 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(color="blue", size=14, face="bold"),
    axis.title.y = element_text(color="red", size=14, face="bold"))+
  annotate(geom="label", x = 2, y = 2, label = lm_eqn(mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$logtransform.maxlength, mtdna_maxlengthandfecunditymean_Pi_no.na_FULL$logtransform.fecundity, mtdna_maxlengthandfecunditymean_Pi_no.na_FULL), 
           color="blue", size = 5, parse=TRUE, alpha=0.8) + #add regression line equation
  scale_colour_manual(values=c("blue")) +
  scale_shape(solid = FALSE)

########################################### T-Tests US Combined Marker Character Data ########################################### 

##### msat #####

##He##
#Fertilization#

external.msat <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="external"] #create vector for one aspect of t-test
internal.msat <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="internal (oviduct)"] #create vector

fertilization_ttest.msat <- t.test(external.msat, internal.msat, var.equal=FALSE) #combine created vectors & perform t-test
fertilization_ttest.msatex <- t.test(external.msat, var.equal=FALSE)
fertilization_ttest.msatin <- t.test(internal.msat, var.equal=FALSE)

#Reproduction Mode#

dioecism.msat <- msat_final_reproductionmode_He_no.na$He[msat_final_reproductionmode_He_no.na$final_reproductionmode=="Dioecious"] #create vector for one aspect of t-test
hermaphrodite.msat <- msat_final_reproductionmode_He_no.na$He[msat_final_reproductionmode_He_no.na$final_reproductionmode=="Hermaphrodite"] #create vector

reproductionmode_ttest.msat <- t.test(dioecism.msat, hermaphrodite.msat, var.equal=FALSE) #combine created vectors & perform t-test
reproductionmode_ttest.msatdio <- t.test(dioecism.msat, var.equal=FALSE)
reproductionmode_ttest.msatherm <- t.test(hermaphrodite.msat, var.equal=FALSE)

##### mtDNA #####

##He##
#Fertilization#

external.mtdnaHe <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="external"] #create vector for one aspect of t-test
internal.mtdnaHe <- msat_final_fertilization_He_no.na$He[msat_final_fertilization_He_no.na$final_fertilization=="internal (oviduct)"] #create vector

fertilization_ttest.mtdnaHe <- t.test(external.mtdnaHe, internal.mtdnaHe, var.equal=FALSE) #combine created vectors & perform t-test
fertilization_ttest.mtdnaHeex <- t.test(external.mtdnaHe, var.equal=FALSE)
fertilization_ttest.mtdnaHein <- t.test(internal.mtdnaHe, var.equal=FALSE)

#Reproduction Mode#

dioecism.mtdnaHe <- mtdna_final_reproductionmode_He_no.na$He[mtdna_final_reproductionmode_He_no.na$final_reproductionmode=="Dioecious"] #create vector for one aspect of t-test
hermaphrodite.mtdnaHe <- mtdna_final_reproductionmode_He_no.na$He[mtdna_final_reproductionmode_He_no.na$final_reproductionmode=="Hermaphrodite"] #create vector

reproductionmode_ttest.mtdnaHe <- t.test(dioecism.mtdnaHe, hermaphrodite.mtdnaHe, var.equal=FALSE) #combine created vectors & perform t-test
reproductionmode_ttest.mtdnaHedio <- t.test(dioecism.mtdnaHe, var.equal=FALSE)
reproductionmode_ttest.mtdnaHeherm <- t.test(hermaphrodite.mtdnaHe, var.equal=FALSE)

##Pi##
#Fertilization#

external.mtdnaPi <- mtdna_final_fertilization_Pi_no.na$Pi[mtdna_final_fertilization_Pi_no.na$final_fertilization=="external"] #create vector for one aspect of t-test
internal.mtdnaPi <- mtdna_final_fertilization_Pi_no.na$Pi[mtdna_final_fertilization_Pi_no.na$final_fertilization=="internal (oviduct)"] #create vector

fertilization_ttest.mtdnaPi <- t.test(external.mtdnaPi, internal.mtdnaPi, var.equal=FALSE) #combine created vectors & perform t-test
fertilization_ttest.mtdnaPiex <- t.test(external.mtdnaPi, var.equal=FALSE) 
fertilization_ttest.mtdnaPiin <- t.test(internal.mtdnaPi, var.equal=FALSE) 

#Reproduction Mode#

dioecism.mtdnaPi <- mtdna_final_reproductionmode_Pi_no.na$Pi[mtdna_final_reproductionmode_Pi_no.na$final_reproductionmode=="Dioecious"] #create vector for one aspect of t-test
hermaphrodite.mtdnaPi <- mtdna_final_reproductionmode_Pi_no.na$Pi[mtdna_final_reproductionmode_Pi_no.na$final_reproductionmode=="Hermaphrodite"] #create vector

reproductionmode_ttest.mtdnaPi <- t.test(dioecism.mtdnaPi, hermaphrodite.mtdnaPi, var.equal=FALSE) #combine created vectors & perform t-test
reproductionmode_ttest.mtdnaPidio <- t.test(dioecism.mtdnaPi, var.equal=FALSE)
reproductionmode_ttest.mtdnaPiherm <- t.test(hermaphrodite.mtdnaPi, var.equal=FALSE)

########################################### ANOVA: US Combined Marker Character Data ########################################### 


############## Fertilization Method ############## 
#msat: He
fertilizationanova.msatHe <- aov(He ~ final_fertilization, data = msat_final_fertilization_He_no.na) #perform anova test for combined data
TukeyHSD(fertilizationanova.msatHe) #perform TukeyHSD to see full table of results

#mtDNA: He
fertilizationanova.mtDNAHe <- aov(He ~ final_fertilization, data = mtdna_final_fertilization_He_no.na) #perform anova test for combined data
TukeyHSD(fertilizationanova.mtDNAHe)

#mtDNA: Pi
fertilizationanova.mtDNAPi <- aov(Pi ~ final_fertilization, data = mtdna_final_fertilization_Pi_no.na) #perform anova test for combined data
TukeyHSD(fertilizationanova.mtDNAPi)

############## Reproduction Mode ############## 
#msat: He
reproductionmodeanova.msatHe <- aov(He ~ final_reproductionmode, data = msat_final_reproductionmode_He_no.na) #perform anova test for combined data
TukeyHSD(reproductionmodeanova.msatHe) #perform TukeyHSD to see full table of results

#mtDNA: He
reproductionmodeanova.mtDNAHe <- aov(He ~ final_reproductionmode, data = mtdna_final_reproductionmode_He_no.na) #perform anova test for combined data
TukeyHSD(reproductionmodeanova.mtDNAHe) 

#mtDNA: Pi
reproductionmodeanova.mtDNAPi <- aov(Pi ~ final_reproductionmode, data = mtdna_final_reproductionmode_Pi_no.na) #perform anova test for combined data
TukeyHSD(reproductionmodeanova.mtDNAPi) 

########################################### Wilcoxon Tests: Numerical Data ########################################### 
############## Max Length ############## 

#msat: He
wilcox.test( msat_maxlength_He_no.na[ ,'maxlength'] , msat_maxlength_He_no.na[ , 'He'], paired=F) #run Wilcoxon test on max length & He

#mtDNA: He
wilcox.test( mtdna_maxlength_He_no.na[ ,'maxlength'] , mtdna_maxlength_He_no.na[ , 'He'], paired=F) #run Wilcoxon test on max length & He

#mtDNA: Pi
wilcox.test( mtdna_maxlength_Pi_no.na[ ,'maxlength'] , mtdna_maxlength_Pi_no.na[ , 'Pi'], paired=F) #run Wilcoxon test on max length & He

############## Fecundity ############## 
#msat: He
wilcox.test( msat_fecundity_He_no.na[ ,'fecundity_mean'] , msat_fecundity_He_no.na[ , 'He'], paired=F) #run Wilcoxon test on fecundity & He

#mtDNA: He
wilcox.test( mtdna_fecundity_He_no.na[ ,'fecundity_mean'] , mtdna_fecundity_He_no.na[ , 'He'], paired=F) #run Wilcoxon test on max length & He

#mtDNA: Pi
wilcox.test( mtdna_fecundity_Pi_no.na[ ,'fecundity_mean'] , mtdna_fecundity_Pi_no.na[ , 'Pi'], paired=F) #run Wilcoxon test on max length & He

########################################### Shapiro-Wilk Tests: Numerical Data ########################################### 

############## Max Length ############## 
#msat: He
shapiro.test(msat_maxlength_He_no.na$maxlength) #run Shapiro-Wilk test on max length & He

#mtDNA: He
shapiro.test(mtdna_maxlength_He_no.na$maxlength) #run Shapiro-Wilk test on max length & He

#mtDNA: Pi
shapiro.test(mtdna_maxlength_Pi_no.na$maxlength) #run Shapiro-Wilk test on max length & He

############## Fecundity ############## 
#msat: He
shapiro.test(msat_fecundity_He_no.na$fecundity_mean) #run Shapiro-Wilk test on fecundity mean & He

#mtDNA: He
shapiro.test(mtdna_fecundity_He_no.na$fecundity_mean) #run Shapiro-Wilk test on max length & He

#mtDNA: Pi
shapiro.test(mtdna_fecundity_Pi_no.na$fecundity_mean) #run Shapiro-Wilk test on max length & He

