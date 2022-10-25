################################################### Script for US Data Graphs vs. Full Data Graphs  ########################################################

#graph US mtDNA & msat datasets 

##########################################################################################################################################

######### Set-up #########

#clear environment
remove(list = ls())

#load libraries
library(tidyverse)
library(dplyr)
library(patchwork)

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
theme_update(plot.title = element_text(hjust = 0), plot.subtitle = element_text(hjust = 0.5)) #centered plot title

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
  ggtitle("(A)") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(face="bold", size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30,margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom") +
  labs(fill="Marker Type")+
  scale_fill_manual(values=c("#444444", "#999999"))

### mtDNA: Pi
Fertplot2 <- ggplot(mtdna_final_fertilization_Pi_no.na) + geom_boxplot(aes(x = final_fertilization, y = Pi, fill = final_fertilization)) + #final fertilization & He box plot
  ggtitle("(B)") + #add plot title
  xlab("Fertilization Method") + ylab("Pi") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(face="bold", size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  labs(fill = "Fertilization") +
  scale_fill_manual(values=c("#444444", "#444444"))

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
  ggtitle("(A)") + #add plot title
  xlab("Reproduction Mode") + ylab("He") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(face="bold", size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30,margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  labs(fill="Marker Type")+
  scale_fill_manual(values=c("#444444", "#999999"))

### mtDNA: Pi
Reproplot2 <- ggplot(mtdna_final_reproductionmode_Pi_no.na) + geom_boxplot(aes(x = final_reproductionmode, y = Pi, fill = final_reproductionmode)) + #final fertilization & He box plot
  ggtitle("(B)") + #add plot title
  xlab("Reproduction Mode") + ylab("Pi") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(face="bold", size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  labs(fill = "Reproduction Mode") +
  scale_fill_manual(values=c("#444444", "#444444"))

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

Maxplot1msat <- ggplot(data=subset(final_maxlength_all, markertype == 'msat' ), aes(x=logtransform.maxlength, y=He)) +
  geom_point(aes(shape=markertype, fill=NULL,color="#444444")) +    # Point aesthetics
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "#444444", size = 2, fill = NA) +
  labs(title= "(A)", x= "Max Length (log(cm))", y= "He", fill="Marker Type") + 
  theme_bw() + 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=28, face="bold"),
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  scale_fill_manual(values=c("#444444")) +
  scale_shape_manual(values = c(msat=1)) +
  scale_color_manual(values = c("#444444"), guide="none")

Maxplot1mtDNA <- ggplot(data=subset(final_maxlength_all, markertype == 'mtDNA' ), aes(x=logtransform.maxlength, y=He)) +
  geom_point(aes(shape=markertype, fill=NULL,color="#999999")) +    # Point aesthetics
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 2.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "#999999", size = 2, fill = NA) +
  labs(title= "(B)", x= "Max Length (log(cm))", y= "He", fill="Marker Type") + 
  theme_bw() + 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=28, face="bold"),
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  scale_fill_manual(values=c("#999999")) +
  scale_shape_manual(values = c(mtDNA=9)) +
  scale_color_manual(values = c("#999999"), guide="none")

### mtDNA: Pi
Maxplot2 <- ggplot(mtdna_maxlength_Pi_no.na, aes(x=logtransform.maxlength, y=Pi, fill=factor(logtransform.maxlength))) + #max length & He scatter plot
  geom_point(aes(fill=NULL), shape=9, color="#999999") +    # Point aesthetics
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 2.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "#999999", size = 2, fill = NA) +
  ggtitle("(C)") + #add plot title
  xlab("Maxlength (log(cm))") + ylab("Pi") + #add axis labels 
  theme_bw()+
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=28, face="bold"),
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom") +
  scale_fill_manual(values=c("#999999")) +
  scale_color_manual(values = c("#999999"), guide="none") +
  labs(fill="Marker Type")


#Graph msat vs. mtDNA: He & mtDNA: Pi side-by-side
Maxplot1msat + Maxplot1mtDNA + Maxplot2 

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

Fecplot1msat <- ggplot(data=subset(final_fecunditymean_all, markertype == 'msat' ), aes(x=logtransform.fecundity, y=He)) +
  geom_point(aes(shape=markertype, fill=NULL,color="#444444")) +    # Point aesthetics
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "#444444", size = 2, fill = NA) +
  labs(title= "(A)", x= "Fecundity Mean", y= "He", fill="Marker Type") + 
  theme_bw() + 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=28, face="bold"),
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  scale_fill_manual(values=c("#444444")) +
  scale_shape_manual(values = c(msat=1)) +
  scale_color_manual(values = c("#444444"), guide="none")

Fecplot1mtDNA <- ggplot(data=subset(final_fecunditymean_all, markertype == 'mtDNA' ), aes(x=logtransform.fecundity, y=He)) +
  geom_point(aes(shape=markertype, fill=NULL,color="#999999")) +    # Point aesthetics
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 2.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "#999999", size = 2, fill = NA) +
  labs(title= "(B)", x= "Fecundity Mean", y= "He", fill="Marker Type") + 
  theme_bw() + 
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=28, face="bold"),
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  scale_fill_manual(values=c("#999999")) +
  scale_shape_manual(values = c(mtDNA=9)) +
  scale_color_manual(values = c("#999999"), guide="none")

#mtDNA: Pi
Fecplot2 <- ggplot(mtdna_fecundity_Pi_no.na, aes(x=logtransform.fecundity, y=Pi, fill=factor(logtransform.fecundity))) + #max length & He scatter plot
  geom_point(aes(fill=NULL), shape=9, color="#999999") +    # Point aesthetics
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "black", size = 2.5, fill = NA) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=TRUE, color = "#999999", size = 2, fill = NA) +
  ggtitle("(C)") + #add plot title
  xlab("Fecundity Mean") + ylab("Pi") + #add axis labels 
  theme_bw()+
  theme(                                 #specifying characteristics of the plot 
    plot.title = element_text(size=28, face="bold"),
    axis.title.x = element_text(face="bold", size=30, margin = margin(t = 20)),
    axis.title.y = element_text(face="bold", size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom") +
  scale_fill_manual(values=c("#999999")) +
  scale_color_manual(values = c("#999999"), guide="none") +
  labs(fill="Marker Type")

#Graph msat vs. mtDNA: He & mtDNA: Pi side-by-side
Fecplot1msat + Fecplot1mtDNA + Fecplot2