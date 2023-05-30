################################################### Script for US Data Graphs  ########################################################

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
mtdna_data_new <- read.csv("Datasets/new_mtdna_full_US_data.csv", stringsAsFactors = FALSE) #read in 
msat_data <- read.csv("Datasets/new_msat_full_US_data.csv", stringsAsFactors = FALSE) #read in

### Needed Equation ###
lm_eqn = function(x, y, df){ #set up formula for regression line equation
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == b %.% italic(x) + a,
                   list(a = format(coef(m)[[1]], digits = 2), 
                        b = format(coef(m)[[2]], digits = 2)))
  as.character(as.expression(eq));                 
}

### Move Subtitle for Graphs ###
theme_update(plot.subtitle = element_text(hjust = 0.5)) #move plot subtitle

############################################ US DATA ############################################ 

##### mtDNA: He #####

#Final mtDNA fertilization
mtdna_data_new$final_fertilization <- NA #create new column to categorize fertilization methods

mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="external"]  <- "External"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="internal (oviduct)"] <- "Internal"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="in brood pouch or similar structure"] <- "Internal" #convert "in brood pouch or similar structure" to internal fertilization

mtdna_final_fertilization_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_fertilization) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

#Final mtDNA Reproduction mode

mtdna_data_new$final_reproductionmode <- NA #create new column to categorize reproduction mode

mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

mtdna_final_reproductionmode_He_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_reproductionmode ) & !is.na(mtdna_data_new$He),] #create new table that excludes NA's from columns of interest

##### mtDNA: Pi ##### 

#Final mtDNA Fertilization 

mtdna_data_new$final_fertilization <- NA #create new column to categorize fertilization methods

mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="external"]  <- "External"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="internal (oviduct)"] <- "Internal"
mtdna_data_new$final_fertilization [mtdna_data_new$fertilization =="in brood pouch or similar structure"] <- "Internal" #convert "in brood pouch or similar structure" to internal fertilization

mtdna_final_fertilization_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_fertilization) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

#Final mtDNA Reproduction mode

mtdna_data_new$final_reproductionmode <- NA #create new column to categorize reproduction mode

mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="dioecism"]  <- "Dioecious"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
mtdna_data_new$final_reproductionmode  [mtdna_data_new$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

mtdna_final_reproductionmode_Pi_no.na <- mtdna_data_new[!is.na(mtdna_data_new$final_reproductionmode ) & !is.na(mtdna_data_new$Pi),] #create new table that excludes NA's from columns of interest

##### msat ##### 

#Final msat fertilization
msat_data$final_fertilization <- NA #create new column to categorize fertilization methods

msat_data$final_fertilization [msat_data$fertilization =="external"]  <- "External"
msat_data$final_fertilization [msat_data$fertilization =="internal (oviduct)"] <- "Internal"
msat_data$final_fertilization [msat_data$fertilization =="in brood pouch or similar structure"] <- "Internal" #convert "in brood pouch or similar structure" to internal fertilization

msat_final_fertilization_He_no.na <- msat_data[!is.na(msat_data$final_fertilization) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

#Final msat Reproduction mode
msat_data$final_reproductionmode <- NA #create new column to categorize reproduction mode

msat_data$final_reproductionmode  [msat_data$reproductionmode =="dioecism"]  <- "Dioecious"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protogyny"] <- "Hermaphrodite" #for protogyny, label as "Hermpahrodites"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="protandry"] <- "Hermaphrodite" #for protandry, label as "Hermpahrodites"
msat_data$final_reproductionmode  [msat_data$reproductionmode =="true hermaphroditism"] <- "Hermaphrodite" #for true hermaphroditism, label as "Hermpahrodites"

msat_final_reproductionmode_He_no.na <- msat_data[!is.na(msat_data$final_reproductionmode ) & !is.na(msat_data$He),] #create new table that excludes NA's from columns of interest

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

### msat & mtDNA: He
Fertplot1 <- ggplot(final_fertilization_all) +
  geom_boxplot(aes(x = final_fertilization, y = He, fill=markertype)) + #final fertilization & He box plot
  ggtitle("(A)") + #add plot title
  xlab("Fertilization Method") + ylab("He") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(size=30, margin = margin(t = 20)),
    axis.title.y = element_text(size=30,margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom") +
  labs(fill="Marker Type") +
  scale_fill_manual(values=c("#444444", "#999999"), labels=c('Nuclear', 'Mitochondrial')) 

### mtDNA: Pi
Fertplot2 <- ggplot(mtdna_final_fertilization_Pi_no.na) + 
  geom_boxplot(aes(x = final_fertilization, y = Pi, fill = final_fertilization)) + #final fertilization & He box plot
  ggtitle("(B)") + #add plot title
  xlab("Fertilization Method") + ylab("Pi") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(size=30, margin = margin(t = 20)),
    axis.title.y = element_text(size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  labs(fill = "Fertilization") +
  scale_fill_manual(values=c("#999999", "#999999"))

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

### msat & mtDNA: He
Reproplot1 <- ggplot(reproductionmode_all) + 
  geom_boxplot(aes(x = final_reproductionmode, y = He, fill= markertype)) + #final fertilization & He box plot
  ggtitle("(A)") + #add plot title
  xlab("Reproduction Mode") + ylab("He") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(size=30, margin = margin(t = 20)),
    axis.title.y = element_text(size=30,margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  labs(fill="Marker Type")+
  scale_fill_manual(values=c("#444444", "#999999"), labels=c('Nuclear', 'Mitochondrial'))

### mtDNA: Pi
Reproplot2 <- ggplot(mtdna_final_reproductionmode_Pi_no.na) + 
  geom_boxplot(aes(x = final_reproductionmode, y = Pi, fill = final_reproductionmode)) + #final fertilization & He box plot
  ggtitle("(B)") + #add plot title
  xlab("Reproduction Mode") + ylab("Pi") + #add axis labels
  theme_bw() +
  theme(                                 #specify characteristics of the plot 
    plot.title = element_text(size=28, margin = margin(b = 15)), 
    axis.title.x = element_text(size=30, margin = margin(t = 20)),
    axis.title.y = element_text(size=30, margin = margin(r = 20)),
    text = element_text(size = 28),
    legend.position="bottom")+
  labs(fill = "Reproduction Mode") +
  scale_fill_manual(values=c("#999999", "#999999"))

#####################################################################################
############### Map of where data was collected ############### 
##### msat #####

geogr_data <- map_data('world')

latlon_msat <-msat_data[!duplicated(msat_data[, c('lat','lon')]), c('lat','lon')] #grab unique lat/lon combos

latlon_msat <- latlon_msat[order(latlon$lat,latlon$lon),] #order by lat then lon

write.csv(latlon_msat, file = paste('latlon_msat',Sys.Date(), '.csv', sep= ''), row.names = FALSE) #build map

msat_plot <- ggplot(geogr_data, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightgray", colour = "white") + 
  geom_point(data = latlon_msat, aes(x = lon, y = lat), size = 5.5, inherit.aes = FALSE, alpha = 0.5) +
  geom_rug(data=latlon_msat, aes(x = lon, y = lat), inherit.aes = F, alpha = 0.6, length = unit(15,"pt"))

#msat_plot

msat_plot_annotated <- msat_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + theme_bw() + 
  ggtitle("(A)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 30, color = "black"), axis.title = element_text(size = 30, face = "bold"), 
        plot.title = element_text(size = 28, face = "bold", hjust = 0)) +
  coord_cartesian(xlim = c(-127, -64), ylim = c(15, 51))

##### mtDNA #####

latlon_mtdna <-mtdna_data_new[!duplicated(mtdna_data_new[, c('lat','lon')]), c('lat','lon')] #grab unique lat/lon combos

latlon_mtdna <- latlon_mtdna[order(latlon$lat,latlonn$lon),] #order by lat then lon

write.csv(latlon_mtdna, file = paste('latlon_mtdna',Sys.Date(), '.csv', sep= ''), row.names = FALSE) #build map

mtdna_plot <- ggplot(geogr_data, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightgray", colour = "white") + 
  geom_point(data = latlon_mtdna, aes(x = lon, y = lat), size = 5.5, inherit.aes = FALSE, alpha = 0.5) +
  geom_rug(data=latlon_mtdna, aes(x = lon, y = lat), inherit.aes = F, alpha = 0.6, length = unit(15,"pt"))

#mtdna_plot

mtdna_plot_annotated <- mtdna_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + theme_bw() + 
  ggtitle("(B)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 30, color = "black"), axis.title = element_text(size = 30, face = "bold"), 
        plot.title = element_text(size = 28, face = "bold", hjust = 0))+
  coord_cartesian(xlim = c(-127, -64), ylim = c(15, 51))
