#Here we take a look at the single cell YFP classification data.
#For cells alive at 5 h after doxycycline removal we compute the fraction repaired by experiment and save this in a table.  This will be used look at the experimental variability in repair fraction.
#For each experiment, there is a column for number of cells, number of cells repaired based on appearance of YFP,  

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(ggplot2)
library(reshape2)


outputfolder = './results/forpresentation/'
info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
yfp = read.csv('./YFPclassification_includinglatermeasurementsforsomecells/yfpclasses_final2.csv')
yfpcells = read.csv('./ProcessedFl/yfpcellsnbadj.csv')
rfpcells = read.csv('./ProcessedFl/rfpcellsnbadj.csv')


#Only looking at the old cells that are age 15 or older at the time of doxycycline addition
bt = getbudtimes(info)
ageatdox = ageattimes(bt,info$doxtime)
condition = (ageatdox<=5 & info$doxtime==43) | (ageatdox>=15 & info$doxtime==169)
info = filter(info,condition)
yfp =filter(yfp,condition)
yfpcells = filter(yfpcells,condition)
rfpcells = filter(rfpcells,condition)

#Only looking at cells that are observable at 4+5 = 9 h after the addition of doxycycline
#We ignore cells that both burst within this the time window or are lost to observation.
condition = info$lastobservationtime >= (info$doxtime+9*6)
info = filter(info,condition)
yfp = filter(yfp, condition)
yfpcells = filter(yfpcells,condition)
rfpcells = filter(rfpcells,condition)


#Removing cells that are always on
yfpat5h = yfpclassattime(yfp,24)
condition = (yfpat5h$yfpclass != 'alwayson')
info = filter(info,condition)
yfp = filter(yfp,condition)
yfpat5h = yfpclassattime(yfp,24)
yfpcells = filter(yfpcells,condition)
rfpcells = filter(rfpcells,condition)


#Counting number of repaired cells in each experiment
#Easier to read labels for condition of the experiment
colnames(yfpat5h) = c('yfpclass9hpdx','lastofftime')
yfpclass = yfp$yfpclass
cells = data.frame(yfpat5h,yfpclass,info)

#We remove any cell that was observed to turn on before the addition of doxycycline
#For such cells, the last offtime is less than 6.
#Dox was added right before t= 15, which corresponds to a flindex = 1 + 14/3 = 5.66
#The first fl measurement after dox has fl index 6
#As of 8/8/18, the only cell that turned on before the addition of dox was 
#The cell with id 708, agestrain: SSA Dnl4ko young, date: 7_5_18 
#After removing this cell, only the SSA Dnl4ko results will be changed.
hist(yfpat5h$lastofftime)
condition = cells$lastofftime >= 6;
cells = filter(cells, condition);
yfpcells = filter(yfpcells, condition);
rfpcells = filter(rfpcells, condition);

#Removing cells from the replicate with fewer than 10 cells
countsbyreplicate <- cells %>% group_by(date,strain,replabel) %>% summarize(total = n())
experimentnames = paste(countsbyreplicate$date,countsbyreplicate$strain,countsbyreplicate$replabel)
countsbyreplicate = data.frame(countsbyreplicate,experimentnames)
countsbyreplicate = filter(countsbyreplicate,total<10)
cellexpnames = paste(cells$date,cells$strain,cells$replabel)
matches = match(cellexpnames,countsbyreplicate$experimentnames)
condition = is.na(matches)
cells = filter(cells,condition)
yfpcells = filter(yfpcells,condition)
rfpcells = filter(rfpcells,condition)

write.csv(cells,'./CombinedData/info_offatdox_alive5hafter.csv',row.names=FALSE)
write.csv(yfpcells,'./CombinedData/yfpcells_offatdox_alive5hafter.csv',row.names=FALSE)
write.csv(rfpcells,'./CombinedData/rfpcells_offatdox_alive5hafter.csv',row.names=FALSE)



#Counts by replicate;
stats <- cells %>% group_by(strain,doxtime,replabel) %>% summarize(nrepair = sum(yfpclass9hpdx=='turnedon'),total = n())
stats <- mutate(stats,fracrepaired = nrepair/total)
write.csv(stats,'./Tables/yfprepair_outoftotal.csv')







