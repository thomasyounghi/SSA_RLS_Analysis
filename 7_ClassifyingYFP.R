#Classifying YFP trajectories for cells in in './ProcessedFl/yfpcellsnbadj.csv'
#For each cell and a given YFP cutoff:
#         Last measuremement below this cutoff -> 'lastofftime'
#         Last measurement above this cutoff -> 'lastontime'
#Classification:
#         If there is no lastofftime, the cell is 'alwayson'
#         If there is no lastontime, the cell is 'alwaysoff'
#         If lastontime>lastofftime, the cell is 'turnedon'
#         If lastontime < lastofftime, the cell is 'alwaysoff'


#Classification is performed with respect to 2 cutoffs: 300 AU and 7 AU
#The lower cutoff is more sensitive to cells with low YFP expression but more susceptible to false positives


#Results are saved in'./YFPclassification/yfpclasses.csv'
#A list of cells to check manually is saved in './YFPclassification yfpclass_tocheck.csv'.  This is based on the class wrt the 7AU cutoff differeing from the class wrt the 300AU cutoff.
#Another list of cells to check manually is saved in './YFPclassification/yfpclass_tocheckdelay.csv'  This is for cells that are 'turnedon' on wrt both cutoffs with a > 4 fluorescent measurement delay between the lastofftimes


setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)

outputfolder = './results/yfpclassification/'

info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
yfp = read.csv('./ProcessedFl/yfpcellsnbadj.csv')



#Applying the strong yfp cutoff of 300 
strongco = 300
yfponly = yfp[,6:ncol(yfp)]
yfponly = mattolol(yfponly)
greaterthanco = lapply(yfponly,greaterthan,strongco)

#identify the index of the last measurment that is below or above the strong cutoff
lastofftime = unlist(lapply(greaterthanco,lastfalseindex))
lastontime = unlist(lapply(greaterthanco,lasttrueindex))

#Classifying single-cell YFP trajectories based on the lastofftime and lastontime
yfpclasshc = rep('',nrow(info))
yfpclasshc[is.na(lastofftime)] = 'alwayson'
yfpclasshc[is.na(lastontime)] = 'alwaysoff'
yfpclasshc[lastontime>lastofftime] = 'turnedon'
yfpclasshc[lastontime<lastofftime] = 'alwaysoff'
lastofftimehc = lastofftime;



#Repeating the YFP classfication with respect to a weak YFP cutoff of 7
weakco = 7;
greaterthanco = lapply(yfponly,greaterthan,weakco)

#identify the index of the last measurment that is above or below the strong cutoff
lastofftime = unlist(lapply(greaterthanco,lastfalseindex))
lastontime = unlist(lapply(greaterthanco,lasttrueindex))

##Classifying single-cell YFP trajectories based on the lastofftime and lastontime
yfpclasslc = rep('',nrow(info))
yfpclasslc[is.na(lastofftime)] = 'alwayson'
yfpclasslc[is.na(lastontime)] = 'alwaysoff'
yfpclasslc[lastontime>lastofftime] = 'turnedon'
yfpclasslc[lastontime<lastofftime] = 'alwaysoff'
lastofftimelc = lastofftime;


#Saving the results of the classification
yfpclassresults = data.frame(info$date, info$xy, info$trap,yfpclasshc,lastofftimehc,yfpclasslc,lastofftimelc)
colnames(yfpclassresults)[1:3]=c('date','xy','trap')
write.csv(yfpclassresults,'./YFPclassification/yfpclasses.csv',row.names=FALSE);

#Saving cells for which the high cutoff class differs from the low cutoff class
yfpclass_tocheck = yfpclassresults[yfpclasshc!=yfpclasslc,]
write.csv(yfpclass_tocheck,'./YFPclassification/yfpclass_tocheck.csv',row.names=FALSE)


#Check the cells classified as 'turned on' where the time for fluoresence to increase from below the low cutoff to above the high cutoff was > 4 fluorence intervals
lastofftimediff = lastofftimehc-lastofftimelc
hist(unlist(lastofftimediff[yfpclasshc=='turnedon']))
meanlastofftimediff = mean(unlist(lastofftimediff[yfpclasshc=='turnedon' & !is.na(lastofftimediff)]))
sdlastofftimediff = sd(unlist(lastofftimediff[yfpclasshc=='turnedon' & ! is.na(lastofftimediff)]))
meanlastofftimediff
sdlastofftimediff

condition = (yfpclasshc == 'turnedon') & (lastofftimediff>4)
yfpclass_tocheckdelay = yfpclassresults[condition,]
yfpclass_tocheckdelay = yfpclass_tocheckdelay[!is.na(yfpclass_tocheckdelay$date),]
write.csv(yfpclass_tocheckdelay,'./YFPclassification/yfpclass_tocheckdelay.csv',row.names=FALSE)



