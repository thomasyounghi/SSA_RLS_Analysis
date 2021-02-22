#Look at the YFP measurements for single cells in './ProcessedFl/yfpcellsnbadj.csv'
#For each cell and a given YFP cutoff, we find the last measuremement below this cutoff -> lastofftime
#For each cell and a given YFP cutoff, we find the last measurement above this cutoff -> lastontime
#If there is no lastofftime, the cell is 'alwayson'
#If there is no lastontime, the cell is 'always off'
#If there lastontime>lastofftime, the cell is 'turnedon' wrt cutoff
#If lastontime < lastofftime, the cell is 'alwaysoff' wrt cutoff.  The on measurements are probably noise
#Two cutoffs are applied: 300 AU and 7 AU
#Results are saved in'./YFPclassification/yfpclasses.csv'
#A list of cells to check manually is saved in './YFPclassification/yfpclass_tocheck.csv'.  This is based on a cell being 'turnedon' wrt the 7AU cutoff but not the 300AU cutoff
#Another list of cells to check manually is saved in './YFPclassification/yfpclass_tocheckdelay.csv'  This is for cells that are 'turned' on wrt both cutoffs with a > 4 flmeasurement delay between the lastofftimes

#Then we consider cells that are above the lower cutoff of 7 but below the cutoff of 300.
#We output a list of the cells to manually check





#we need to specify how the columns of the two 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(ggplot2)
library(reshape2)

outputfolder = './results/yfpclassification/'


info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
yfp = read.csv('./ProcessedFl/yfpcellsnbadj.csv')


#Applying the strong yfp cutoff.  We identify a measurement index such that it and all later measurements are above the cutoff and all measurements before are below the cutoff.
#We look for trajectories where goes from below the strong co to above the strong co and stays there.
#We find the index of the last measurement below the strong co. 
#We find the index of the first measurement above the strong co.
strongco = 300;
yfponly = yfp[,6:ncol(yfp)]
yfponly = mattolol(yfponly)
greaterthanco = lapply(yfponly,greaterthan,strongco)


#identify the index of the last measurment that is below the strong cutoff
lastofftime = unlist(lapply(greaterthanco,lastfalseindex))
lastontime = unlist(lapply(greaterthanco,lasttrueindex))

#What evidence is sufficient to call a cell truly repaired?
-A trailing string of measurement greater than the strong co.
-show that the last measurement below the strong co is 
yfpclasshc = rep(nrow(info))

#If there are no off measurements, the cell is always on.  Due to missing data, we may not necessarilly know if the cell was off for the unobserved indices.  If there are no on measurements, the cell is always off.
yfpclasshc[is.na(lastofftime)] = 'alwayson'
yfpclasshc[is.na(lastontime)] = 'alwaysoff'
yfpclasshc[lastontime>lastofftime] = 'turnedon'
yfpclasshc[lastontime<lastofftime] = 'alwaysoff'
lastofftimehc = lastofftime;


#Now we apply the same analysis for the weak cutoff
weakco = 7;
greaterthanco = lapply(yfponly,greaterthan,weakco)

#identify the index of the last measurment that is below the strong cutoff
lastofftime = unlist(lapply(greaterthanco,lastfalseindex))
lastontime = unlist(lapply(greaterthanco,lasttrueindex))

#What evidence is sufficient to call a cell truly repaired?
-A trailing string of measurement greater than the strong co.
-show that the last measurement below the strong co is 
yfpclasslc = rep(nrow(info))

#If there are no off measurements, the cell is always on.  Due to missing data, we may not necessarilly know if the cell was off for the unobserved indices.  If there are no on measurements, the cell is always off.
yfpclasslc[is.na(lastofftime)] = 'alwayson'
yfpclasslc[is.na(lastontime)] = 'alwaysoff'
yfpclasslc[lastontime>lastofftime] = 'turnedon'
yfpclasslc[lastontime<lastofftime] = 'alwaysoff'
lastofftimelc = lastofftime;

yfpclassresults = data.frame(info$date, info$xy, info$trap,yfpclasshc,lastofftimehc,yfpclasslc,lastofftimelc)
colnames(yfpclassresults)[1:3]=c('date','xy','trap')
write.csv(yfpclassresults,'./YFPclassification/yfpclasses.csv',row.names=FALSE);

yfpclass_tocheck = yfpclassresults[yfpclasshc!=yfpclasslc,]
write.csv(yfpclass_tocheck,'./YFPclassification/yfpclass_tocheck.csv',row.names=FALSE)

#Check the cells in which the difference between the last off time with the low cutoff and last off time with the high cutoff differ by more than 5 and the cell is repaired based on the high cutoff
lastofftimediff = lastofftimehc-lastofftimelc
hist(unlist(lastofftimediff[yfpclasshc=='turnedon']))
meanlastofftimediff = mean(unlist(lastofftimediff[yfpclasshc=='turnedon']))
sdlastofftimediff = sd(unlist(lastofftimediff[yfpclasshc=='turnedon']))
colastofftimediff = meanlastofftimediff + sdlastofftimediff;

##The lower cutoff last off time might not exist.
condition = (yfpclasshc == 'turnedon') & (lastofftimediff>colastofftimediff) & !is.na(lastofftimediff)

condition = (yfpclasshc == 'turnedon') & (lastofftimediff>4)
yfpclass_tocheckdelay = yfpclassresults[condition,]
write.csv(yfpclass_tocheckdelay,'./YFPclassification/yfpclass_tocheckdelay.csv',row.names=FALSE)



