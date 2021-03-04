#Classifying YFP trajectories based on observation of cells across the full movie 40 fluorescent snapshots. Data in './ProcessedFl/yfpcellsnbadj.csv'
#Previously, (7_ClassifyingYFP.R), only measurements up to the 25th fluorescent snapshot were considered. 

#As before:
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

#Results are saved in './YFPclassification_includinglatermeasurementsforsomecells/'

#Cells to manually check were also detected. Since these overlap with previously manually checked cells in './YFPclassfication/', they need to be filtered out of the new list.
#A list of cells to check manually is saved in 'yfpclass_tocheck2.csv'.  This is based on the class wrt the 7AU cutoff differeing from the class wrt the 300AU cutoff.
#Another list of cells to check manually is saved in 'yfpclass_tocheckdelay2.csv'  This is for cells that are 'turnedon' on wrt both cutoffs with a > 4 fluorescent measurement delay between the lastofftimes


setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)

outputfolder = './YFPclassification_includinglatermeasurementsforsomecells/'


info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
yfp = read.csv('./ProcessedFl/yfpcellsnbadj.csv')


#Applying the strong yfp cutoff of 300 
strongco = 300;
yfponly = yfp[,6:ncol(yfp)]
yfponly = mattolol(yfponly)
greaterthanco = lapply(yfponly,greaterthan,strongco)


#identify the index of the last measurment that is below the strong cutoff
lastofftime = unlist(lapply(greaterthanco,lastfalseindex))
lastontime = unlist(lapply(greaterthanco,lasttrueindex))


#Classifying single-cell YFP trajectories based on the lastofftime and lastontime with the high cutoff
yfpclasshc = rep('',nrow(info))
yfpclasshc[is.na(lastofftime)] = 'alwayson'
yfpclasshc[is.na(lastontime)] = 'alwaysoff'
yfpclasshc[lastontime>lastofftime] = 'turnedon'
yfpclasshc[lastontime<lastofftime] = 'alwaysoff'
lastofftimehc = lastofftime;


#Repeating the YFP classfication with respect to a weak YFP cutoff of 7
weakco = 7;
greaterthanco = lapply(yfponly,greaterthan,weakco)

#identify the index of the last measurment that is below the strong cutoff
lastofftime = unlist(lapply(greaterthanco,lastfalseindex))
lastontime = unlist(lapply(greaterthanco,lasttrueindex))

#Classifying single-cell YFP trajectories based on the lastofftime and lastontime with the low cutoff
yfpclasslc = rep('',nrow(info))
yfpclasslc[is.na(lastofftime)] = 'alwayson'
yfpclasslc[is.na(lastontime)] = 'alwaysoff'
yfpclasslc[lastontime>lastofftime] = 'turnedon'
yfpclasslc[lastontime<lastofftime] = 'alwaysoff'
lastofftimelc = lastofftime;


#Saving the YFP classifications along with all other cell info
yfpclassresults = data.frame(info$date, info$xy, info$trap,yfpclasshc,lastofftimehc,yfpclasslc,lastofftimelc)
colnames(yfpclassresults)[1:3]=c('date','xy','trap')
outfn = paste(outputfolder,'yfpclasses1.csv',sep="")
write.csv(yfpclassresults,outfn,row.names=FALSE);


#The list of cells to check based on the high-cutoff yfp class differing from the low-cutoff yfp class
yfpclasstocheck1 = yfpclassresults[yfpclasshc!=yfpclasslc,]


#Check the cells classified as 'turned on' where the time for fluoresence to increase from below the low cutoff to above the high cutoff was > 4 fluorence intervals
lastofftimediff = lastofftimehc-lastofftimelc
hist(unlist(lastofftimediff[yfpclasshc=='turnedon']))
meanlastofftimediff = mean(unlist(lastofftimediff[yfpclasshc=='turnedon']))
sdlastofftimediff = sd(unlist(lastofftimediff[yfpclasshc=='turnedon']))
colastofftimediff = meanlastofftimediff + sdlastofftimediff;

#Selecting the cells that meet the criteria to be checked based on a slow increase in YFP
condition = (yfpclasshc == 'turnedon') & (lastofftimediff>4) & !(is.na(lastofftimehc)) | (is.na(lastofftimelc) & !is.na(lastofftimehc))
yfpclassdelay1 = yfpclassresults[condition,]




#Narrowing down the cells needing to be manually checked to validate their YFP status after including all time points in the movie.
#Some cells had already been manually checked for YFP fluorescence (but only up to the 25th fluorescent snapshot). Among these cells, only the non-repaired ones need to be further checked. The remainder are removed.
#The prior manually checked data was saved in './YFPclassification/'


#Manual checks that were performed prior to considering all times in the movie
#These were previously flagged for checking by applying the YFP threshold only up to a fluorescent index of 25 (The 25th fluorscent snapshot)
yfpclassmanualnodeg = read.csv('./YFPclassification/yfpclass_tocheck_modified_noRFPdeg_before10-2018.csv')
yfpclassmanualdeg = read.csv('./YFPclassification/yfpclass_tocheck_modified_RFPdeg_before10-2018.csv')
yfpclassmanual = rbind(yfpclassmanualnodeg,yfpclassmanualdeg)
write.csv(yfpclassmanual,'./YFPclassification/yfpclass_tocheck_modified_before10-2018.csv',row.names=FALSE)

#Cells that have already been manually verified to be repaired (when only considering time points before 25), do not need to be rechecked
verifiedon = (yfpclassmanual$lclastoffnexton == 'yes') | (yfpclassmanual$laterlegitontime != 'no')
yfpclassmanualon = yfpclassmanual[verifiedon,]


#Finding the flagged cells to check that have not already been manually verified to be repaired
prevcheckedids  = paste(yfpclassmanualon$date,yfpclassmanualon$xy,yfpclassmanualon$trap)
ids1 = paste(yfpclasstocheck1$date,yfpclasstocheck1$xy,yfpclasstocheck1$trap)

matches = match(prevcheckedids,ids1)
matches = matches[!is.na(matches)]
tocheck = yfpclasstocheck1[-matches,]

#Saving the list of cells to manually check for repair.
tocheck = filter(tocheck,lastofftimehc>=25)
write.csv(tocheck,'./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheck2.csv',row.names=FALSE)




#Manual checks (based on a long delay for cells to pass the fluorescent repair threshold) that were performed prior to considering all times in the movie
#These were previously flagged for checking by applying the YFP threshold only up to a fluorescent index of 25 (The 25th fluorscent snapshot)
yfpclassdelaynodeg = read.csv('./YFPclassification/yfpclass_tocheckdelay_modified_noRFPdeg_before10-2018.csv')
yfpclassdelaydeg = read.csv('./YFPclassification/yfpclass_tocheckdelay_modified_RFPdeg_before10-2018.csv')
yfpclassdelay = rbind(yfpclassdelaynodeg,yfpclassdelaydeg)
write.csv(yfpclassdelay,'./YFPclassification/yfpclass_tocheckdelay_modified_before10-2018.csv',row.names=FALSE)


#Identify those cells flagged to manually check that were not in the previously checked data.  These are cells for which a delay between low and high YFP expression was not previously detected.
delay1ids = paste(yfpclassdelay1$date,yfpclassdelay1$xy,yfpclassdelay1$trap)
delayids = paste(yfpclassdelay$date,yfpclassdelay$xy,yfpclassdelay$trap)
notprevmanuallychecked = which(is.na(match(delay1ids,delayids)))

#Also need to include those cells that were not observed to be YFP positive in the manually checked data (since only times up to a fluorescence index of 25 were considered)
prevcheckedbutYFPneg = which(is.na(yfpclassdelay$firstlegitontime))

#Saving the list of cells to manually check for repair.
tocheckdelay = rbind(yfpclassdelay1[notprevmanuallychecked,],yfpclassdelay[prevcheckedbutYFPneg,1:7])
write.csv(tocheckdelay,'./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheckdelay2.csv',row.names=FALSE)



