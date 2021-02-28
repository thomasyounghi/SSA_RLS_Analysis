#Determining the cells for which there are fluorescence measurements. 
#This is based on the time that the traps are fully occupied and the last observation time
#Cells for which there are no fluorescence measurements, young cells born
#before time 50, and cells that fail to divide prior to dox removal and 
#not removed
#Remaining cell data is saved in "./CombinedData/infoall2.csv", "./#CombinedData/yfpall2.csv",CombinedData/rfpall2.csv"

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

info = read.csv('./CombinedData/infoall.csv')
yfp = read.csv('./CombinedData/yfpall.csv')
rfp = read.csv('./CombinedData/rfpall.csv')

##Determining whether there are any fluoresence measurements to extract
##We are being strict by using info$fullyoccupied as the lower cutoff
flstart = mapply(firstflindexafter,info$fullyoccupied,info$flfreq,info$lastibeforepart2+1)
flstart[flstart<1] = 1;
flend = mapply(lastflindexbefore,info$lastobservationtime,info$flfreq,info$lastibeforepart2+1)

##condition for there being fluoresence measurements for a cell
numflmeasurements = max(ncol(yfp),ncol(rfp))-2
condition = ((flend-flstart)>=3) & (flend >= 4)
info = filter(info, condition)
yfp = filter(yfp, condition)
rfp = filter(rfp, condition)

#Also make sure young cells are born before time point 50.
condition = (info$birth < 50 & info$doxtime<50) | (info$doxtime>100);
info = filter(info,condition)
yfp = filter(yfp,condition)
rfp = filter(rfp,condition)

#Make sure young cells divide at least once before dox removal
bt = getbudtimes(info)
ageatdoxremoval = ageattimes(bt,info$doxtime+24)
condition = ageatdoxremoval >= 1;
info = filter(info,condition);
yfp = filter(yfp,condition);
rfp = filter(rfp,condition);

write.csv(info,"./CombinedData/infoall2.csv",row.names=FALSE);
write.csv(yfp,"./CombinedData/yfpall2.csv",row.names=FALSE);
write.csv(rfp,"./CombinedData/rfpall2.csv",row.names=FALSE);
