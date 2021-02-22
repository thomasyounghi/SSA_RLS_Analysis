#We have two sets of YFP repair classification based on 2 sets of of YFP measurements
#Here we make sure that we only manually check the cases that are necessary
#Only cells that were always off can be revised to on.
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(ggplot2)
library(reshape2)


yfpclassmanualnodeg = read.csv('./YFPclassification/yfpclass_tocheck_modified_noRFPdeg_before10-2018.csv')
yfpclassmanualdeg = read.csv('./YFPclassification/yfpclass_tocheck_modified_RFPdeg_before10-2018.csv')
yfpclassmanual = rbind(yfpclassmanualnodeg,yfpclassmanualdeg)
write.csv(yfpclassmanual,'./YFPclassification/yfpclass_tocheck_modified_before10-2018.csv',row.names=FALSE)


yfpclasstocheck1 = read.csv('./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheck1.csv')

#YFPclassmanual only consists of cells that were observed to the flindex 25.
#Cells that were classified as unrepaired may actually be observed to be repaired if we follow them for longer.

#Some of the manually checked cells that were classified as unrepaired may turn out to be repaired at later times.
#Since some cells were already manually assessed as repair or not, they do not have to be reassessed.
verifiedon = (yfpclassmanual$lclastoffnexton == 'yes') | (yfpclassmanual$laterlegitontime != 'no')
yfpclassmanualon = yfpclassmanual[verifiedon,]



prevcheckedids  = paste(yfpclassmanualon$date,yfpclassmanualon$xy,yfpclassmanualon$trap)
ids1 = paste(yfpclasstocheck1$date,yfpclasstocheck1$xy,yfpclasstocheck1$trap)

matches = match(prevcheckedids,ids1)
matches = matches[!is.na(matches)]
tocheck = yfpclasstocheck1[-matches,]

#Some of the matches include cells that aren't measured beyond flindex 25. Since these cells are not observed beyond flindex 25, they must have died or been censored
tocheck = filter(tocheck,lastofftimehc>=25)
write.csv(tocheck,'./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheck2.csv',row.names=FALSE)





#Now identifying the cells 
yfpclassdelaynodeg = read.csv('./YFPclassification/yfpclass_tocheckdelay_modified_noRFPdeg_before10-2018.csv')
yfpclassdelaydeg = read.csv('./YFPclassification/yfpclass_tocheckdelay_modified_RFPdeg_before10-2018.csv')
yfpclassdelay = rbind(yfpclassdelaynodeg,yfpclassdelaydeg)
write.csv(yfpclassdelay,'./YFPclassification/yfpclass_tocheckdelay_modified_before10-2018.csv',row.names=FALSE)

yfpclassdelay1 = read.csv('./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheckdelay1.csv')


delay1ids = paste(yfpclassdelay1$date,yfpclassdelay1$xy,yfpclassdelay1$trap)
delayids = paste(yfpclassdelay$date,yfpclassdelay$xy,yfpclassdelay$trap)

#Idenitfy those delay1ids not found in delayids.  These are cells for which a delay between low and high YFP expression was not previously detected.
notprevmanuallychecked = which(is.na(match(delay1ids,delayids)))

#Also identify those cells that are in both lists but for which there are no legitimate YFP positive measurements in the old list (not including later measurements).  We check 
#To see if there are legitimate YFP positive measurements at later time
prevcheckedbutYFPneg = which(is.na(yfpclassdelay$firstlegitontime))

#Saving the cells to check
tocheckdelay = rbind(yfpclassdelay1[notprevmanuallychecked,],yfpclassdelay[prevcheckedbutYFPneg,1:7])
write.csv(tocheckdelay,'./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheckdelay2.csv',row.names=FALSE)

