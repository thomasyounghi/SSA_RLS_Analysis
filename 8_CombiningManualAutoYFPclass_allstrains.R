#Using the manual YFP classification to correct the automatic YFP classification (cells satisfying the strong cutoff) in './YFPclassification/yfpclasses.csv'
#The manually updated assessments are in './YFPclassification/yfpclass_tocheck_modified.csv' and 'YFPclassification/yfpclass_tocheckdelay_modified.csv'
#Results are saved in './YFPclassification/yfpclasses_final.csv'


#we need to specify how the columns of the two 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(ggplot2)
library(reshape2)

outputfilename= './YFPclassification_includinglatermeasurementsforsomecells/yfpclasses_final2.csv'


manual1 = read.csv('./YFPclassification/yfpclass_tocheck_modified_before10-2018_with3-2019.csv',stringsAsFactors=FALSE)
manual2 = read.csv('./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheck2_modified.csv',stringsAsFactors=FALSE)
auto = read.csv('./YFPclassification_includinglatermeasurementsforsomecells/yfpclasses1.csv',stringsAsFactors=FALSE)


#First update the rows of manual1 with manual2
#Manual2 consists of manual verification of YFP repair for those cells of manual1 that were not observed to be repaired
ids1 = paste(manual1$date,manual1$xy,manual1$trap)
ids2 = paste(manual2$date,manual2$xy,manual2$trap)
rowstochange = match(ids2,ids1)
manual = manual1
manual[rowstochange[!is.na(rowstochange)],]=manual2[!is.na(rowstochange),]

#There are also rows of manual2 that don't appear in manual1. Append these rows
manual = rbind(manual,manual2[is.na(rowstochange),])



#We now modify the auto yfpclassification file with the manually checked information
yfpclassmanual = manual$yfpclasslc
newlastofftimelc = manual$lastofftimelc


#If yfpclasslc is 'turnedon' but lclastoffnexton is 'no' & later legitontime is 'no'  the cell is not repaired
#We  also change the last offtime  to the last index measured.
yfpclassmanual[(manual$yfpclasslc == 'turnedon')&(manual$lclastoffnexton == 'no') & (manual$laterlegitontime == 'no')] = 'alwaysoff';

newlastofftimelc[(manual$yfpclasslc == 'turnedon')&(manual$lclastoffnexton == 'no') & (manual$laterlegitontime == 'no')] = manual$lastofftimehc[(manual$yfpclasslc == 'turnedon')&(manual$lclastoffnexton == 'no') & (manual$laterlegitontime == 'no')]

#If yfpclasslc is 'turnedon' but lclastoffnexton is 'no' but later legitontime is a number, the cell is still repaired but the lastofftimelc needs to be changed to this later number
newlastofftimelc[(manual$yfpclasslc=='turnedon')&(manual$lclastoffnexton =='no')&(manual$laterlegitontime!='no')] = manual$laterlegitontime[(manual$yfpclasslc=='turnedon')&(manual$lclastoffnexton =='no')&(manual$laterlegitontime!='no')]

#if yfpclasslc is 'alwayson' but lclastoffnexton is an integer, or laterlegitontime is an integer, the cell is always on
yfpclassmanual[(manual$yfpclasslc == 'alwayson') & (manual$laterlegitontime==1)] = 'alwayson'
newlastofftimelc[(manual$yfpclasslc == 'alwayson') & (manual$laterlegitontime==1)] = NA



#Now take the manually modified yfpclass and lastofftimes and combine with the automatically measured data
manual = data.frame(manual,yfpclassmanual,newlastofftimelc)

#We incorporate the manually checked yfp measurements into the automatic classification
autoids = paste(auto$date,auto$xy,auto$trap);
manualids = paste(manual$date,manual$xy,manual$trap);


matches = match(manualids,autoids)
table(is.na(matches))
condition = !is.na(matches)
yfpclassmanual = yfpclassmanual[condition]
matches = matches[condition]
newlastofftimelc = newlastofftimelc[condition]

yfpclass = auto$yfpclasslc;
yfpclass[matches] = yfpclassmanual;
lastofftime = auto$lastofftimelc;
lastofftime[matches] = newlastofftimelc;









#Now open the file containing cells in which a delay of > 4 measurements was observed
#yfp between 7 and 300 AU for more than 4 measurements

yfpdelayed = read.csv('./YFPclassification/yfpclass_tocheckdelay_modified_before10-2018_with3-2019.csv');
yfpdelayedids = paste(yfpdelayed$date,yfpdelayed$xy,yfpdelayed$trap)

yfpdelayed2 = read.csv('./YFPclassification_includinglatermeasurementsforsomecells/yfpclass_tocheckdelay2_modified.csv')
yfpdelayed2ids = paste(yfpdelayed2$date,yfpdelayed2$xy,yfpdelayed2$trap)

#We update the rows of yfpdelayed with the rows of yfpdelayed2
matches = match(yfpdelayed2ids,yfpdelayedids)
table(is.na(matches))
rowstoreplace = matches[!is.na(matches)]
if(length(rowstoreplace)>0){
	yfpdelayed[rowstoreplace,]=yfpdelayed2[!is.na(matches),]
}
rowstoappend = yfpdelayed2[is.na(matches),];
yfpdelayed = rbind(yfpdelayed,rowstoappend)



yfpdelayedids = paste(yfpdelayed$date,yfpdelayed$xy,yfpdelayed$trap)
matches = match(yfpdelayedids,autoids)
table(is.na(matches))
condition = !is.na(matches)
lastofftimemanual = yfpdelayed$firstlegitontime-1;
lastofftimemanual = lastofftimemanual[condition]
matches = matches[condition];
lastofftime[matches] = lastofftimemanual;


lastofftime= as.numeric(lastofftime)
yfpclassfinal = data.frame(auto,yfpclass,lastofftime)
write.csv(yfpclassfinal,outputfilename,row.names=FALSE)

