#Purpose of the this script is to merge the fluorescence measurements obtained by manually circling cells and automatically circling cells

#we need to specify how the columns of the two 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
yfpmanual = read.csv('./CombinedData/yfpbgcorrected_manuallycircled.csv')
rfpmanual = read.csv('./CombinedData/rfpbgcorrected_manuallycircled.csv')
yfpauto = read.csv('./CombinedData/yfpbgcorrected.csv')
rfpauto = read.csv('./CombinedData/rfpbgcorrected.csv')

date = str_match(yfpmanual$filenames,"./circledcellfl/(\\w+)_circled_meanC2.csv")
date = date[,2]
manualids = cbind(date,yfpmanual[,2:3])
manualids = paste(manualids$date,manualids$xy,manualids$trap)

#The data frame containing the merged manually circled and automatically measured measurements
firstfcauto = 6;
rfpcombined = rfpauto[,firstfcauto:ncol(rfpauto)]
yfpcombined = yfpauto[,firstfcauto:ncol(rfpauto)]
autoids = paste(yfpauto$date,yfpauto$xy,yfpauto$trap)

#Step through each row of manual ids, find the corresponding row of info
matches = NA;
firstfcmanual = 6;
yfpmanual = yfpmanual[,firstfcmanual:ncol(yfpmanual)]
rfpmanual = rfpmanual[,firstfcmanual:ncol(rfpmanual)]

matches = match(manualids,autoids)

#If a manually measured cell has no corresponding automatically measured cell, we can't replace anything
positioninmanual = which(!is.na(matches))
matches = matches[!is.na(matches)]


for(i in 1:length(matches)){
	yfpvalues = yfpmanual[positioninmanual[i],]
	rfpvalues = rfpmanual[positioninmanual[i],]
	measuredi = which(!is.na(yfpvalues))
	matchedrow = matches[i];
	yfpcombined[matchedrow,measuredi] =  yfpvalues[measuredi]
	rfpcombined[matchedrow,measuredi] = rfpvalues[measuredi]
}


#Now we remove any fluorescent measurements occuring before the birth of the cell.  Such measurements are made in error
#Removing invalid measurements for each cell (based when the cells are marked as burst or interfered with
flstart = mapply(firstflindexafter,info$birth,info$flfreq,info$lastibeforepart2+1)




rfp = mattolol(rfpcombined)
rfp = mapply(removeleftvalues,rfp,flstart)
rfp = t(rfp)

yfp = mattolol(yfpcombined)
yfp = mapply(removeleftvalues,yfp,flstart)
yfp = t(yfp)




#YFP combined is formed by filling in the automatic fluorescence measurements with the manually measured ones.
yfpcombined = cbind(yfpauto[,1:(firstfcauto-1)],yfp)
rfpcombined = cbind(rfpauto[,1:(firstfcauto-1)],rfp)


#Saving the merged cell fluorescence data
#Saving the manual fluorescence files with additional cell id information
write.csv(yfpcombined,'./ProcessedFl/yfpcells.csv',row.names=FALSE)
write.csv(rfpcombined,'./ProcessedFl/rfpcells.csv',row.names=FALSE)


