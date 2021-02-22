#Outputs yfpbgcorrected_manuallycircled.csv, rfpbgcorrected_manuallycircled.csv in the ./Combined/ folder
#Subtracts background fl values found in ./AllManualBG/ from the corresponding trap measurements in ./circledcellfl/.  These trap measurements were obtained by manually circling regions of the cells of interest
#Ignore's bg measurement if there are less than 1000 background pixels. After subtraction the bgcorrected table has an NA
#Removes cells for which there are less than 6 fluorescent measurements




setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
library(dplyr)
library(reshape2)
library(stringr)

#first column corresponding to measurement of fl
firstfcbg = 4;
firstfccirc = 6;

#Open up all the files in data, to get one big file of info, YFP, RFP
bgfilestoopen = list.files("./AllManualBG/",pattern="[[:alnum:]]*_all.csv")
bgfilestoopen = paste('./AllManualBG/',bgfilestoopen,sep="")
rawbg = concatenatefiles(bgfilestoopen)

#Remove entries for which there is no cell
condition = rawbg$xy != 'NaN'
rawbg = filter(rawbg,condition)

#Remove entries for which there are less than 1500 bg pixels (In this case, there isn't enough bg to measure')
condition = rawbg$numbgpix >=  1500
rawbg = filter(rawbg,condition)



#Some of the files have empty columns in the data frame
#Open up all the files containing circled cell fl to get one big file of yfp and rfp
circyfptoopen = list.files("./circledcellfl/",pattern="[[:alnum:]]*_circled_meanC2.csv")
circyfptoopen = paste('./circledcellfl/',circyfptoopen,sep="")
circyfp = concatenatefiles(circyfptoopen)

circrfptoopen = list.files("./circledcellfl/",pattern="[[:alnum:]]*_circled_meanC3.csv")
circrfptoopen = paste('./circledcellfl/',circrfptoopen,sep="")
circrfp = concatenatefiles(circrfptoopen)

#Convert the NaNs to NAs
circyfp[is.na(circyfp)]=NA
circrfp[is.na(circrfp)]=NA


#Reshape the rawbg measurements into an array
numpixbg = dcast(rawbg,filenames+xy+trap~time,value.var='numbgpix')
c2bg = dcast(rawbg,filenames+xy+trap~time,value.var='c2avgbgpix')
c3bg = dcast(rawbg,filenames+xy+trap~time,value.var='c3avgbgpix')
bgdate = str_match(numpixbg$filenames,"./AllManualBG/manualbg_aliveafterdox_(\\w+)_all.csv")
bgdate = bgdate[,2]




#Compare the set of cells for which bg and mean roi intensity are measured.  Why are they different

manualdateC2 = str_match(circyfp$filenames,"./circledcellfl/(\\w+)_circled_meanC2.csv")
manualdateC2 = manualdateC2[,2]

manualdateC3 = str_match(circrfp$filenames,"./circledcellfl/(\\w+)_circled_meanC3.csv")
manualdateC3 = manualdateC3[,2]


manualdate = manualdateC2;


#find positions of each row of the info matrix in the bg matrix.
#We remove manually measured cells for which no background has been measured
manualids = paste(manualdate,circyfp$xy,circyfp$trap)
bgids = paste(bgdate,numpixbg$xy,numpixbg$trap)
matches = match(manualids,bgids)
condition = !is.na(matches)
circyfp = filter(circyfp,condition)
circrfp = filter(circrfp,condition)
matches = matches[condition]
c2bg = c2bg[matches,]
c3bg = c3bg[matches,]
numpixbg = numpixbg[matches,]


#Now subtract bg.  If bg is not valid or not measured, the difference will not be known.  

c2only = circyfp[,firstfccirc:ncol(circyfp)]
c3only = circrfp[,firstfccirc:ncol(circrfp)]
flids = circyfp[,1:(firstfccirc-1)]

c2bgonly = c2bg[,firstfcbg:ncol(c2bg)]
c3bgonly = c3bg[,firstfcbg:ncol(c3bg)] 
c2bgonly = c2bgonly[,1:ncol(c2only)];
c3bgonly = c3bgonly[,1:ncol(c3only)]



#We only care about indices up to 28, so get rid of any later measurements
#We also subtract the bg measurements.
yfp = c2only - c2bgonly
rfp = c3only - c3bgonly

yfp = cbind(flids,yfp)
rfp = cbind(flids,rfp)
write.csv(yfp,file='./CombinedData/yfpbgcorrected_manuallycircled.csv',row.names=FALSE)
write.csv(rfp,file='./CombinedData/rfpbgcorrected_manuallycircled.csv',row.names=FALSE)

