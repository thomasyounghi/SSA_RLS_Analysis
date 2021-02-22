#Based on the measurements identified in './ProblematicNeighborMeasurements', generates a correction matrix
#for the cells in './ProcessedFl/yfpcells.csv' or './ProcessedFl/rfpcells.csv'
#Columns of the correction matrix correspond to time, rows to cells
#Subtracted values are the means of test background rois in the './ProblematicNeighborMeasurements/' file
#If it was not possible to measure a test background roi, no correction is made, even when a bright neighbor is near the cell of interest.

#we need to specify how the columns of the two 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)


#Adjustment for bright neighboring cells
problemyfps = read.csv('./ProblematicNeighborMeasurements/problematicyfps.csv');
problemrfps = read.csv('./ProblematicNeighborMeasurements/problematicrfps.csv');

yfpcells = read.csv('./ProcessedFl/yfpcells.csv');
rfpcells = read.csv('./ProcessedFl/rfpcells.csv');


#creating a 'dummy' cell with entries for all times, so that the correction matrix has columns for all times
numtimes = 28;
dummycell = data.frame(matrix(ncol=5,nrow=numtimes))
colnames(dummycell) = c('filenames','xy','trap','time','flcorrection')
dummycell$time = 1:numtimes;
dummycell$filenames = 'dummycell'



#We make a correction matrix. If it was not possible to measure the fluorescence of a test roi, the value is set to 0.  Otherwise we add back the avg background and subtract the average of the test roi bgs
yfpcorrections = problemyfps$c2avgbgpix-problemyfps$avgtestyfp;
yfpcorrections[is.na(yfpcorrections)] = 0;
problemyfps = cbind(problemyfps,yfpcorrections)
problemyfps = problemyfps %>% group_by(filenames,xy,trap,time) %>% summarize(flcorrection = min(yfpcorrections))
problemyfps = data.frame(problemyfps)


#We remove any correction  values less than -100, since these probably correspond to situations in which a cell was present in the test roi region.
hist(problemyfps$flcorrection)
problemyfps = filter(problemyfps,problemyfps$flcorrection >= -100)

##Adding a dummy cell with all possible times so the correction matrix has a column for each time point.
problemyfps = rbind(problemyfps,dummycell)

yfpcorrectionmat = dcast(problemyfps,filenames+xy+trap~time,value.var='flcorrection')
yfpcorrectionmat = filter(yfpcorrectionmat,filenames!='dummycell')


yfpcorrectionmat[is.na(yfpcorrectionmat)]=0
date = str_match(yfpcorrectionmat$filenames,"./NeighborMeasurements/(\\w+)_yTY(\\w+).csv")
date = date[,2]
yfpcorrectionmat = cbind(date,yfpcorrectionmat)
yfpcellscorrected = addbgcorrection(yfpcells,yfpcorrectionmat,6,5,28)
write.csv(yfpcellscorrected,'./ProcessedFl/yfpcellsnbadj.csv',row.names=FALSE)


#Now correcting the rfp measurements
rfpcorrections = problemrfps$c3avgbgpix-problemrfps$avgtestrfp;
rfpcorrections[is.na(rfpcorrections)] = 0;
problemrfps = cbind(problemrfps,rfpcorrections)
problemrfps = problemrfps %>% group_by(filenames,xy,trap,time) %>% summarize(flcorrection = min(rfpcorrections))
problemrfps = data.frame(problemrfps)
hist(problemrfps$flcorrection)
problemrfps = filter(problemrfps,problemrfps$flcorrection >= -50)

problemrfps = rbind(problemrfps,dummycell)


rfpcorrectionmat = dcast(problemrfps,filenames+xy+trap~time,value.var='flcorrection')
rfpcorrectionmat = filter(rfpcorrectionmat,filenames !='dummycell')
rfpcorrectionmat[is.na(rfpcorrectionmat)]=0
date = str_match(rfpcorrectionmat$filenames,"./NeighborMeasurements/(\\w+)_yTY(\\w+).csv")
date = date[,2]
rfpcorrectionmat = cbind(date,rfpcorrectionmat)
rfpcellscorrected = addbgcorrection(rfpcells,rfpcorrectionmat,6,5,28);
write.csv(rfpcellscorrected,'./ProcessedFl/rfpcellsnbadj.csv',row.names=FALSE)


