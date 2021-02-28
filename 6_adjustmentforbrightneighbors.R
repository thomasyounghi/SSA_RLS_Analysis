#If there are bright neighboring cells near the trapped cell of interest, correct the cell of interest fluorescence measurement by the average fluorescence in a test roi nearby the bright neighboring cell.
#The test roi is a region adjacent to the bright neighboring cell that is not occupied by any cell.
#Problematic neighbor measurements are in './ProblematicNeighborMeasurements'
#If it was not possible to measure a test background roi, no correction is made, even when a bright neighbor is near the cell of interest.
#Corrected trapped cell measurments are saved in './ProcessedFl' as 'yfpcellsnbadj.csv' and 'rfpcellsnbadj.csv'

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

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



#Making a correction matrix to add to the background corrected fluoresence measurements of the cells of interest. The values stored in this correction matrix are (average bg in large rectangular bg region) - (average fluoresence of test roi)
#The correction is only added for measurements where there is a bright cell next to the cell of interest.
#If it was not possible to measure the fluorescence of a test roi, and a bright neighbor is present, no adjustment is made.
yfpcorrections = problemyfps$c2avgbgpix-problemyfps$avgtestyfp;
yfpcorrections[is.na(yfpcorrections)] = 0;
problemyfps = cbind(problemyfps,yfpcorrections)
problemyfps = problemyfps %>% group_by(filenames,xy,trap,time) %>% summarize(flcorrection = min(yfpcorrections))
problemyfps = data.frame(problemyfps)


#Remove any correction  values less than -100, since these correspond to situations in which a bright cell was present in the test roi region.
hist(problemyfps$flcorrection)
problemyfps = filter(problemyfps,problemyfps$flcorrection >= -100)

#Adding a dummy cell with all possible times so the correction matrix has a column for each time point.
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


