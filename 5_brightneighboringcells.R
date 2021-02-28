#This script identifies problematic neighbor cell measurements based on the files in './NeighborMeasurements'
#Neighbor measurements are considered problematic if their mean pixel intensity is > 200 AU above the average pixel intensity measured over the entire local rectangle. The brightness of these neighbors adds to the signal measured for the cell of interest and is not accounted for by background subtraction


#Outputs 2 lists with the problematic measurements in the folder './ProblematicNeighborMeasurements/': 'problematicyfps.csv' and 'problematicrfps.csv' 
#The lists contain the time, and xy location, and trap number for problematic measurements. 


#we need to specify how the columns of the two 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

neighborfilestoopen = list.files("./NeighborMeasurements/",pattern="[[:alnum:]]*neighbormeasurements.csv")
neighborfilestoopen = paste('./NeighborMeasurements/',neighborfilestoopen,sep="")


localbgfilestoopen = list.files("./NeighborMeasurements/",pattern="[[:alnum:]]*neighbgmeasurements.csv")
localbgfilestoopen = paste('./NeighborMeasurements/',localbgfilestoopen,sep="")

neighbors = concatenatefiles(neighborfilestoopen)
localbg = concatenatefiles(localbgfilestoopen)

#extract the date and strain portion of the filename
neighbordates = str_match(neighbors$filenames,"./NeighborMeasurements/(\\w+)_yTY(\\d+)_neighbormeasurements.csv")
neighbordates = neighbordates[,2]
localbgdates = str_match(localbg$filenames,"./NeighborMeasurements/(\\w+)_yTY(\\d+)_neighbgmeasurements.csv")
localbgdates= localbgdates[,2]
localbg$date = localbgdates
neighbors$date = neighbordates


#Take average of local test bgs. Obtain counts as well
bgtestroistats = localbg %>% group_by(date,xy,trap,time,neighborid) %>% summarize(avgtestyfp = mean(yfpmean),avgtestrfp = mean(rfpmean),sdtestyfp = sd(yfpmean),sdtestrfp = sd(rfpmean),testcount = n())

allneighborstats = left_join(neighbors,bgtestroistats,by=c("date","xy","trap","time","neighborid"))


#Read in the manualbgmeasurements (averaged over a larger rectangular area)
#left join to the neighboring cell measurements at the same location and time
rawbg = read.csv('./CombinedData/cleanedconcatrawbg.csv')
allneighborstats = left_join(allneighborstats,rawbg,by=c("date","xy","trap","time"))


#Background subtracting the neighboring cell measurements
allneighborstats = mutate(allneighborstats,yfpabovebg = yfpmean-c2avgbgpix,rfpabovebg = rfpmean-c3avgbgpix)
allneighborstats = allneighborstats[!is.na(allneighborstats$numbgpix),]



#Neighboring cells are problematic if their average fluorescence is > 200 above the bg fluroescence and they are a distance of 15 pixels or less from the center of the cell of interest.
#A measurement that is problematic in yfp might not be problematic in rfp and vice versa.
brightyfpneighbor = (allneighborstats$yfpabovebg>200) & (allneighborstats$disttotrapcenter<15)
brightrfpneighbor = (allneighborstats$rfpabovebg>200) & (allneighborstats$disttotrapcenter<15)

problematicyfps = allneighborstats[brightyfpneighbor,]
problematicrfps = allneighborstats[brightrfpneighbor,]



#Saving the problematic neighbor measurements
write.csv(problematicyfps,'./ProblematicNeighborMeasurements/problematicyfps.csv',row.names=FALSE)
write.csv(problematicrfps,'./ProblematicNeighborMeasurements/problematicrfps.csv',row.names=FALSE)






