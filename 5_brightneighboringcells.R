#This script identifies problematic neighbor measurements based on the files in './NeighborMeasurements'
#In the folder './ProblematicNeighborMeasurements/' outputs a 2 lists containing problematic yfp and rfp measurements: 'problematicyfps.csv' and 'problematicrfps.csv'.  These list the time, and xy location, and trap for problematic measurements. If it was possible to measure a test region nearby the bright neighbor for background, the mean of these test regions is calculated
#Neighbor measurements are considered problematic if their mean pixel intensity is > 200 AU above the average pixel intensity measured over the entire local rectangle.

#we need to specify how the columns of the two 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
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


#Read in the manualbgmeasurements (averaged over a larger area)
#We want to compare these with the brightness of the neighboring cells to determine whether #The local fluorescence will have 
rawbg = read.csv('./CombinedData/cleanedconcatrawbg.csv')
allneighborstats = left_join(allneighborstats,rawbg,by=c("date","xy","trap","time"))


allneighborstats = mutate(allneighborstats,yfpabovebg = yfpmean-c2avgbgpix,rfpabovebg = rfpmean-c3avgbgpix)
allneighborstats = allneighborstats[!is.na(allneighborstats$numbgpix),]


#We will only concern ourselves with the neighbor if it is very bright > 200 above the bg fluroescence.
#A measurement that is problematic in yfp might not be problematic in rfp and vice versa.
#We define a rule to call a bright neighbor problematic (>200 above the empirically measured bg, and a distance of 15 or less from the center of the roi).


#Even if a neighbor is problematic and there is no test fluorescence measurement, it doesn't mean that we can't still make a decent evaluation of the cell of interest.  If the cell of interest has far higher fluoresence than the neighbor, then it is probably positive for the fluorescence.

brightyfpneighbor = (allneighborstats$yfpabovebg>200) & (allneighborstats$disttotrapcenter<15)
brightrfpneighbor = (allneighborstats$rfpabovebg>200) & (allneighborstats$disttotrapcenter<15)

problematicyfps = allneighborstats[brightyfpneighbor,]
problematicrfps = allneighborstats[brightrfpneighbor,]


#If the trapped cell of interest is just as bright, its not that big a deal, If the trapped cell of interest is middling in brightness, that could be due to the bright neighbor.
write.csv(problematicyfps,'./ProblematicNeighborMeasurements/problematicyfps.csv',row.names=FALSE)
write.csv(problematicrfps,'./ProblematicNeighborMeasurements/problematicrfps.csv',row.names=FALSE)






