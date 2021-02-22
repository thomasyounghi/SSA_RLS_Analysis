#Comparing cutting at 7h between young and old cells using a single time point
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(ggplot2)
library(reshape2)

info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
rfp = read.csv('./ProcessedFl/rfpcells.csv')


#filter out the control cells
#condition = info$strain != 'yTY126';
#info = filter(info,condition);
#rfp = filter(rfp,condition);


#filter out young cells less than 5 generations and old cells < 15 generations at the time of doxycycline exposure
bt = getbudtimes(info)
ageatdox = ageattimes(bt,info$doxtime)
condition = (ageatdox<=5 & info$doxtime==43) | (ageatdox>=15 & info$doxtime==169)
info = filter(info,condition)
rfp = filter(rfp,condition)
rfponly = rfp[,-c(1:4)]

#filter out cells that


#Find the maximum RFP level for each cell
maxrfp = apply(rfponly,1,max,na.rm=TRUE);

#Last RFP measurement 
rfpat7h = rfponly[,28];

#Difference from max
rfpdrop = maxrfp-rfpat7h;
rfpcut = (rfpdrop > 120);





#We also only consider cells that are still observable/alive at the 7 h time point
condition = (info$lastobservationtime >= info$doxtime + 4*6+7*6 )
info = filter(info,condition)
rfp = filter(rfp,condition)
plotfilename = './SSAresults/SSAhetrepairfractions_cellsalive7hpostdox.tiff'
print(plotfilename)
tiff(plotfilename, units="in", width=6, height=5, res=300)
print(p1)
dev.off()


