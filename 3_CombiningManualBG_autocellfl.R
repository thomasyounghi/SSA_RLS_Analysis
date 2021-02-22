#Outputs yfpbgcorrected.csv, rfpbgcorrected.csv,infosameasmanualbgcorrected.csv in the ./Combined/ folder
#Subtracts background fl values found in ./AllManualBG/ from the corresponding trap measurements in yfpall2.csv,rfpall2.csv
#Ignore's bg measurement if there are less than 1000 background pixels. After subtraction the bgcorrected table has an NA
#Removes cells for which there are less than 6 fluorescent measurements



setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
library(dplyr)
library(reshape2)
library(stringr)
library(cowplot)

#Open up all the files in data, to get one big file of info, YFP, RFP
filestoopen = list.files("./AllManualBG/",pattern="[[:alnum:]]*_all.csv")
filestoopen = paste('./AllManualBG/',filestoopen,sep="")
rawbg = concatenatefiles(filestoopen)

#Add a date column 
rawbgdates = str_match(rawbg$filenames,"./AllManualBG/manualbg_aliveafterdox_(\\w+)_all.csv")
rawbgdates = rawbgdates[,2]
rawbg$date=rawbgdates

#Remove entries for which there is no cell
condition = rawbg$xy != 'NaN'
rawbg = filter(rawbg,condition)
id = paste(rawbg$filenames,rawbg$xy,rawbg$trap)
rawbg = cbind(id,rawbg)

#Identify entries in which the number of bg pixels is equal to the number of pixels in the rectangle.  In these situations, there was mostlikely an error with saving the nonbg mask.
table(rawbg$numbgpix == rawbg$numlocalpix)
rawbg[rawbg$numbgpix == rawbg$numlocalpix,]
condition = (rawbg$numbgpix != rawbg$numlocalpix)
rawbg = filter(rawbg,condition)
hist(unlist(rawbg$numbgpix))
p1 = ggplot(rawbg,aes(time,c2avgbgpix)) + geom_point() + ylim(0,1500)
#A problematic cell measurement after subtracting background (< -500 AU measurement)
rawbg[rawbg$date =='7_11_18' & rawbg$xy == 9 & rawbg$trap==3,]


#Remove entries for which there are less than 1000 bg pixels (In this case, there isn't enough bg to measure')
condition = rawbg$numbgpix >=  1000
rawbg = filter(rawbg,condition)


#Remove background values for times that were not verified.
trapsmeasuredbeyond25 = read.csv('./TrapsToMeasure_oldover15_unrepairedatt24/alloldover15_unrepairedatt24.csv')
namesbeyond25 = paste(trapsmeasuredbeyond25$date,trapsmeasuredbeyond25$xy,trapsmeasuredbeyond25$trap)
namesrawbg = paste(rawbg$date,rawbg$xy,rawbg$trap)
rawbgnotbeyond25 = is.na(match(namesrawbg,namesbeyond25))
toremove = (rawbgnotbeyond25 & rawbg$time>25);
rawbg = filter(rawbg,!toremove)


#save the concatenated rawbg file
write.csv(rawbg[,-c(1,2)],'./CombinedData/cleanedconcatrawbg.csv',row.names=FALSE)


#Open the info, rfp, and yfp folders
info = read.csv('./CombinedData/infoall2.csv')
yfp = read.csv('./CombinedData/yfpall2.csv')
rfp = read.csv('./CombinedData/rfpall2.csv')

#Reshape the rawbg measurements into an array
numpixbg = dcast(rawbg,filenames+xy+trap~time,value.var='numbgpix')
c2bg = dcast(rawbg,filenames+xy+trap~time,value.var='c2avgbgpix')
c3bg = dcast(rawbg,filenames+xy+trap~time,value.var='c3avgbgpix')

bgdate = str_match(numpixbg$filenames,"./AllManualBG/manualbg_aliveafterdox_(\\w+)_all.csv")
bgdate = bgdate[,2]


#Compare the set of cells for which bg and mean roi intensity are measured.  Why are they different
lastobsflindex = mapply(lastflindexbefore,info$lastobservationtime,info$flfreq,info$lastibeforepart2+1)
firstflindexafterdox = ceiling((15-1)/3+1)
condition = lastobsflindex >= firstflindexafterdox;
info = filter(info, condition)
yfp = filter(yfp,condition)
rfp = filter(rfp,condition)

p1 = ggplot(rawbg,aes(time,c2avgbgpix)) + geom_point() + ylim(0,1500)

#find positions of each row of the info matrix in the bg matrix
#For some cells in info, there are no background measurements
infoids = paste(info$date,info$xy,info$trap)
bgids = paste(bgdate,numpixbg$xy,numpixbg$trap)
matches = match(infoids,bgids)

#Removing cells for which there are no background measurements
condition = !is.na(matches)
info = filter(info,condition)
yfp = filter(yfp,condition)
rfp = filter(rfp,condition)
matches = matches[condition]
c2bg = c2bg[matches,]
c3bg = c3bg[matches,]
numpixbg = numpixbg[matches,]


#Checking that the rows of cells and trap backgrounds have been matched appropriately
table(is.na(matches))
table(paste(c2bg$xy,c2bg$trap)==paste(yfp$xy,yfp$trap))
table(paste(c2bg$xy,c2bg$trap)==paste(rfp$xy,rfp$trap))


#Now subtract bg.  If bg is not valid or not measured, the difference will not be known.  And we subtract a row of NAs
firstfcbg = 4
firstfccell = 6
c2bgonly = c2bg[,firstfcbg:ncol(c2bg)]
c3bgonly = c3bg[,firstfcbg:ncol(c3bg)]
c2only = yfp[,firstfccell:ncol(yfp)]
c3only = rfp[,firstfccell:ncol(rfp)]
flids = yfp[,1:(firstfccell-1)]


#We only care about indices up to 28, so get rid of any later measurements
#We also subtract the bg measurements.
maxfli =38
c2bgonly = c2bgonly[,1:maxfli]
c3bgonly = c3bgonly[,1:maxfli]
c2only = c2only[,1:maxfli]
c3only = c3only[,1:maxfli]
yfp = c2only - c2bgonly
rfp = c3only - c3bgonly



#Removing invalid measurements for each cell (based when the cells are marked as burst or interfered with
flstart = mapply(firstflindexafter,info$fullyoccupied,info$flfreq,info$lastibeforepart2+1)
flstart[flstart<1] = 1;
flend = mapply(lastflindexbefore,info$lastobservationtime,info$flfreq,info$lastibeforepart2+1)


rfp = mattolol(rfp)
rfp = mapply(removerightvalues,rfp,flend)
rfp = t(rfp)
rfp = mattolol(rfp)
rfp = mapply(removeleftvalues,rfp,flstart)
rfp = t(rfp)

yfp = mattolol(yfp)
yfp = mapply(removerightvalues,yfp,flend)
yfp = t(yfp)
yfp = mattolol(yfp)
yfp = mapply(removeleftvalues,yfp,flstart)
yfp = t(yfp)


#Find the number of fluorescent measurements for each cell
nummeas = apply(yfp,1,f<-function(x){sum(!is.na(x))})

#add back cell id information
yfp = cbind(flids,yfp)
rfp = cbind(flids,rfp)

#Get rid of cells for which there are less than 1 fluorescent measurements
#This also gets rid of cells for which no background measurements were taken

condition = nummeas > 1
info = filter(info,condition)
yfp = filter(yfp,condition)
rfp = filter(rfp,condition)

#Checking a situation where the bg was an overestimate. After subtraction -63.  The problem is that there is
#a cell present in the YFP channel but not the phaseg channel.
#The bg estimate from taking a test bg is even further off resulting in a measurement of -600 AU.
yfp[info$date=='7_11_18' & info$xy==9 &info$trap ==3,]


write.csv(yfp,file='./CombinedData/yfpbgcorrected.csv',row.names=FALSE)
write.csv(rfp,file='./CombinedData/rfpbgcorrected.csv',row.names=FALSE)
write.csv(info,'./CombinedData/info_sameasmanualbgflcorrected.csv',row.names=FALSE)

