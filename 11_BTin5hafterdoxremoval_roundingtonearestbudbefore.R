#For cells that were unrepaired (YFP-) at the time of doxycycline treatment, calculate the average budding time in the 4 hour time window preceding doxycycline treatment and the 9 hour time window after doxycycline addition.
#Resulting data is saved in './CombinedData/info_offatdox_alive5hafter_withavgbt5hpostdox_overlappingbisbefore.csv'

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

#figure settings:
themes = theme(axis.text=element_text(size=18), axis.title=element_text(size=24),strip.text.x = element_text(size = 18))
ebwidth = 0.25;
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))


info = read.csv('./CombinedData/info_offatdox_alive5hafter.csv');
bt = getbudtimes(info)
bt = mattolol(bt)


#If the cell fails to produce a bud after the end observation, we use the last observation time to calculate the  final budding interval.  
avgbtininterval <- function(bts,start,end,lastobst){
	bts = unlist(bts)
	firsti = tail(which(bts<start),1)
	lasti = tail(which(bts<end),1)
	print(bts)
	print(is.na(firsti))
	print(firsti)
	print(lasti)
	if(length(firsti)==0){
		return(NA)
	}
	endpoints = c(bts[firsti:lasti],lastobst)
	return(mean(diff(endpoints),na.rm=TRUE))
}

countininterval <- function(bts,start,end){
	bts = unlist(bts)
	count = sum((bts >= start) & (bts<end),na.rm=TRUE)
	return(count)
}

highestvaluebelow <- function(bts,end){
	bts = unlist(bts)
	valuesbelow = bts[bts<end];
	if(length(valuesbelow)==0){
		return(NA)
	}
	return(max(valuesbelow,na.rm=TRUE))
}


avgbtonlyininterval <- function(bts,start,end,lastobst){
	bts = unlist(bts)
	firsti = head(which(bts>=start & bts<end),1)
	lasti = tail(which(bts<end),1)
	print(bts)
	print(is.na(firsti))
	print(firsti)
	print(lasti)
	if(length(firsti)==0){
		return(NA)
	}
	if(length(lasti)==0){
		return(NA)
	}
	endpoints = c(bts[firsti:lasti],lastobst)
	return(mean(diff(endpoints),na.rm=TRUE))
}



#Computing the average budding time for the 5 h post dox time window, only rounding to the nearest bud before
#Also computing the number of buds appearing within the time window for reference
start = info$doxtime + 24;
end = info$doxtime + 54;
avgbt = mapply(avgbtininterval,bt,start,end,end)
avgbt4to9 = avgbt*10
nbud4to9 = mapply(countininterval,bt,start,end)


#Computing the average budding time for the 4 h dox time window, only rounding to the nearest bud before
start = info$doxtime;
end = info$doxtime + 24;
avgbt = mapply(avgbtininterval,bt,start,end,end)
avgbt0to4 = avgbt*10
nbud0to4 = mapply(countininterval,bt,start,end)

#Computing the length of partial budding interval ending at the 5 hour post-dox time point. Converted to minutes
end = info$doxtime+ 54
partbibefore9h = end - mapply(highestvaluebelow,bt,end)
partbibefore9h = partbibefore9h*10

#Computing the average budding interval for the 9 h period after dox addition, including only intervals starting after the dox addition #time
start = info$doxtime
end = info$doxtime + 9*6;
avgbt = mapply(avgbtininterval,bt,start,end,end)
avgbt0to9 = avgbt*10

#Computing the average budding interval for the 4 h period prior to dox addition
start = info$doxtime-24
end = info$doxtime;
avgbt = mapply(avgbtininterval,bt,start,end,end)
avgbtneg4to0 = avgbt*10
nbudneg4to0 = mapply(countininterval,bt,start,end)

tosummarize = data.frame(avgbtneg4to0,avgbt0to9,partbibefore9h,nbudneg4to0,nbud0to4,nbud4to9,avgbt0to4,avgbt4to9,info)
write.csv(tosummarize,'./CombinedData/info_offatdox_alive5hafter_withavgbt5hpostdox_overlappingbisbefore.csv',row.names=FALSE)


