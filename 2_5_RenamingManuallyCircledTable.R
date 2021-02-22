#Correcting columns names for files in the './circcells_headertofix/' folder so that the column names corresponding to different time measurements are integers.
#Results are saved in the './circledcellfl'

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis')

circfilestoopen = c('7_5_18_circled_meanC2.csv','7_5_18_circled_meanC3.csv','3_23_18_circled_meanC2.csv','3_23_18_circled_meanC3.csv','3_28_19_circled_meanC2.csv','3_28_19_circled_meanC3.csv')

for(i in circfilestoopen){
	currenttab = read.csv(paste('./circcells_headertofix/',i,sep=""),header=TRUE)
	colnames(currenttab)[3:ncol(currenttab)]= 1:(ncol(currenttab)-2)
	write.csv(currenttab,paste('./circledcellfl/',i,sep=""),row.names=FALSE)
}