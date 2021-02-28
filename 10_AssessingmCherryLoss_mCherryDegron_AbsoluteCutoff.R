#Assessing RFP loss in the SSA strain containing the degron tagged mCherry
#Results are saved in './figures/rfpdegron_cuttingefficiency/'


setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)

#figure settings:
theme_set(theme_cowplot())
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))
themes = theme(axis.text=element_text(size=19), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))

ylims = ylim(-0,25)
jitterw = 0.12
ebwidth = 0.35;
ylabel = ylab("age (generations)");
xscale = scale_x_discrete(expand = c(0.6,0))


#folder to save in
outfolder = './figures/rfpdegron_cuttingefficiency/'





#As our cutoff, we use the 95th percentile of rfp measurements for cells that are always observed to be yfp+
info = read.csv('./CombinedData/info_offatdox_alive5hafter.csv')
info$replabel = factor(info$replabel)
rfp = read.csv('./CombinedData/rfpcells_offatdox_alive5hafter.csv')
yfp = read.csv('./CombinedData/yfpcells_offatdox_alive5hafter.csv')
rfpcutoffs = read.csv('./figures/rfpdeg_trajectories_yfpalwayson/pooledrfpdeg_meansd.csv')
rfpco = rfpcutoffs[rfpcutoffs$expage=='all',]$qtl95

#We only look at the degron containing strains
condition = (info$strain=='yTY146') | (info$strain=='yTY147')
info = filter(info,condition)
rfp = filter(rfp,condition)
yfp = filter(yfp,condition)


info$agestrain = factor(info$agestrain,levels = c('SSAdeg young','SSAdeg old','SSAdegcontrol old'),labels=c('SSA+\nRFPdegron\nYoung','SSA+\nRFPdegron\nOld','SSAcontrol+\nRFPdegron\nOld'))


#filter out young cells less than 5 generations and old cells < 15 generations at the time of doxycycline exposure
bt = getbudtimes(info)
ageatdox = ageattimes(bt,info$doxtime)
condition = (ageatdox<=5 & info$doxtime==43) | (ageatdox>=15 & info$doxtime==169)
info = filter(info,condition)
rfp = filter(rfp,condition)
yfp = filter(yfp,condition)


#Look at only the cells that are visible for 9 h after addition of doxycycline
condition = info$lastobservationtime >= (info$doxtime+ 6*9)
info = filter(info,condition)
rfp = filter(rfp,condition)
yfp = filter(yfp,condition)


#Get only the RFP measurements
rfponly = rfp[,6:ncol(rfp)]
yfponly = yfp[,6:ncol(yfp)]
rfponly = rfponly[,1:24]
yfponly = yfponly[,1:24]


#For each cell we count the fraction of trailing measurements that are below the cutoff.  
belowco = rfponly<rfpco


#find the last FALSE measurement in a string of binary measurements
lastaboverfpco = apply(belowco,1,lastfalseindex)
info = data.frame(lastaboverfpco,info)

#We measure the number of trailing TRUE or NAs  
lastnonnai <- function(values){
	return(tail(which(!is.na(values)),1))
}

lastmeasurementi = apply(belowco,1,lastnonnai)

#Time between last false measurement and last measurement
numtrailingoff = lastmeasurementi - lastaboverfpco;
timetrailingoff = numtrailingoff/2;
info = data.frame(timetrailingoff,info)
tograph = data.frame(timetrailingoff,lastaboverfpco,numtrailingoff,info) 
p1 = ggplot(tograph,aes(agestrain,timetrailingoff)) + geom_boxplot(width=ebwidth,outlier.shape=NA)+geom_jitter(height=0.1,width=jitterw,colour="red") + ylab('time RFP- (h)') + themes + xlab('')
plotfilename = paste(outfolder,'jitter_timerfptrailingbelow95qtlco.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)
removeggplottext(p1,plotfilename,4,600,6,4.5)



summary <- tograph %>% group_by(agestrain) %>% summarize(meantimetrailingoff = mean(timetrailingoff,na.rm=TRUE),sdtimetrailingoff=sd(timetrailingoff,na.rm=TRUE))
outfilename = paste(outfolder,'meantimetrailingoff.csv',sep="")
write.csv(summary,outfilename)



#We compute statistics to assess the the quality of the 'timetrailingoff' approach
#We are not measuring cutting efficiency here.  Just comparing the distributions 
#of timetrailingoff betwen the two strains
#Looking only at the fraction of cells with a trailing RFP- (at 12 h) longer than 4.5 h
rfpnegover4pt5h_atfli24 =  timetrailingoff > 4.5;
rfpnegleq4pt5h_atfli24 = timetrailingoff <= 4.5;
alwaysbelowco_atfli24 = is.na(numtrailingoff)
notrailingRFP_atfli24 = timetrailingoff == 0;
yestrailingRFP_atfli24 = timetrailingoff >0;
alltrue = rep(TRUE,length(timetrailingoff))

agestrain = info$agestrain

rfptrailingdf = data.frame(agestrain,rfpnegover4pt5h_atfli24,rfpnegleq4pt5h_atfli24,alwaysbelowco_atfli24,notrailingRFP_atfli24, yestrailingRFP_atfli24,alltrue)
rfptrailingdfmelted  = melt(rfptrailingdf, measure.vars = c('rfpnegover4pt5h_atfli24','rfpnegleq4pt5h_atfli24','alwaysbelowco_atfli24','notrailingRFP_atfli24', 'yestrailingRFP_atfli24','alltrue'))
stats = rfptrailingdfmelted%>%group_by(agestrain,variable) %>% summarize(count=sum(value,na.rm=TRUE))
tabrfptrailing = dcast(stats, agestrain~variable)
colnames(tabrfptrailing)[colnames(tabrfptrailing)=='alltrue'] = 'totalcells'
outfn = paste(outfolder,'rfptrailingstats_byagestrain_fli24.csv',sep="")
write.csv(tabrfptrailing,outfn)



##Cutting Efficiency calculation.
#First ignore those cells that are always below the cutoff. According to our definition, these cells have already been 'cut' and lost RFP expression
condition = !is.na(lastaboverfpco)
info = filter(info,condition)
rfponly = filter(rfponly,condition)
yfponly = filter(yfponly,condition)


#Getting the fraction of cells initially YFP- that show either YFP repair, or a drop in YFP for a period greater than 4.5 h.
rfplost = info$timetrailingoff>4.5;
info = data.frame(rfplost,info)
reporcut = info %>% group_by(agestrain) %>% summarize(nrfplost = sum(rfplost),nrepaired = sum(yfpclass=='turnedon'),nrfplost_or_yfprepaired = sum(rfplost | yfpclass9hpdx=='turnedon'),total=n(),nsomerfpaboveandyfpinitbelow = sum(!is.na(timetrailingoff) & yfpclass9hpdx !='alwayson'))
outfilename = paste(outfolder,'reporcut_timetrailingoffover4pt5horyfprep.csv',sep="")
write.csv(reporcut,outfilename,row.names=TRUE)


#Comparing RFP levels before and after the drop.  We partition the measurements to the trailing RFPoff measurements and the preceding time.  Compute the mean for each time frame.  
condition = info$lastaboverfpco < ncol(rfponly);
before = subsetrows(rfponly,rep(1,nrow(rfponly)),info$lastaboverfpco)
after = subsetrows(rfponly[condition,],info[condition,]$lastaboverfpco+1,lastmeasurementi)
meanbefore = apply(before,1,mean,na.rm=TRUE)
meanaftertemp = apply(after,1,mean,na.rm=TRUE)

meanafter = rep(NA,nrow(info))
meanafter[condition] = meanaftertemp
rfpchange = meanafter-meanbefore;


info = data.frame(rfpchange,meanbefore,meanafter,info)
p1 = ggplot(info,aes(agestrain,rfpchange)) + geom_boxplot(width=ebwidth,outlier.shape=NA)+ geom_jitter(height=0,width=jitterw,colour="red") + ylab('Change in RFP (AU)') + themes + xlab('')
plotfilename = paste(outfolder,'jitter_rfpchange_wrtabsoluterfpcotime.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)
removeggplottext(p1,plotfilename,4,600,6,4.5)



#To address the concern that RFP is low in the SSA mCherry strain irrespective of when doxycycline is added. If this is low, we can't be sure that the long length of trailing background level measurements is do to loss of mCherry, or mCherry signal simply being too weak to begin with.
p1 = ggplot(info,aes(agestrain,meanbefore)) + geom_boxplot(width=ebwidth,outlier.shape=NA)+ geom_jitter(height=0,width=jitterw,colour="red") + ylab('RFP (AU)') + themes + xlab('') + geom_hline(yintercept=rfpco,colour="blue") + ylim(0,80) + themes
plotfilename = paste(outfolder,'jitter_rfpbefore_wrtabsoluterfpcotime.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)
removeggplottext(p1,plotfilename,4,600,6,4.5)


#Now plotting mean RFP after (of course only for cells which show an RFP drop to background)
p1 = ggplot(info,aes(agestrain,meanafter)) + geom_boxplot(width=ebwidth,outlier.shape=NA)+ geom_jitter(height=0,width=jitterw,colour="red") + ylab('RFP (AU)') + themes + xlab('') + themes
plotfilename = paste(outfolder,'jitter_rfpafter_wrtabsoluterfpcotime.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)
removeggplottext(p1,plotfilename,4,600,6,4.5)


summary <- info %>% group_by(agestrain) %>% summarize(meanrfpbefore = mean(meanbefore,na.rm=TRUE),sdrfpbefore=sd(meanbefore,na.rm=TRUE),meanrfpafter=mean(meanafter,na.rm=TRUE),sdrfpafter=sd(meanafter,na.rm=TRUE))
outfilename = paste(outfolder,'meanrfp_relativetodropbelowco.csv',sep="")
write.csv(summary,outfilename)

