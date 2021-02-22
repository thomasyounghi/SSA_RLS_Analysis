#we need to specify how the columns of the two 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(cowplot)
library(reshape2)

#figure settings:
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
themes = theme(axis.text=element_text(size=18), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
ebwidth = 0.25;
bwidth = 0.3
yscale = scale_y_continuous(expand = c(0,0),limits = c(0,10.2))
fillc = "gray66"

info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
rfp = read.csv('./ProcessedFl/rfpcellsnbadj.csv')
yfp = read.csv('./ProcessedFl/yfpcellsnbadj.csv')

#Looking at how often we get large negative values
yfponly = yfp[,6:ncol(yfp)]
rfponly = rfp[,6:ncol(rfp)]
yfpmins = apply(yfponly,1,min,na.rm=TRUE)
negcounts = table(yfpmins<=-10,info$experi)
fracneg = negcounts[2,]/apply(negcounts,2,sum)
barplot(fracneg,ylim=c(0,1),cex.names = 0.5,ylab='fraction',main='Fraction of cells with large negative YFP (< -10) by experiment')
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))




condition = info$strain == 'yTY126';
info = filter(info,condition);
rfp = filter(rfp,condition);
yfp = filter(yfp,condition);
yfponly = yfp[,-c(1:4)]
rfponly = rfp[,-c(1:4)]



#Look at the distribution of YFP measurements prior to doxycycline. Calculate mean and sd for separate experients
yfppredox = yfponly[,1:6]
yfpbyexp = split(yfppredox,info$experiment)
unlistedyfpbyexp = lapply(yfpbyexp,unlist)
meanyfpbyexp = unlist(lapply(unlistedyfpbyexp,mean,na.rm=TRUE))
sdyfpbyexp = unlist(lapply(unlistedyfpbyexp,sd,na.rm=TRUE))
expnames = names(yfpbyexp)
cutoffs = meanyfpbyexp+5*sdyfpbyexp
yfpcontrolpredox = rbind(meanyfpbyexp,sdyfpbyexp,cutoffs)
write.csv(yfpcontrolpredox,'./SSAcontrolresults/yfpcontrolpredox.csv')


#Making a plot of the yfp predox histogram for the presentation
yfpcontrolpredox = yfpcontrolpredox[,5:8]
tograph = data.frame(t(yfpcontrolpredox))
agestrain = c('old r1','old r2','young r1','young r2')
tograph = data.frame(tograph,agestrain)
tograph$agestrain = factor(tograph$agestrain,levels = c('young r1','young r2','old r1','old r2'))
p1 = ggplot(tograph,aes(agestrain,meanyfpbyexp)) + geom_point() + geom_bar(width = bwidth,stat="identity",fill=fillc) + geom_errorbar(aes(ymin=meanyfpbyexp-sdyfpbyexp,ymax= meanyfpbyexp+sdyfpbyexp),width=ebwidth) + themes + ylab('YFP (AU)') + xlab(' replicate') + ylim(0,10) + geom_hline(yintercept = 7) + yscale+border
plotfilename = './figures/SSAcontrolYFP.pdf'
removeggplottext(p1,plotfilename,4,600,6,4.5)

p1 = p1 + annotate("text",x=2,y=7.4,label = "Cutoff = 7 AU")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)





tograph = data.frame(meanyfpbyexp,sdyfpbyexp,names(meanyfpbyexp))
names(tograph) = c('mean','sd','experiment')
tograph = tograph[!is.na(tograph$sd),]
explabels = c('old dox, r1', 'old dox, r2', 'young dox, r1', 'young dox, r2')
p1 = ggplot(tograph,aes(experiment,mean)) + geom_point() + ylim(0,8) + xlab('Experiment') + ylab('Mean +/- SD') + geom_errorbar(aes(ymin = mean-sd,ymax=mean+sd)) + ggtitle('Control SSA strain, YFP for 3 hours prior\n to doxycycline addition')+ scale_x_discrete(labels = explabels)
plotfilename = './SSAcontrolresults/YFPmeansd.tiff'
print(plotfilename)
tiff(plotfilename, units="in", width=6, height=5, res=300)
print(p1)
dev.off()


#Look at the distribution of RFP measurements in the period between removal of doxycycline and 7 hours after removal of doxycycline.  This corresponds to timepoints 13-28
rfp7hpostdox = rfponly[,13:28];
rfpbyexp = split(rfp7hpostdox,info$experiment)
unlistedrfpbyexp = lapply(rfpbyexp,unlist)
meanrfpbyexp = unlist(lapply(unlistedrfpbyexp,mean,na.rm=TRUE))
sdrfpbyexp = unlist(lapply(unlistedrfpbyexp,sd,na.rm=TRUE))
rfpcontrolpostdox = rbind(meanrfpbyexp,sdrfpbyexp)
write.csv(rfpcontrolpostdox,'./SSAcontrolresults/rfpcontrol7hpostdox.csv')

#Plot the mean rfp +/- sd for each control experiment.  This is computed by pooling all measurements in the 7 hours after doxycycline is removed.
tograph = data.frame(meanrfpbyexp,sdrfpbyexp,names(meanrfpbyexp))
names(tograph) = c('meanrfp','sdrfp','experiment')
tograph = tograph[!is.na(tograph$sdrfp),]
p1 = ggplot(tograph,aes(experiment,meanrfp)) + geom_point() + ylim(0,900) + xlab('experiment') + ylab('RFP') + geom_errorbar(aes(ymin = meanrfp-sdrfp,ymax=meanrfp+sdrfp)) + ggtitle('Mean +/- SD of RFP measured for 7 hours after doxycycline removal \n (pooled by experiment)')
p1


#The sd is quite large.  Using a cutoff of several sd would result in a a negative cutoff value to say a strain is repaired.
#Instead, we'll try computing the sd for each cell trajectory individually
rfpsdbycell = apply(rfponly[,13:28],1,sd,na.rm=TRUE);
meanofrfpsdbycell = tapply(rfpsdbycell,info$experiment,mean,na.rm=TRUE);
meanofrfpsdbycell = meanofrfpsdbycell[!is.na(meanofrfpsdbycell)]

sdofrfpsdbycell = tapply(rfpsdbycell,info$experiment,sd,na.rm=TRUE);
sdofrfpsdbycell = sdofrfpsdbycell[!is.na(sdofrfpsdbycell)]
tograph = data.frame(meanofrfpsdbycell,sdofrfpsdbycell)
experiment = row.names(tograph)
tograph = cbind(tograph,experiment)
explabels = c('old dox, r1', 'old dox, r2', 'young dox, r1', 'young dox, r2')
p1 = ggplot(tograph,aes(experiment,meanofrfpsdbycell)) + geom_point() + geom_errorbar(width = 0.1,aes(ymin=meanofrfpsdbycell-sdofrfpsdbycell,ymax = meanofrfpsdbycell+sdofrfpsdbycell))+ylim(0,200)+ ylab("Mean +/- SD") +ggtitle('Control SSA strain, Single cell SD of mCherry over the 7 hours \nafter dox removal, grouped by experiment') + xlab("Experiment") + scale_x_discrete(labels = explabels)
plotfilename = './SSAcontrolresults/mCherrySDbycell.tiff'
print(plotfilename)
tiff(plotfilename, units="in", width=6, height=5, res=300)
print(p1)
dev.off()

write.csv(tograph,'./SSAcontrolresults/singlecellmcsds_summary.csv',row.names=FALSE)