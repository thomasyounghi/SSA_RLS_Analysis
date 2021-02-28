#Plot YFP and RFP in the SSAcontrol strain prior to doxycycline treatment
#Compare the mean and sd of the YFP in the SSAcontrol strain to the YFP cutoff used to assess SSA repair
#Figures and tables are saved in ./figures/SSAcontrol/


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
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
themes = theme(axis.text=element_text(size=18), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
ebwidth = 0.25;
bwidth = 0.3
yscale = scale_y_continuous(expand = c(0,0),limits = c(0,10.2))
fillc = "gray66"

info = read.csv('./CombinedData/info_sameasmanualbgflcorrected.csv')
rfp = read.csv('./ProcessedFl/rfpcellsnbadj.csv')
yfp = read.csv('./ProcessedFl/yfpcellsnbadj.csv')


condition = info$strain == 'yTY126';
info = filter(info,condition);
rfp = filter(rfp,condition);
yfp = filter(yfp,condition);
yfponly = yfp[,-c(1:4)]
rfponly = rfp[,-c(1:4)]


#Calculate mean and sd for pre-dox SSAcontrol YFP measurements
yfppredox = yfponly[,1:6]
yfpbyexp = split(yfppredox,info$experiment)
unlistedyfpbyexp = lapply(yfpbyexp,unlist)
meanyfpbyexp = unlist(lapply(unlistedyfpbyexp,mean,na.rm=TRUE))
sdyfpbyexp = unlist(lapply(unlistedyfpbyexp,sd,na.rm=TRUE))
expnames = names(yfpbyexp)
cutoffs = meanyfpbyexp+5*sdyfpbyexp
yfpcontrolpredox = rbind(meanyfpbyexp,sdyfpbyexp,cutoffs)
write.csv(yfpcontrolpredox,'./figures/SSAcontrol/yfpcontrolpredox.csv')


#Plotting YFP mean +/- SD for the SSAcontrol experiments.
yfpcontrolpredox = yfpcontrolpredox[,5:8]
tograph = data.frame(t(yfpcontrolpredox))
agestrain = c('old r1','old r2','young r1','young r2')
tograph = data.frame(tograph,agestrain)
tograph$agestrain = factor(tograph$agestrain,levels = c('young r1','young r2','old r1','old r2'))
p1 = ggplot(tograph,aes(agestrain,meanyfpbyexp)) + geom_point() + geom_bar(width = bwidth,stat="identity",fill=fillc) + geom_errorbar(aes(ymin=meanyfpbyexp-sdyfpbyexp,ymax= meanyfpbyexp+sdyfpbyexp),width=ebwidth) + themes + ylab('YFP (AU)') + xlab(' replicate') + ylim(0,10) + geom_hline(yintercept = 7) + yscale+border
plotfilename = './figures/SSAcontrol/SSAcontrolYFP.pdf'
removeggplottext(p1,plotfilename,4,600,6,4.5)

p1 = p1 + annotate("text",x=2,y=7.4,label = "Cutoff = 7 AU")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)

