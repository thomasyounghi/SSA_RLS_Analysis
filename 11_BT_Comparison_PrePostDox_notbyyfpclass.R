


setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(cowplot)
library(reshape2)
library(gridExtra)
library(grid)
library(scales)


#figure settings:
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 16,margin = margin(.3,0,0.3,0,"cm")),strip.background=element_rect(fill="grey"))
themes = theme(axis.text=element_text(size=15), axis.title=element_text(size=16),strip.text.x = element_text(size = 16,margin = margin(.3,0,0.3,0,"cm")),strip.background=element_rect(fill="grey"))

ylims = ylim(-0,25)
ebwidth = 0.25;
bwidth = 0.3;
jitterw=0.15
yscaleforbar = scale_y_continuous(expand = c(0,0),limits = c(0,275))
yscale = scale_y_continuous(expand = c(0,0),limits = c(0,550))
xlabel = xlab('age')
ylabel = ylab('budding time (min)')
fillc = "gray66"
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))


outfolder = './figures/avgbt_notbyyfp/'
outfolder1 = './figures/avgbt_notbyyfp_forpaper/'
info=read.csv('./CombinedData/info_offatdox_alive5hafter_withavgbt5hpostdox_overlappingbisbefore.csv')

#Making sure the last observation time of all the cells occurs after 9 hours post dox addition
condition = (info$lastobservationtime >= (info$doxtime+9*6) | (info$lastobservation == 'burst'))
info = filter(info,condition)

#ignore the Rad52, Rad51 ko cells since we never did an old experiment with those
condition = (info$strain != 'yTY164a') & (info$strain != 'yTY165a')
info = filter(info,condition)

info$expage = factor(info$expage,levels=c('young','old'))
info$expstrain = factor(info$expstrain,levels = c('SSAcontrol','SSA','SSAhet','SSA Dnl4ko','SSAdeg','SSAdegcontrol','SSA3xCln2'))
#levels(info$expstrain) = =c('SSAcontrol','SSA','SSAheterology','SSA Dnl4ko')
info$agestrain = factor(info$agestrain,levels= c('SSAcontrol young','SSAcontrol old','SSA young','SSA old','SSAdegcontrol young','SSAdegcontrol old','SSAdeg young','SSAdeg old','SSA Dnl4ko young','SSA Dnl4ko old','SSAhet young','SSAhet old','SSA3xCln2 young','SSA3xCln2 old'))

summary = info %>% group_by(expstrain,expage,strain,replabel,agestrain) %>% summarize(meanbt0to9= mean(avgbt0to9,na.rm=TRUE),meanbtneg4to0 = mean(avgbtneg4to0,na.rm=TRUE))
meanandsem = summary %>% group_by(expstrain,expage,agestrain) %>% summarize(avgbt0to9 = mean(meanbt0to9),sembt0to9 = sd(meanbt0to9)/sqrt(n()),avgbtneg4to0 = mean(meanbtneg4to0),sembtneg4to0 = sd(meanbtneg4to0)/sqrt(n()))

outtable = paste(outfolder,'avgbtinbeforeandafterdoxaddition.csv',sep="")
write.csv(meanandsem,outtable,row.names=FALSE)



##Making separate plots of avgbt0to9 for different strains
expstrain1 = meanandsem$expstrain
levels(expstrain1) = c('SSAcontrol','SSA','SSAheterology',expression(paste('SSA ',italic('dnl4'),Delta))) 
meanandsem1 = data.frame(expstrain1,meanandsem)

expstrain1 = info$expstrain
levels(expstrain1) = c('SSAcontrol','SSA','SSAheterology',expression(paste('SSA ',italic('dnl4'),Delta))) 
info1 = data.frame(expstrain1,info)


strainstoplot = c('SSA','SSAhet','SSA Dnl4ko')
ylabel = ylab('Budding Time (min)')
xlabel = xlab('')
curroutfolder = paste(outfolder,'./avgbt0to9/',sep="")
for(currstrain in strainstoplot[1:3]){
	#Boxplots with raw data overlaid
	currinfo = filter(info1,(expstrain=='SSAcontrol' | expstrain == currstrain))
	plotfilename = paste(curroutfolder,'avgbt0to9_',currstrain,'_boxplot.pdf',sep="")
	p1 = ggplot(currinfo,aes(expage,avgbt0to9))+facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + geom_boxplot(width=bwidth,outlier.shape=NA) + geom_jitter(width = jitterw,colour = "red",alpha=0.7) + themes + labs(colour = 'Number of buds in interval') + xlabel + ylabel+yscale
	ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


	currmeanandsem = filter(meanandsem1,(expstrain=='SSAcontrol' | expstrain == currstrain))
	plotfilename = paste(curroutfolder,'avgbt0to9_',currstrain,'_youngandold_avgwithinexp_meanandsembar.pdf',sep="")
	p1 = ggplot(currmeanandsem,aes(expage,avgbt0to9))+facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + geom_point() + geom_bar(stat = "identity", width = bwidth,fill=fillc) + geom_errorbar(aes(ymin=avgbt0to9-sembt0to9,ymax=avgbt0to9+sembt0to9),width=ebwidth) + ylabel +themes + xlabel+yscaleforbar
	ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
	removeggplottext(p1,plotfilename,4,600,5,4.5)
	
}



#Plotting average budding time for the time period before doxycycline addition
for(currstrain in strainstoplot[1:3]){
	#Boxplots with raw data overlaid
	currinfo = filter(info1,(expstrain=='SSAcontrol' | expstrain == currstrain))
	plotfilename = paste(curroutfolder,'avgbtneg4to0_',currstrain,'_boxplot.pdf',sep="")
	p1 = ggplot(currinfo,aes(expage,avgbtneg4to0))+facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + geom_boxplot(width=bwidth,outlier.shape=NA) + geom_jitter(width = jitterw,colour = "red",alpha=0.7) + themes + labs(colour = 'Number of buds in interval') + xlabel + ylabel+yscale
	ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


	currmeanandsem = filter(meanandsem1,(expstrain=='SSAcontrol' | expstrain == currstrain))
	plotfilename = paste(curroutfolder,'avgbtneg4to0_',currstrain,'_youngandold_avgwithinexp_meanandsembar.pdf',sep="")
	p1 = ggplot(currmeanandsem,aes(expage,avgbtneg4to0))+facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + geom_point() + geom_bar(stat = "identity", width = bwidth,fill=fillc) + geom_errorbar(aes(ymin=avgbtneg4to0-sembtneg4to0,ymax=avgbtneg4to0+sembtneg4to0),width=ebwidth) + ylabel +themes + xlabel+yscaleforbar
	ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
	removeggplottext(p1,plotfilename,4,600,5,4.5)
	
}



##Plotting predox avgbt alongside postdox budding time
meanandsembt0to9 = select(meanandsem1,agestrain,expstrain,expage,avgbt0to9,sembt0to9,expstrain1)
meanandsembtneg4to0 = select(meanandsem1,agestrain,expstrain,expage,avgbtneg4to0,sembtneg4to0,expstrain1)
colnames(meanandsembt0to9)[4:5] = c('mean','sem')
colnames(meanandsembtneg4to0)[4:5] = c('mean','sem')
avgbtwindow = c(rep('0to9',nrow(meanandsembt0to9)),rep('neg4to0',nrow(meanandsembtneg4to0)))
avgbtwindow = factor(avgbtwindow,levels= c('neg4to0','0to9'),labels=c('pre-dox','post-dox'))
meanandsem2 = rbind(meanandsembt0to9,meanandsembtneg4to0)
meanandsem2 = data.frame(avgbtwindow,meanandsem2)


#Plotting pre and post-dox side by side for SSA and SSA control strain.  Average across cells in an experiment
currmeanandsem = filter(meanandsem2,(expstrain=='SSAcontrol'| expstrain =='SSA'))
p1 = ggplot(currmeanandsem,aes(x=expage,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=expage,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
plotfilename = paste(outfolder,'prepostdox_avgbt_SSA_SSAcontrol_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


#Same plot with green and grey as the colors
greengrey = scale_fill_manual(values=c('gray66','green3'))

hidefacet = theme(strip.background = element_blank(),strip.text.x= element_blank())
newyscale = scale_y_continuous(expand = c(0,0),limits = c(0,300))
p1 = ggplot(currmeanandsem,aes(x=agestrain,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=agestrain,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
p1 = p1 + greengrey + border
plotfilename = paste(outfolder1,'prepostdox_avgbt_SSA_SSAcontrol_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=3.5)
p1 = p1 + hidefacet + newyscale
removeggplottext(p1,plotfilename,4,600,6,3.5)






#Plotting pre and post-dox side by side for SSAheterology and SSA strains.  Average across cells in an experiment
currmeanandsem = filter(meanandsem2,(expstrain=='SSAhet'| expstrain =='SSA'))
p1 = ggplot(currmeanandsem,aes(x=expage,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=expage,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
plotfilename = paste(outfolder,'prepostdox_avgbt_SSA_SSAhet_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)

p1 = ggplot(currmeanandsem,aes(x=agestrain,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=agestrain,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
p1 = p1 + greengrey + border
plotfilename = paste(outfolder1,'prepostdox_avgbt_SSA_SSAhet_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=3.5)
p1 = p1 + hidefacet + newyscale
removeggplottext(p1,plotfilename,4,600,6,3.5)




#Plotting pre and post-dox side by side for SSA Dnl4ko and SSA strains.  Average across cells in an experiment
currmeanandsem = filter(meanandsem2,(expstrain=='SSA Dnl4ko'| expstrain =='SSA'))
p1 = ggplot(currmeanandsem,aes(x=expage,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=expage,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
plotfilename = paste(outfolder,'prepostdox_avgbt_SSA_SSADnl4ko_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)

p1 = ggplot(currmeanandsem,aes(x=agestrain,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=agestrain,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
p1 = p1 + greengrey + border
plotfilename = paste(outfolder1,'prepostdox_avgbt_SSA_SSADnl4ko_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=3.5)
p1 = p1 + hidefacet + newyscale
removeggplottext(p1,plotfilename,4,600,6,3.5)

	
	
#Plotting pre and post-dox side by side for SSA Dnl4ko and SSA strains.  Average across cells in an experiment
currmeanandsem = filter(meanandsem2,(expstrain=='SSAdeg'| expstrain =='SSA'))
p1 = ggplot(currmeanandsem,aes(x=expage,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=expage,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
plotfilename = paste(outfolder,'prepostdox_avgbt_SSA_SSAdeg_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)

p1 = ggplot(currmeanandsem,aes(x=agestrain,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=agestrain,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
p1 = p1 + greengrey + border
plotfilename = paste(outfolder1,'prepostdox_avgbt_SSA_SSAdeg_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=3.5)
p1 = p1 + hidefacet + newyscale
removeggplottext(p1,plotfilename,4,600,6,3.5)




#Plotting pre and post-dox side by side for SSA Dnl4ko and SSA strains.  Average across cells in an experiment
currmeanandsem = filter(meanandsem2,(expstrain=='SSAdeg'| expstrain =='SSA'))
p1 = ggplot(currmeanandsem,aes(x=expage,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=expage,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
plotfilename = paste(outfolder,'prepostdox_avgbt_SSA_SSAdeg_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)

p1 = ggplot(currmeanandsem,aes(x=agestrain,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=agestrain,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
p1 = p1 + greengrey + border
plotfilename = paste(outfolder1,'prepostdox_avgbt_SSA_SSAdeg_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=3.5)
p1 = p1 + hidefacet + newyscale
removeggplottext(p1,plotfilename,4,600,6,3.5)




#Plotting pre and post-do side by side for SSAcontrol and SSAdegcontrol strains.
currmeanandsem = filter(meanandsem2,(expstrain=='SSAdegcontrol'| expstrain =='SSAcontrol'))
currmeanandsem = filter(currmeanandsem,expage=='old')
p1 = ggplot(currmeanandsem,aes(x=expage,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=expage,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + facet_grid(.~expstrain1,scales="free_x",labeller = label_parsed) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
plotfilename = paste(outfolder,'prepostdox_avgbt_SSAcontrol_SSAdegcontrol_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)
removeggplottext(p1,plotfilename,4,600,6,4.5)

p1 = ggplot(currmeanandsem,aes(x=agestrain,y=mean,fill=avgbtwindow)) + geom_bar(stat="identity",width=bwidth,position="dodge") + geom_errorbar(aes(x=agestrain,ymin=mean-sem,ymax=mean+sem),width=ebwidth,position=position_dodge(width=0.3)) + yscaleforbar + labs(fill = "time window") + ylabel +themes + xlabel
p1 = p1 + greengrey + border
plotfilename = paste(outfolder1,'prepostdox_avgbt_SSAcontrol_SSAdegcontrol_meansembar.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4)
p1 = p1 + hidefacet + newyscale
removeggplottext(p1,plotfilename,4,600,5,4)







