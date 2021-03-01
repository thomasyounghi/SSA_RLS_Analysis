#Calculating and plotting fraction of cells repaired (within replicate) by strain
#Counts of repaired and unrepaired cells for each experiment are saved in a table
#Results are saved in './figures/repairfracbyreplicate/'

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(grid)

outfolder = './figures/repairfracbyreplicate/'

#figure settings:
theme_set(theme_cowplot())
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))
themes = theme(axis.text=element_text(size=17), axis.title=element_text(size=16),strip.text.x = element_text(size = 14),panel.border = element_rect(colour="black",linetype="solid",fill=NA,size=1))
ylims = ylim(-0,25)
ebwidth = 0.25;
#cartcoordx= coord_cartesian(ylim = c(0,1),clip='off')
xscale = scale_x_discrete(expand = c(0.6,0))
yscale = scale_y_continuous(breaks = seq(0,1,0.25),labels=seq(0,100,25),limits = c(0,1.1))
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))



#Counting repaired cells and calculating repair fractions for each replicate
info = read.csv('./CombinedData/info_offatdox_alive5hafter.csv')
info$expstrain = factor(info$expstrain,levels = c('SSAcontrol','SSA','SSAhet','SSA Dnl4ko','SSAdeg','SSAdegcontrol','SSA3xCln2','SSA Rad51ko','SSA Rad52ko'))
info$expage = factor(info$expage, levels = c('young','old'))

yfpbyexp <- info %>% group_by(expstrain,expage,agestrain,date,lane) %>% summarize(nrepair = sum(yfpclass9hpdx == 'turnedon'),total=n())
yfpbyexpfn = paste(outfolder,'repaircounts_9hpdx.csv',sep='')
write.csv(yfpbyexp,yfpbyexpfn,row.names=FALSE)
yfpbyexp <- mutate(yfpbyexp,frepair = nrepair/total)
tograph <- yfpbyexp %>% group_by(expstrain,expage) %>% summarize(meanfrepair = mean(frepair),semfrepair = sd(frepair)/sqrt(n()))



#For including a p-value in the plot and a labeling line
ptext = geom_text(x=1.5,y=1.05,label = "p = 0.003")
hsegment = geom_segment(aes(x=1,y=1,xend=2,yend=1))
#p1 <- p1 + ptext + hsegment
#ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)



#Creating separate plots of repair efficiency for all strains
p1 = ggplot(tograph[tograph$expstrain=='SSA',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes  + xscale + yscale
plotfilename = paste(outfolder,'SSAfracyfprepaired.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


p1 = ggplot(tograph[tograph$expstrain=='SSAhet',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes  + xscale + yscale
plotfilename = paste(outfolder,'SSAhetfracyfprepaired.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


p1 = ggplot(tograph[tograph$expstrain=='SSA Dnl4ko',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSA_Dnl4kofracyfprepaired.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


#A narrower version of the same Dnl4 ko plot
p1 = ggplot(tograph[tograph$expstrain=='SSA Dnl4ko',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSA_Dnl4kofracyfprepaired_narrower.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#SSA with Cln2degron tagged to the RFP
p1 = ggplot(tograph[tograph$expstrain=='SSAdeg',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSAdeg_kofracyfprepaired.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


#A narrower version of the SSAdeg plot
p1 = ggplot(tograph[tograph$expstrain=='SSAdeg',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSAdeg_kofracyfprepaired_narrower.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#SSA with 3xCln2
p1 = ggplot(tograph[tograph$expstrain=='SSA3xCln2',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSA3xCln2_fracyfprepaired.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


#A narrower version of the SSA3xCln2 plot
p1 = ggplot(tograph[tograph$expstrain=='SSA3xCln2',],aes(expage,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSA3xCln2_fracyfprepaired_narrower.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#SSA with Rad52 ko or Rad51 ko
p1 = ggplot(tograph[tograph$expstrain=='SSA Rad52ko' | tograph$expstrain=='SSA Rad51ko',],aes(expstrain,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSARad52_Rad51ko_fracyfprepaired.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)

plotfilename = plotfilename = paste(outfolder,'SSARad52_Rad51ko_fracrepaired_narrower.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#A narrower version of the SSA3xCln2 plot
p1 = ggplot(tograph[tograph$expstrain=='SSA Rad52ko' | tograph$expstrain=='SSA Rad51ko',],aes(expstrain,meanfrepair)) + geom_point()+ geom_errorbar(aes(ymin = meanfrepair-semfrepair,ymax=meanfrepair+semfrepair),width=ebwidth) + ylim(0,1)  +xlab('') + ylab('Repair Efficiency (%)') + ggtitle('')+themes + xscale + yscale
plotfilename = plotfilename = paste(outfolder,'SSARad52_Rad51ko_fracyfprepaired_narrower.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)

