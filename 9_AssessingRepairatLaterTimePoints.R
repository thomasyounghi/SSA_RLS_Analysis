

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')library(dplyr)

library(ggplot2)
library(cowplot)
library(reshape2)
library(grid)


#original ebwidth is 0.25
#figure settings:
theme_set(theme_cowplot())
themes = theme(axis.text=element_text(size=18), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
ylims = ylim(-0,25)
jitterw = 0.12
ebwidth = 0.25;
bwidth = 0.4;
ylabel = ylab("Fraction Repaired");
xscale = scale_x_discrete(expand = c(0.6,0))
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))
fillc = "gray66"
yscale = scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = seq(0,1,0.25))

#folder to save in
outfolder = './figures/repairbeyond9hpostdoxaddition/'
info = read.csv('./CombinedData/info_offatdox_alive5hafter_withavgbt5hpostdox_overlappingbisbefore.csv')
info$expstrain2 = factor(info$expstrain,levels = c('SSAcontrol','SSA','SSA Dnl4ko','SSAdeg','SSAdegcontrol','SSAhet'),labels = c("SSAcontrol",'SSA',expression(paste('SSA ',italic('dnl4'),Delta)),'SSA+RFPdegron','SSAcontrol+RFPdegron',"SSAheterology"))
info$expage = factor(info$expage,levels=c('young','old'))



#ignore control strains
info = filter(info,(strain!='yTY126') & (strain!='yTY146'))

#We want to look at whether cells that were assessed as yfp unrepaired at 9 h post-dox addition are eventually repaired
#If so, that would mean lack of repair in 9 h is really slow repair
#Here we consider both cells that burst in the 5 h window, and cells that are alive at the end of the window.  We include
#the burst cells in case repair correlates with bursting
urinfo = filter(info,yfpclass9hpdx == 'alwaysoff')
urinfo = filter(urinfo,(lastobservationtime > (doxtime + 14*6)) | (lastobservation == 'burst'))

toplot = urinfo %>% group_by(agestrain,expage,expstrain) %>% summarize(nrepaired = sum(yfpclass =='turnedon' & (lastofftime < 5.7 + 14*2) ),total = n())
toplot = mutate(toplot,fracrepaired  = nrepaired/total)
fraclabels = paste(toplot$nrepaired,'/',toplot$total,sep="")
toplot = data.frame(fraclabels,toplot)

p1 = ggplot(toplot,aes(expage,fracrepaired)) + geom_bar(stat="identity",width=bwidth,fill=fillc) + facet_grid(.~expstrain) + border  + ylabel + yscale + themes + geom_text(aes(label=fraclabels,y=fracrepaired+0.05),size=6)
plotfilename = paste(outfolder,'fractionofcellsoffat9hpdx_repairedinnext5h.pdf',sep="") 
ggsave(file = plotfilename,p1,dpi=600,width=9,height=4.5)
removeggplottitles(p1,plotfilename,4,600,9,4.5)


#Now we consider only the cells that are alive at the end of the 5 hour window.  Including burst cells might not make sense, since we would be including cells that are dying by other means. 
urinfo = filter(info,yfpclass9hpdx == 'alwaysoff')
urinfo = filter(urinfo,(lastobservationtime > (doxtime + 14*6)))

toplot = urinfo %>% group_by(agestrain,expage,expstrain) %>% summarize(nrepaired = sum(yfpclass =='turnedon' & (lastofftime < 5.7 + 14*2) ),total = n())
toplot = mutate(toplot,fracrepaired  = nrepaired/total)
fraclabels = paste(toplot$nrepaired,'/',toplot$total,sep="")
toplot = data.frame(fraclabels,toplot)

p1 = ggplot(toplot,aes(expage,fracrepaired)) + geom_bar(stat="identity",width=bwidth,fill=fillc) + facet_grid(.~expstrain) + border  + ylabel + yscale + themes + geom_text(aes(label=fraclabels,y=fracrepaired+0.05),size=6)
plotfilename = paste(outfolder,'fractionofcellsoffat9hpdx_repaired_andalive_innext5h.pdf',sep="") 
ggsave(file = plotfilename,p1,dpi=600,width=9,height=4.5)
removeggplottitles(p1,plotfilename,4,600,9,4.5)

