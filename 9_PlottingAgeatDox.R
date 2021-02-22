setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(cowplot)
library(reshape2)
library(grid)


#original ebwidth is 0.25
#figure settings:
themes = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
themes = theme(axis.text=element_text(size=17), axis.title=element_text(size=20),strip.text.x = element_text(size = 20))
ylims = ylim(-0,25)
jitterw = 0.12
ebwidth = 0.35;
ylabel = ylab("Age(generations)");
xscale = scale_x_discrete(expand = c(0.6,0))
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))



#folder to save in
outfolder = './figures/ageatdox/'


info = read.csv('./CombinedData/info_offatdox_alive5hafter.csv')

info$agestrain = factor(info$agestrain,levels = c('SSAcontrol young','SSAcontrol old','SSA young','SSA old','SSAhet young','SSAhet old','SSA Dnl4ko young','SSA Dnl4ko old','SSAdegcontrol young','SSAdegcontrol old','SSAdeg young','SSAdeg old','SSA3xCln2 old'))
info$expage = factor(info$expage,levels = c('young','old'))


#Plotting mean+/-sd of the age at dox by experiment  
stats <- info %>% group_by(strain,doxtime,replabel,agestrain,expage) %>% summarize(meanageatdox = mean(ageatdox),total = n())
tograph <- stats %>% group_by(strain,doxtime,agestrain,expage) %>% summarize(meandoxage = mean(meanageatdox),sem = sd(meanageatdox)/sqrt(n()))
write.csv(tograph,paste(outfolder,'ageatdox_meanandsemofreplicateaverages.csv',sep=""))



#Graph for the SSA strain
p1 = ggplot(tograph[tograph$strain =='yTY125',],aes(expage,meandoxage)) + geom_point() + geom_errorbar(aes(ymin = meandoxage-sem,ymax = meandoxage+sem),width=ebwidth) + ylims +xlab('age group') + ggtitle('') + themes + ylabel
plotfilename = paste(outfolder,'meanagebyexp_SSA.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)

##SSAhet graph
p1 = ggplot(tograph[tograph$strain =='yTY133',],aes(expage,meandoxage)) + geom_point() + geom_errorbar(aes(ymin = meandoxage-sem,ymax = meandoxage+sem),width=ebwidth) + ylim(0,25)  +xlab('age group') + ylims + ggtitle('')+themes + ylabel
plotfilename = paste(outfolder,'meanagebyexp_SSAhet.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)



##SSA Dnl4ko graph
p1 = ggplot(tograph[tograph$strain =='yTY149',],aes(expage,meandoxage)) + geom_point() + geom_errorbar(aes(ymin = meandoxage-sem,ymax = meandoxage+sem),width=ebwidth) + ylims  +xlab('age group') + ylab('age(generations)') + ggtitle('')+themes + ylabel
plotfilename = paste(outfolder,'meanagebyexp_SSA_Dnl4ko.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)



##SSA control graph
p1 = ggplot(tograph[tograph$strain =='yTY126',],aes(expage,meandoxage)) + geom_point() + geom_errorbar(aes(ymin = meandoxage-sem,ymax = meandoxage+sem),width=ebwidth) + ylims  +xlab('age group')  + ggtitle('')+themes + ylabel
plotfilename = paste(outfolder,'meanagebyexp_SSAcontrol.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)

#Graph for the SSA3xCln2 strain
p1 = ggplot(tograph[tograph$strain =='yTY161a',],aes(expage,meandoxage)) + geom_point() + geom_errorbar(aes(ymin = meandoxage-sem,ymax = meandoxage+sem),width=ebwidth) + ylims +xlab('age group') + ggtitle('') + themes + ylabel
plotfilename = paste(outfolder,'meanagebyexp_SSA3xCln2.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)


###
#Now pooling across replicates to get mean and sd
tograph <- info %>% group_by(strain,doxtime,agestrain) %>% summarize(meanageatdox = mean(ageatdox),sd = sd(ageatdox))
p1 = ggplot(tograph[tograph$strain =='yTY125',],aes(agestrain,meanageatdox)) + geom_point() + geom_errorbar(aes(ymin = meanageatdox-sd,ymax = meanageatdox+sd),width=0.5) + ylim(0,25)  +xlab('age group') + ylab('age') + ggtitle('')+theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),strip.text.x = element_text(size = 12
))
plotfilename = paste(outfolder,'meanagepooled_SSA.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)

p1 = ggplot(tograph[tograph$strain =='yTY133',],aes(agestrain,meanageatdox)) + geom_point() + geom_errorbar(aes(ymin = meanageatdox-sd,ymax = meanageatdox+sd),width=0.5) + ylim(0,25)  +xlab('') + ylab('age') + ggtitle('')+theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),strip.text.x = element_text(size = 12
))
plotfilename = paste(outfolder,'meanagepooled_SSAhet.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)

p1 = ggplot(tograph[tograph$strain =='yTY149',],aes(agestrain,meanageatdox)) + geom_point() + geom_errorbar(aes(ymin = meanageatdox-sd,ymax = meanageatdox+sd),width=0.5) + ylim(0,25)  +xlab('') + ylab('age') + ggtitle('')+theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),strip.text.x = element_text(size = 12
))
plotfilename = paste(outfolder,'meanagepooled_SSA_Dnl4ko.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)

p1 = ggplot(tograph[tograph$strain =='yTY126',],aes(agestrain,meanageatdox)) + geom_point() + geom_errorbar(aes(ymin = meanageatdox-sd,ymax = meanageatdox+sd),width=0.5) + ylim(0,25)  +xlab('') + ylab('age') + ggtitle('')+theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15),strip.text.x = element_text(size = 12
))
plotfilename = paste(outfolder,'meanagepooled_SSAcontrol.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)






#Jitter plots of cell ages at dox with replicate based mean +/-sem overlaid
#Now plotting the raw age at dox data pooled across replicates sharing the same strain and age
stats <- info %>% group_by(strain,doxtime,replabel,agestrain,expage,expstrain) %>% summarize(meanageatdox = mean(ageatdox),total = n())
tograph <- stats %>% group_by(strain,doxtime,agestrain,expage,expstrain) %>% summarize(ageatdox= mean(meanageatdox),sem = sd(meanageatdox)/sqrt(n()))
outtable = paste(outfolder,'mean+-sem_ageatdox.csv',sep="")
write.csv(tograph,outtable,row.names=FALSE)


#SSA
xlabel = xlab('')
tograph1= tograph[tograph$strain =='yTY125',]
p1 = ggplot(info[info$strain == 'yTY125',],aes(expage,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ ylims +xlabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expage,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth) + ylabel + xscale + border
plotfilename = paste(outfolder,'rawageatdox_meansem_SSA.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


tograph1= tograph[tograph$strain =='yTY126',]
p1 = ggplot(info[info$strain == 'yTY126',],aes(expage,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ ylims +xlabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expage,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth)+ ylabel + xscale
plotfilename = paste(outfolder,'rawageatdox_meansem_SSAcontrol.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


tograph1= tograph[tograph$strain =='yTY133',]
p1 = ggplot(info[info$strain == 'yTY133',],aes(expage,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ ylims +xlabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expage,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth) + ylabel
plotfilename = paste(outfolder,'rawageatdox_meansem_SSAhet.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)


tograph1= tograph[tograph$strain =='yTY149',]
p1 = ggplot(info[info$strain == 'yTY149',],aes(expage,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ ylims +xlabel + ylabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expage,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth)
plotfilename = paste(outfolder,'rawageatdox_meansem_SSA_Dnl4ko.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)

tograph1= tograph[tograph$strain =='yTY147',]
p1 = ggplot(info[info$strain == 'yTY147',],aes(expage,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ ylims +xlabel + ylabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expage,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth)
plotfilename = paste(outfolder,'rawageatdox_meansem_SSAdeg.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)

#Old SSAcontrol side by side with old SSAdegcontrol. For paper figure 4f.
tograph1= tograph[(tograph$strain =='yTY126'|tograph$strain=='yTY146')&(tograph$expage=='old'),]
p1 = ggplot(info[(info$strain == 'yTY126' | info$strain == 'yTY146') & (info$expage=='old'),],aes(expstrain,ageatdox)) + geom_jitter(width=jitterw*1.25,height=0,colour='red',size=1)+ ylims +xlabel + ylabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expstrain,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth*1.25)+border
plotfilename = paste(outfolder,'rawageatdox_meansem_SSAcontrol_degornodeg_old.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=4,height=4.5)
removeggplottext(p1,plotfilename,4,600,4,4.5)


#SSA3xCln2
xlabel = xlab('')
tograph1= tograph[tograph$strain =='yTY161a',]
p1 = ggplot(info[info$strain == 'yTY161a',],aes(expage,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ ylims +xlabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expage,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth) + ylabel + xscale + border
plotfilename = paste(outfolder,'rawageatdox_meansem_SSA3xCln2.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)
plotfilename = plotfilename = paste(outfolder,'rawageatdox_meansem_SSA3xCln2_narrower.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#SSA3 Rad52 and Rad52 KO
xlabel = xlab('')
tograph1= tograph[tograph$strain =='yTY164a' | tograph$strain=='yTY165a',]
p1 = ggplot(info[info$strain == 'yTY164a' | info$strain == 'yTY165a',],aes(expstrain,ageatdox)) + geom_jitter(width=jitterw,height=0,colour='red',size=1)+ ylims +xlabel + ggtitle('')+themes+geom_point(data = tograph1,aes(expstrain,ageatdox)) + geom_errorbar(data = tograph1,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth) + ylabel + xscale + border
plotfilename = paste(outfolder,'rawageatdox_meansem_SSA_Rad52_51ko.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=5,height=4.5)
removeggplottext(p1,plotfilename,4,600,5,4.5)
plotfilename = plotfilename = paste(outfolder,'SSARad52_Rad51ko_ageatdox_narrower.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=3.5,height=4.5)
removeggplottext(p1,plotfilename,4,600,3.5,4.5)


#A plot for all strains side by side
info$expstrain = factor(info$expstrain,levels = c('SSA','SSAcontrol','SSAdeg','SSAdegcontrol','SSA Dnl4ko','SSA3xCln2','SSAhet','SSA Rad51ko','SSA Rad52ko'))
tograph$expstrain = factor(tograph$expstrain,levels = c('SSA','SSAcontrol','SSAdeg','SSAdegcontrol','SSA Dnl4ko','SSA3xCln2','SSAhet','SSA Rad51ko','SSA Rad52ko'))
p1 = ggplot(info,aes(expage,ageatdox))+ geom_jitter(width=jitterw,height=0,colour='red',size=1) + ylims + xlabel+ ylabel + ggtitle('') + themes + geom_point(data = tograph,aes(expage,ageatdox)) + geom_errorbar(data = tograph,aes(ymin = ageatdox-sem,ymax = ageatdox+sem),width=ebwidth) + facet_grid(.~expstrain,scales="free") +border
plotfilename = paste(outfolder,'rawageatdox_meansem_allstrains.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=8,height=4.5)
removeggplottext(p1,plotfilename,4,600,8,4.5)


geom_text(data=tograph,aes(x=expage,y=ageatdox,label=floor(ageatdox)),nudge_x=0.1)







#Plotting individual data points, age at dox for each strain separately

#SSA control
p1 = ggplot(info[info$strain == 'yTY126',],aes(expage,ageatdox)) + geom_jitter()+ ylims+xlab('age group') + ylab('age') + ggtitle('')+themes
plotfilename = paste(outfolder,'rawageatdox_SSAcontrol.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)


#SSA heterology
p1 = ggplot(info[info$strain == 'yTY133',],aes(expage,ageatdox)) + geom_jitter()+ ylims +xlab('age group') + ylab('age') + ggtitle('')+themes
plotfilename = paste(outfolder,'rawageatdox_SSAhet.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)


#SSA Dnl4ko
p1 = ggplot(info[info$strain == 'yTY149',],aes(expage,ageatdox)) + geom_jitter()+ ylims+xlab('age group') + ylab('age') + ggtitle('')+themes
plotfilename = paste(outfolder,'rawageatdox_SSA_Dnl4ko.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)



