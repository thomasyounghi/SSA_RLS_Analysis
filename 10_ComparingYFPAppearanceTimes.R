#Comparing the time taken for YFP to appear across all experiments.
#The mean and sem for replicate averaged YFP appearance times are plotted.
#Results are saved in './figures/YFPappearancetimes/'

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
#library(plotmath)


#original ebwidth is 0.25
#figure settings:
theme_set(theme_cowplot())
themes = theme(axis.text=element_text(size=18), axis.title=element_text(size=16),strip.text.x = element_text(size = 16))
ylims = ylim(-0,25)
jitterw = 0.25
ebwidth = 0.35;
ylabel = ylab("time(h)");
xscale = scale_x_discrete(expand = c(0.6,0))
yscale = scale_y_continuous(expand = c(0,0),limits = c(0,16.3))
fillc = "gray66"
pointalpha = 0.7;
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))


strainnames = list('SSA' = 'SSA','SSA Dnl4ko' = expression(paste('SSA ',italic('dnl4'),Delta)),'SSAhet'='SSAhet','SSA3xCln2' = 'SSA3xCln2')

strain_labeller <- function(variable,value){
	return(strainnames[value])
}


#folder to save in
outfolder = './figures/YFPappearancetimes/'


info = read.csv('./CombinedData/info_offatdox_alive5hafter.csv')
info = filter(info,expstrain!='SSAcontrol' & expstrain!='SSAdegcontrol' & expstrain != 'SSA Rad51ko' & expstrain != 'SSA Rad52ko' & !(expstrain =='SSA3xCln2'&expage=='young'))
info$expstrain2 = factor(info$expstrain,levels = c('SSAcontrol','SSA','SSA Dnl4ko','SSAdeg','SSAdegcontrol','SSA3xCln2','SSAhet'),labels = c("SSAcontrol",'SSA',expression(paste('SSA ',italic('dnl4'),Delta)),'SSA+RFPdegron','SSAcontrol+RFPdegron','SSA3xCln2',"SSAheterology"))
info$expage = factor(info$expage,levels=c('young','old'))


#Among the cells in which YFP appears, is there any difference in when it appears between young and old cells
#If it takes longer to appear in old cells, that would suggest that delay in cutting, repair or YFP expression could be a reason why a lower % of cells are detected as repaired.
condition = info$yfpclass == 'turnedon';
info = filter(info,condition)
firstontimeh = info$lastofftime/2;
info = data.frame(firstontimeh,info)


summary <- info %>% group_by(agestrain,expage,expstrain2,replabel) %>% summarize(meanfirstonth = mean(firstontimeh),sdfirstonth = sd(firstontimeh),total = n())
meanandsem = summary %>% group_by(agestrain,expage,expstrain2) %>% summarize(firstontimeh = mean(meanfirstonth),sem = sd(meanfirstonth)/sqrt(n()))



p1 = ggplot(info,aes(expage,firstontimeh))+geom_jitter(alpha=pointalpha,height=0,width = jitterw,aes(colour=factor(replabel)))+facet_grid(.~expstrain2,labeller=label_parsed,scales="free")+ylim(0,40) + ylabel + xlab('age') + geom_hline(yintercept=140/60) + geom_errorbar(data = meanandsem,aes(ymin = firstontimeh-sem,ymax=firstontimeh+sem),width=ebwidth)+ geom_point(data=meanandsem,aes(expage,firstontimeh))+labs(colour = 'Replicate') + themes + yscale + border + guides(colour=FALSE)
plotfilename = paste(outfolder,'yfpappearancetimes_all.pdf',sep="") 
ggsave(file = plotfilename,p1,dpi=600,width=9,height=4.5)
removeggplottext(p1,plotfilename,4,600,9,4.5)



p1 = ggplot(info,aes(expage,lastobservationtime))+geom_boxplot(width=ebwidth)+geom_jitter(height=0,width = jitterw,aes(colour=factor(replabel)))+facet_grid(.~expstrain2,labeller=label_parsed)+geom_hline(aes(yintercept=doxtime+14*6)) + ylabel + xlab('age') + geom_hline(yintercept=140/60) + geom_point(data=meanandsem,aes(expage,firstontimeh))+labs(colour = 'Burst?') + themes + border
plotfilename = paste(outfolder,'lastobservationtimes_all.pdf',sep="") 
ggsave(file = plotfilename,p1,dpi=600,width=9,height=4.5)
removeggplottext(p1,plotfilename,4,600,9,4.5)

