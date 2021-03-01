#Using the YFP based repair phenotype, compute old repair as a fraction of young repair for each strain
#For each strain, compute (old repair efficiency)/(young repair efficiency). repair efficiency = (number of repaired cells over replicates)/(total cells over replicates)
#The goal is to normalize for different young cell repair efficiencies to compare age related declines between strains. 
#Results are saved in './figures/oldtoyoung_repairratio/'

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(grid)

outfolder = './figures/oldtoyoung_repairratio/'

#figure settings:
theme_set(theme_cowplot())
themes = theme(axis.text=element_text(size=16), axis.title=element_text(size=20),strip.text.x = element_text(size = 18))
ylims = ylim(-0,25)
ebwidth = 0.25;
fillc = "gray66"
yscale = scale_y_continuous(expand = c(0,0),limits = c(0,1.8),breaks = seq(0,1.8,0.3))
ylabel = labs(y=expression(paste('%',repair[old]," / %",repair[young]))) 
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))


#Summing over replicates of the same age and strain to compute total repaired cells and total cells. Calculating repair efficiency for each age+strain combination
yfpbyexp = read.csv('./Tables/yfprepair_outoftotal.csv')
yfppooled = yfpbyexp %>% group_by(strain,doxtime) %>% summarize(totalrepaired = sum(nrepair),totalcells = sum(total))
yfppooled = mutate(yfppooled,fracrepair = totalrepaired/totalcells)

yfpold = filter(yfppooled,doxtime==169 & strain!= 'yTY126' & strain !='yTY146' & strain != 'yTY164a' & strain != 'yTY165a')
yfpyoung = filter(yfppooled,doxtime==43 & strain!= 'yTY126' & strain !='yTY146' & strain != 'yTY164a' & strain != 'yTY165a')


#Calculating the ratio of old repair efficiency to young repair efficiency for each strain
#Calculate the delta method estimate for the standard error for the log of this ratio
#Calculate a 95% confidence interval for each ratio by assuming the estimate of the 
#log of the ratio is normally distributed
#Save the confidence intervals in a table
selog = sqrt(1/yfpold$totalrepaired + 1/yfpyoung$totalrepaired - 1/yfpold$totalcells - 1/yfpyoung$totalcells);
ratio = yfpold$fracrepair/yfpyoung$fracrepair;
upperfactor = exp(1.96*selog);
lowerfactor = exp(-1.96*selog)
strain = factor(yfpold$strain, levels=c('yTY125','yTY149','yTY147','yTY161a','yTY133'));
tograph = data.frame(strain,ratio,ratio*upperfactor,ratio*lowerfactor)
colnames(tograph) = c('strain','ratio','ci95upper','ci95lower');


#Plotting the ratio of old repair efficiency to young repair efficiency for each strain
#For the paper
strainnames = scale_x_discrete(labels=c('SSA',expression(atop('SSA',paste(italic('dnl4'),Delta))),expression(atop('SSA+','RFPdegron')),expression(atop('SSA','3xCln2')),expression(atop('SSA','heterology'))))

p1 = ggplot(tograph,aes(strain,ratio)) + geom_point() + geom_bar(stat='identity',width = 0.4,fill=fillc) + geom_errorbar(aes(ymin=ci95lower,ymax=ci95upper),width=ebwidth) + xlab('') + ggtitle('') + themes + strainnames + yscale + ylabel + border
plotfilename = paste(outfolder,'oldtoyoungreprat_barwithe95CIdeltamethod.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=10,height=4.5)
removeggplottext(p1,plotfilename,4,600,6,4.5)
write.csv(tograph,paste(outfolder,'ratio_oldrepair_over_youngrepair.csv',sep=""),row.names=FALSE);


#without the border
themes1 = theme(axis.text=element_text(size=14), axis.title=element_text(size=16),strip.text.x = element_text(size = 14))
p1 = ggplot(tograph,aes(strain,ratio)) + geom_point() + geom_bar(stat='identity',width = 0.4,fill=fillc) + geom_errorbar(aes(ymin=ci95lower,ymax=ci95upper),width=ebwidth) + xlab('') + ggtitle('')+themes1 + strainnames+yscale+ylabel 
plotfilename = paste(outfolder,'oldtoyoungreprat_barwithe95CIdeltamethod_noborder.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=6.5,height=4.5)




#Not in the paper. Try estimating the standard error of the old to young repair ratio by resampling.  Given n old cells and m young cells with repair fractions eff_old, and eff_young, repeat many times: sample n old cells from a binomial with probability of success eff_old. Sample m young cells from a binomial with probability of success eff_young. Take the ratio of success in the old sample to success in the young sample.


SSA = filter(yfppooled,strain=='yTY125')
SSAhet = filter(yfppooled,strain=='yTY133')
SSADnl4ko = filter(yfppooled,strain=='yTY149')
SSA3xCln2 = filter(yfppooled,strain=='yTY161a')
SSAdegron = filter(yfppooled,strain=='yTY147')

bsriskratiorepairtab <- function(repairtable,nsamp){
	#order so that young cells are in row1, old cells are in row2.
	repairtable = repairtable[order(repairtable$doxtime),]	
	bssamples = mapply(rbinom,nsamp,repairtable$totalcells,repairtable$fracrepair)
	youngoveroldtotal = repairtable$totalcells[1]/repairtable$totalcells[2]
	oldoveryoungbs = bssamples[,2]/bssamples[,1]*youngoveroldtotal;
	return(oldoveryoungbs)
}

allpooled = list(SSA,SSADnl4ko,SSAdegron,SSA3xCln2,SSAhet)
bsriskratios = lapply(allpooled,bsriskratiorepairtab,10000)
pctils95 = lapply(bsriskratios,quantile,probs = c(0.05,0.95))
#Ranges are similar to the delta method.
#In the future will try sampling with replacement for the old sample and the young sample.
	
	
