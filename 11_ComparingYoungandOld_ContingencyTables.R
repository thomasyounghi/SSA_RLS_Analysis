#Compute repair fractions (based on YFP trajectory) per age+strain by combining cells over replicates with the same age and strain. Plot these age+strain repair fractions
#Use counts of unrepaired and repaired cells for each age+strain combination to generate contingency tables. Rows correspond to age (young or old experiment), columns correspond to repair outcome (repaired or unrepaired).
#Perform Fisher's Exact Test for each strain's contingency table to assess the likelihood of a contingency table at least as extreme as that observed given no association between repair outcome and age.
#Results saved in './ContingencyTables/'

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(dplyr)
library(cowplot)
library(reshape2)
library(gridExtra)
library(grid)


#figure settings:
themes = theme(axis.text=element_text(size=11), axis.title=element_text(size=16),strip.text.x = element_text(size = 8,margin = margin(.3,0,0.3,0,"cm")),strip.background=element_rect(fill="gray66"))
ylims = ylim(0,25)
ebwidth = 0.1;
jitterw = 0.15;
bwidth = 0.4
cartcoordy= coord_cartesian(ylim = c(0,830),clip='off')
yscale = scale_y_continuous(expand = c(0,0),limits=c(0,1))
removelegend = theme(legend.position="none")
fillc = "gray66"
pointalpha = 0.6;


#Pool replicates if they share the same strain and dox time
#Contingency table by strain for YFP repair out of the total
outfolder = './ContingencyTables/'
yfpbyexp = read.csv('./Tables/yfprepair_outoftotal.csv')


#Making nicer labels for the experimental conditions
expage = factor(yfpbyexp$doxtime)
levels(expage) = c('young','old')
strain = factor(yfpbyexp$strain,levels = c('yTY126','yTY125','yTY133','yTY149','yTY146','yTY147','yTY161a','yTY164a','yTY165a'))
expstrain = strain
levels(expstrain) = c('SSAcontrol','SSA','SSAhet','SSA Dnl4ko','SSAdegcontrol','SSAdeg','SSA3xCln2','SSA Rad52ko','SSA Rad51ko')
agestrain = factor(paste(expstrain,expage),levels = c('SSAcontrol young','SSAcontrol old','SSA young','SSA old','SSAhet young','SSAhet old','SSA Dnl4ko young','SSA Dnl4ko old','SSAdegcontrol young','SSAdegcontrol old','SSAdeg young','SSAdeg old','SSA3xCln2 young','SSA3xCln2 old','SSA Rad52ko young','SSA Rad51ko young'))

#Making the contingency tables
yfpbyexp = data.frame(agestrain,expstrain,expage,yfpbyexp)

yfpbystrain <- yfpbyexp %>% group_by(agestrain,strain,doxtime,expstrain,expage) %>% summarize(repaired = sum(nrepair),total = sum(total))
yfpbystrain <- mutate(yfpbystrain,notrepaired = total-repaired)
yfpbystrain = data.frame(yfpbystrain,row.names=yfpbystrain$agestrain)

repairtabbystrain = split(yfpbystrain[,c(6,8)],yfpbystrain$expstrain)

#Saving contingency tables
for(i in 1:length(repairtabbystrain)){
	strain = names(repairtabbystrain)[i]
	outfn = paste(outfolder,strain,'_yfprepair.csv',sep='')
	write.csv(repairtabbystrain[i],outfn)
}



##Plotting fraction repaired (pooled across replicates of each age+strain combination)
tograph <- mutate(yfpbystrain,prepaired = repaired/total,SEp = sqrt(prepaired*(1-prepaired)/total))
p1 = ggplot(tograph,aes(expage,prepaired)) + geom_bar(stat = "identity",fill=fillc,width = bwidth) + geom_point()+ geom_errorbar(aes(ymin = prepaired-SEp,ymax=prepaired+SEp),width=0.5) + ylim(0,1)  +xlab('age') + ylab('Repair Efficiency (%)') + ggtitle('')+ themes + facet_grid(.~expstrain) + yscale
plotfilename = paste(outfolder,'allssastrains_pooledfractionrepaired.pdf',sep="")
ggsave(file = plotfilename,p1,dpi=600,width=9,height=4.5)



#Applying Fisher's Exact Test to Each Contingency Table
#Saving the result for each test in a table.
repairtabbystrain = repairtabbystrain[lapply(repairtabbystrain,nrow)==2]
results = lapply(repairtabbystrain,fisher.test)
pvalues <- sapply(results,'[[','p.value')
oddsratconf95 <- sapply(results,'[[','conf.int')
taboffishertests = rbind(pvalues,oddsratconf95)
row.names(taboffishertests)[2:3]=c('oddsratio95confint_lower','oddsratio95confint_higher')
write.csv(taboffishertests,'./ContingencyTables/fisherexacttestresults.csv')


