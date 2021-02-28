#Here we plot 20 single cell YFP traces for all strains that can be cut and repaired by SSA
#Results saved to './figures/fltraj_forpaper/

setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')

library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(scales)

info = read.csv('./CombinedData/info_offatdox_alive5hafter.csv')
rfp = read.csv('./CombinedData/rfpcells_offatdox_alive5hafter.csv')
yfp = read.csv('./CombinedData/yfpcells_offatdox_alive5hafter.csv')
agestrain = info$agestrain
strain = info$strain
yfp = data.frame(strain,agestrain,yfp)
rfp = data.frame(strain,agestrain,rfp)

outfolder = './figures/fltraj_forpaper/'

#figure settings:
theme_set(theme_cowplot())
ylogscale = scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x)),limits=c(0.2,2000))
xscale = scale_x_continuous(breaks = seq(0,12,2))
fonts = theme(axis.text=element_text(size=8), axis.title=element_text(size=8),strip.text.x = element_text(size = 8),plot.title=element_text(size=8))
ylims = ylim(-10,1000)
xlims = xlim(0,12.5)
doxwindow = data.frame(14/6,38/6);
colnames(doxwindow) = c('start','end')
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=0,ymax=Inf),color = "transparent",fill="green3",alpha=0.3)
doxlabel = annotate("text",label ="dox treatment",x = 4.25,y=1500,colour='red',size=5)
ylabel = ylab('YFP (a.u.)')
xlabel = xlab('time(h)')
removelegend = theme(legend.position="none")
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))
linewidth = 0.2
pointsize = 0.2
hlinewidth = 0.35


#filter out young cells less than 5 generations and old cells < 15 generations at the time of doxycycline exposure
bt = getbudtimes(info)
ageatdox = ageattimes(bt,info$doxtime)
condition = (ageatdox<=5 & info$doxtime==43) | (ageatdox>=15 & info$doxtime==169)
info = filter(info,condition)
rfp = filter(rfp,condition)
yfp = filter(yfp,condition)



#Look at only the cells that are visible for 9 hours after the addition of doxycycline
condition = info$lastobservationtime >= (info$doxtime+ 6*9)



#Plot the rfp trajectories
meltedyfp = melt(yfp,id.vars=c('strain','agestrain','id','xy','trap','date','replabel'))
time = as.numeric(meltedyfp$variable)
time = (time-1)/2
yfpexp = paste(meltedyfp$strain,meltedyfp$date,meltedyfp$replabel)
meltedyfp = data.frame(time,yfpexp,meltedyfp)
meltedyfp = filter(meltedyfp,time<=11.5)

meltedrfp = melt(rfp,id.vars=c('strain','agestrain','id','xy','trap','date','replabel'))
time = as.numeric(meltedrfp$variable)
time = (time-1)/2
rfpexp = paste(meltedrfp$strain,meltedrfp$date,meltedrfp$replabel)
meltedrfp = data.frame(time,rfpexp,meltedrfp)
meltedrfp = filter(meltedrfp,time<=11.5)

#Now plotting the results by experiment with hours on the x-axis
yfpbyagestrain = split(meltedyfp,meltedyfp$agestrain)
rfpbyagestrain = split(meltedrfp,meltedrfp$agestrain)

yfpbyagestrain = yfpbyagestrain[lapply(yfpbyagestrain,nrow)!=0]
rfpbyagestrain = rfpbyagestrain[lapply(rfpbyagestrain,nrow)!=0]



#Repeating the plotting without cell labels, and plotting 20 trajectories for one replicate of each condition.  This is for 
#the paper.  For each strain, we plot 20 cells of one old and one young replicate.  We only do this for the non-control strains
#Only plot YFP.  Don't plot trajectories for the control strains yTY126 and yTY146.
expsummary = read.csv('./figures/countsbyreplicate/repcounts_all.csv') 
exptoplot <- expsummary %>% group_by(strain,expage) %>% summarize(repwithmaxcells = replabel[which.max(numforrepcalc)],date = date[which.max(numforrepcalc)],numcells = max(numforrepcalc))
exptoplot = filter(exptoplot,strain != 'yTY126' & strain != 'yTY146')
exptoplotids = paste(exptoplot$strain,exptoplot$date,exptoplot$repwithmaxcells)


for(i in 1:length(exptoplotids)){
	curryfp = meltedyfp[meltedyfp$yfpexp==exptoplotids[i],]
	curryfpsplit = split(curryfp,curryfp$id)	
	listofplots = list();
	print(length(curryfpsplit))
	for(j in 1:min(20,length(curryfpsplit))){
		print(j)
		curryfp = curryfpsplit[[j]]
		title = ''
		p1 <- ggplot(curryfp,aes(time,value))+geom_line(lwd=linewidth) + geom_point(size = pointsize) + ylim(-20,1600)+fonts + xlim(0,13) + ggtitle(title) + ylab('YFP (a.u.)') + xlab('time(h)') + doxrect + border + geom_hline(yintercept=7,colour="red",size=hlinewidth) + ylogscale ##+ ##ggtitle(curryfp$id)
		listofplots[[j]] = p1;
	}
	p1 = plot_grid(plotlist = listofplots,nrow=5,ncol=4)
	plotfilename = paste(outfolder,'yfp_',curryfp[1,]$agestrain,'_',exptoplotids[i],'.pdf',sep="")
	print(plotfilename)
	ggsave(file = plotfilename,p1,dpi=600,width=6,height=7.5)
	
}
