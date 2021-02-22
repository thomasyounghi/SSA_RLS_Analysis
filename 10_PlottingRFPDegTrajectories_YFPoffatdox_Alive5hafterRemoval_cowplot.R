#Her we plot individual RFP trajectories and label by 
setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis/')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/Preprocessing_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/timeseries_func.Rd')
source('/Users/thomasyoung/Dropbox/templates/R_aging_template/functions/func.Rd')
library(dplyr)
library(cowplot)
library(reshape2)
library(scales)

info = read.csv('./CombinedData/info_offatdox_alive5hafter.csv')
rfp = read.csv('./CombinedData/rfpcells_offatdox_alive5hafter.csv')
yfp = read.csv('./CombinedData/yfpcells_offatdox_alive5hafter.csv')
agestrain = info$agestrain
yfp = data.frame(agestrain,yfp)
rfp = data.frame(agestrain,rfp)
rfpcos = read.csv('./figures/rfpdeg_trajectories_yfpalwayson/pooledrfpdeg_meansd.csv');
rfpco = rfpcos[rfpcos$expage=='all',]$qtl95;

outfolder = './figures/fltraj_rfpdeg_cowplot/'
outfolder1 = './figures/fltraj_rfpdeg_cowplotlarger/'


#figure settings:
ylogscale = scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x)),limits=c(0.2,2000))
xscale = scale_x_continuous(breaks = seq(0,12,2))
fonts = theme(axis.text=element_text(size=8), axis.title=element_text(size=8),strip.text.x = element_text(size = 8),plot.title=element_text(size=8))
ylims = ylim(-10,1000)
xlims = xlim(0,12.5)
doxwindow = data.frame(14/6,38/6);
colnames(doxwindow) = c('start','end')
doxrect = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),color = "transparent",fill="blue",alpha=0.3)
doxrect1 = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),color = "transparent",fill="green",alpha=0.3)
doxrect2 = geom_rect(data=doxwindow,inherit.aes=FALSE,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),color = "transparent",fill='green3',alpha=0.3)
border = theme(panel.border=element_rect(colour="black",linetype="solid",fill=NA,size=1))
linewidth = 0.2
pointsize = 0.2



doxlabel = annotate("text",label ="dox treatment",x = 4.25,y=1500,colour='red',size=5)
ylabel = ylab('YFP(AU)')
xlabel = xlab('time(h)')
removelegend = theme(legend.position="none")


#Only look at cassettes containing the mCherry-degron
condition = (info$strain=='yTY146') | (info$strain =='yTY147')
info = filter(info,condition)
rfp = filter(rfp,condition)
yfp = filter(yfp,condition)



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
meltedyfp = melt(yfp,id.vars=c('agestrain','id','xy','trap','date','replabel'))
time = as.numeric(meltedyfp$variable)
time = (time-1)/2
meltedyfp = data.frame(time,meltedyfp)

meltedrfp = melt(rfp,id.vars=c('agestrain','id','xy','trap','date','replabel'))
time = as.numeric(meltedrfp$variable)
time = (time-1)/2
meltedrfp = data.frame(time,meltedrfp)


#Now plotting the results by experiment with hours on the x-axis
yfpbyagestrain = split(meltedyfp,meltedyfp$agestrain)
rfpbyagestrain = split(meltedrfp,meltedrfp$agestrain)

#Only look at the non-emply elements of the list
yfpbyagestrain = yfpbyagestrain[lapply(yfpbyagestrain,nrow)!=0]
rfpbyagestrain = rfpbyagestrain[lapply(rfpbyagestrain,nrow)!=0]


for(i in 1:8){
	curryfp = yfpbyagestrain[[i]]
	curryfpsplit = split(curryfp,curryfp$id)	
	listofplots = list();
	for(j in 1:16){
		curryfp = curryfpsplit[[j]]
		title = paste("Cell #",curryfp$id,sep = "")
		p1 <- ggplot(curryfp,aes(time,value))+geom_line() + ylim(-20,1600)+fonts + xlim(0,13) + ggtitle(title) + ylab('YFP(AU)') + xlab('time(h)') + doxrect
		listofplots[[j]] = p1;
	}
	p1 = plot_grid(plotlist = listofplots,nrow=4,ncol=4)
	plotfilename = paste(outfolder,'yfp_',curryfp[1,]$agestrain,'.pdf',sep="")
	ggsave(file = plotfilename,p1,dpi=600,width=6,height=6)
	
	currrfp = rfpbyagestrain[[i]]
	currrfpsplit = split(currrfp,currrfp$id)	
	listofplots = list();
	for(j in 1:16){
		currrfp = currrfpsplit[[j]]
		currrfp = filter(currrfp,time<12.5)
		title = paste("Cell #",currrfp$id,sep = "")
		p1 <- ggplot(currrfp,aes(time,value))+geom_line() + ylim(-10,80)+fonts + xlim(0,13) + ggtitle(title) + ylab('RFP(AU)') + xlab('time(h)') + doxrect + geom_hline(yintercept=rfpco,colour="red")
		listofplots[[j]] = p1;
	}
	p1 = plot_grid(plotlist = listofplots,nrow=4,ncol=4)
	plotfilename = paste(outfolder,'rfp_',currrfp[1,]$agestrain,'.pdf',sep="")
	ggsave(file = plotfilename,p1,dpi=600,width=6,height=6)
	
	
}




##Plotting Individual Trajectories
fonts1 = theme(axis.text=element_text(size=18), axis.title=element_text(size=20),strip.text.x = element_text(size = 20),plot.title=element_text(size=20))

for(i in 1:3){

	currrfp = rfpbyagestrain[[i]]
	currrfpsplit = split(currrfp,currrfp$id)	
	listofplots = list();
	
	currrfp = currrfpsplit[[2]]
	currrfp = filter(currrfp,time<=11.5)
	print(currrfp)
	title = paste("Cell #",currrfp$id,sep = "")
	p1 <- ggplot(currrfp,aes(time,value))+geom_line() + geom_point() + ylim(0,80)+fonts1 + xlim(0,12) + ylab('RFP (a.u.)') + xlab('time(h)') + doxrect2 + geom_hline(yintercept=rfpco,colour="red")+xscale+border
	p1

	plotfilename = paste(outfolder1,'rfp_',currrfp[1,]$agestrain,'.pdf',sep="")
	print(plotfilename)
	ggsave(file = plotfilename,p1,dpi=600,width=6,height=4.5)
	removeggplottitles(p1,plotfilename,4,600,6,4.5)
	
}

