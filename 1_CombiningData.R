#Make 1 big data frame each for raw YFP data, RFP data, info in the 'Data' folder
#Identifies and labels replicates (experiments sharing the same strain and age of dox induction)
#Add labels to be used for publication (redundant with the official strain names, doxtime, etc)


setwd('/Users/thomasyoung/Dropbox/MovieProcessing/March2018_Analysis_git')
source('./functions/timeseries_func.Rd')
source('./functions/func.Rd')
source('./functions/Preprocessing_func.Rd')
library(dplyr)
library(stringr)

#Open up all the files in data, to get one big file of info, YFP, RFP
infofilestoopen = list.files("./Data/",pattern="[[:alnum:]]*infoclean.csv")
infofilestoopen = paste('./Data/',infofilestoopen,sep="")
yfpfilestoopen = list.files("./Data/",pattern="[[:alnum:]]*yfpmeanroiraw.csv")
yfpfilestoopen = paste('./Data/',yfpfilestoopen,sep="")
rfpfilestoopen = list.files("./Data/",pattern="[[:alnum:]]*rfpmeanroiraw.csv")
rfpfilestoopen = paste('./Data/',rfpfilestoopen,sep="")


info = concatenatefiles(infofilestoopen)
yfp = concatenatefiles(yfpfilestoopen)
rfp = concatenatefiles(rfpfilestoopen)

#Check for gaps in the budding time matrix
bt = getbudtimes(info)
btgaps = apply(bt,1,checkgaps)
btdecreasing = apply(bt,1,checkdecreasing)
btrepeats = apply(bt,1,checkrepeats)
table(btgaps)
table(btdecreasing)
table(btrepeats)
write.csv(info[!btgaps,],'./ProblematicBTs/NAbetweenbts.csv',row.names=FALSE)
write.csv(info[btdecreasing,],'./ProblematicBTs/decreasingbts.csv',row.names=FALSE)
write.csv(info[btrepeats,],'./ProblematicBTs/repeatedbts.csv',row.names=FALSE)


#In situations where there are repeated buddingtime, remove the repeated buddingtimes and shift the later budding times over
condition = btrepeats & !btdecreasing & !btgaps
bttofix = bt[condition,]
fixedbt = t(apply(bttofix,1,removerepeatsinincreasing))
birthcol = which(names(info)=='birth')
info[condition,(birthcol+1):ncol(info)]=fixedbt

#Eliminate all rows for which budding times are decreasing or there are gaps
condition = !(btdecreasing | btgaps)
table(condition)
info = filter(info,condition)
yfp = filter(yfp,condition)
rfp = filter(rfp,condition)


#Plot the birth distributions (boxplot). Also calculate the doxtime, which is important for all experiments.
doxtime = info$lastibeforepart2 + info$firstpart2ipostdoxstart;
info = cbind(doxtime,info);
neworder = order(info$strain,doxtime,info$date,info$lane)
info = info[neworder,]
yfp = yfp[neworder,]
rfp = rfp[neworder,]
yfp = yfp[,-1]
rfp = rfp[,-1]


#Plot the birth distributions (boxplot). Replicates near one another. 
experiment = paste(info$strain,info$doxtime, info$date,info$lane)
experiment = factor(experiment)


#Providing an id to each cell
id = 1:nrow(info)


#identifying replicates. Experiments sharing the same strain and doxycycline treatment time, but derived from different colonies. Order of the labeling specified by (date lane)
replabel = identifyreplicates(info);

#Labels to use for plotting
expage = factor(info$doxtime)
levels(expage) = c('young','old')
expstrain = factor(info$strain,levels = c('yTY126','yTY125','yTY133','yTY149','yTY146','yTY147','yTY161b','yTY161a','yTY164a','yTY165a'))
levels(expstrain) = c('SSAcontrol','SSA','SSAhet','SSA Dnl4ko','SSAdegcontrol','SSAdeg','SSA2xCln2','SSA3xCln2','SSA Rad52ko','SSA Rad51ko')
agestrain = factor(paste(expstrain,expage),levels = c('SSAcontrol young','SSAcontrol old','SSA young','SSA old','SSAhet young','SSAhet old','SSA Dnl4ko young','SSA Dnl4ko old','SSAdegcontrol young','SSAdegcontrol old','SSAdeg young','SSAdeg old','SSA2xCln2 old','SSA3xCln2 young','SSA3xCln2 old','SSA Rad52ko young','SSA Rad51ko young'))

#adding the additional features to the info data frame
#labels that make things easier to understand.  
bt = getbudtimes(info)
ageatdox = ageattimes(bt,info$doxtime)
info = data.frame(ageatdox,info)

info = data.frame(id,agestrain,expage,expstrain,replabel,experiment,info)
date = as.character(info$date)
date = str_match(date,"(\\w+)_")
date = date[,2]
info$date=date;


#Annotating the yfp file with xy and trap information
xy = info$xy
trap = info$trap
yfp = data.frame(id,date,xy,trap,replabel,yfp)
rfp = data.frame(id,date,xy,trap,replabel,rfp)

write.csv(info,"./CombinedData/infoall.csv",row.names=FALSE)
write.csv(yfp,"./CombinedData/yfpall.csv",row.names=FALSE)
write.csv(rfp,"./CombinedData/rfpall.csv",row.names=FALSE)


