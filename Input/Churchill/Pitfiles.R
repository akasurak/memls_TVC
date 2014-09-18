#AKasurak
#File to mangle pitfiles to appropriate layering etc.



########################################
########################################
## Libraries
########################################
########################################
library(lattice)


########################################
########################################
## Options
########################################
########################################
options(stringsAsFactors=F)
lattice.options(default.args = list(as.table = TRUE))
setwd("~/Documents/School/Projects/MEMLS/MEMLS3&a/Input/Churchill")



########################################
########################################
## Sourcefiles
########################################
########################################


source('~/Documents/School/Rfunctions/R functions env.R')
source('~/Documents/School/Projects/Paper-PeatModel/AGU2012/mironovSims/JRT_functions.R')

#Have run the following:
source('~/Documents/School/Projects/Visual-Scatterometry/pitfiles.R')
#To produce:
load('~/Documents/School/Projects/Paper-PeatModel/AGU2012/mironovSims/pitfiles/pitfiles_tundra.rdata')

#Run inputbuilder portion of source('~/Documents/School/Projects/Paper-PeatModel/AGU2012/mironovSims/AGU2012-MironovSims.R')
#This builds the following files:
allsweobs=read.csv("~/Documents/School/Projects/Paper-PeatModel/AGU2012/mironovSims/pitfiles/allsweobs.csv")

########################################
########################################
## Build MEMLS3a inputfiles
########################################
########################################

allsweobs=data.frame(date="s",Date=strptime("19990101-0000","%Y%m%d-%H%M"),sweA=1,sweB=1,sweE=1,sweT=1,tgnd=1,tair=1,tsnow=1,depthpit=1,depthmagna=1,pdsw=1,DHfrac=1,meanGS_dobs=1,meanGS_dobs_upper=1,meanGS_dobs_lower=1,sdGS_dobs_upper=1,sdGS_dobs_lower=1,file="")

for(date in sort(unique(unlist(lapply(pitfiles,function(x)paste(format(x$info$datetime,"%Y-%j"),gsub("(.*?)_([ABS])([12]_.*)","-\\2",x$info$file),"-",toupper(substr(x$info$landcover,1,1)),ifelse(x$info$isDestructive,"D","S"),sep="")))))   ){#){
	pit=pitfiles[[unlist(lapply(pitfiles,function(x)if(paste(format(x$info$datetime,"%Y-%j"),gsub("(.*?)_([ABS])([12]_.*)","-\\2",x$info$file),"-",toupper(substr(x$info$landcover,1,1)),ifelse(x$info$isDestructive,"D","S"),sep="")==date)return(x$info$file)))[1]]]
	cat("\n\n\n",date)
	print(pit$info$file)
	delim=delimitLayers(strat=pit$stratigraphy,nlayers=2,splitStep=NA)
	dindex_1=which(pit$stratigraphy$profile$top>delim)
	dindex_2=which(pit$stratigraphy$profile$top<=delim)
				
	allsweobs[date,]=data.frame(
		date=pit$info$date,
		Date=pit$info$datetime,
		sweA=mean(pit$density$profile$densityA,na.rm=T)*pit$density$info$totaldepth*10,
		sweB=mean(pit$density$profile$densityB,na.rm=T)*pit$density$info$totaldepth*10,
		sweE=mean(pit$esc30$samples$density,na.rm=T)*mean(pit$esc30$samples$depth,na.rm=T)*10,
		sweT=mean(c(pit$density$profile$densityA,pit$density$profile$densityB,pit$esc30$samples$density),na.rm=T)*mean(c(pit$density$info$totaldepth,pit$esc30$samples$depth,pit$magna$meansd),na.rm=T)*10,
		tgnd=pit$temperature$scene$soil,
		tair=pit$temperature$scene$air,
		tsnow=mean(pit$temperature$profile$temperature,na.rm=T),
		depthpit=pit$density$info$totaldepth,
		depthmagna=pit$magna$meansd,
		pdsw=mean(c(pit$density$profile$densityA,pit$density$profile$densityB,pit$esc30$samples$density),na.rm=T),
		DHfrac=sum(pit$stratigraphy$profile$top[pit$stratigraphy$profile$grainType=="Depth Hoar"]-pit$stratigraphy$profile$bottom[pit$stratigraphy$profile$grainType=="Depth Hoar"],na.rm=T)/sum(pit$stratigraphy$profile$top-pit$stratigraphy$profile$bottom,na.rm=T),
		meanGS_dobs=pit$pitfile$OneLYR$do_est,#mean(,na.rm=T)
		meanGS_dobs_upper=pit$pitfile$TwoLYR$do_est[1],
		meanGS_dobs_lower=pit$pitfile$TwoLYR$do_est[2],
		sdGS_dobs_upper=sd(unlist(pit$stratigraphy$profile[dindex_1,3:8]),na.rm=T),
		sdGS_dobs_lower=sd(unlist(pit$stratigraphy$profile[dindex_2,3:8]),na.rm=T),
		file=pit$info$file
		)
}
(allsweobs=allsweobs[-1,])
allsweobs["2010-308-S-TS",]=data.frame( #There is no pitfile as there is no snow!
		date="04/11/2010",
		Date=strptime("04/11/2010-1725","%d/%m/%Y-%H%M"),
		sweA=0,
		sweB=0,
		sweE=0,
		sweT=0,
		tgnd=-1.899,
		tair=-5.414,
		tsnow=-5.414,
		depthpit=0,
		depthmagna=0,
		pdsw=0,
		DHfrac=0,
		meanGS_dobs=0,
		meanGS_dobs_upper=0,
		meanGS_dobs_lower=0,
		sdGS_dobs_upper=0,
		sdGS_dobs_lower=0,
		file="Pit_Tundra_STA_30810_S1_2010_11_04.xls"
		)
allsweobs["2011-056-B-TD",]=data.frame( #There is no pitfile as there is no snow!
		date="25/02/2011",
		Date=strptime("25/02/2011-1600","%d/%m/%Y-%H%M"),
		sweA=0,
		sweB=0,
		sweE=0,
		sweT=0,
		tgnd=-1.899,
		tair=-5.414,
		tsnow=-5.414,
		depthpit=0,
		depthmagna=0,
		pdsw=0,
		DHfrac=0,
		meanGS_dobs=0,
		meanGS_dobs_upper=0,
		meanGS_dobs_lower=0,
		sdGS_dobs_upper=0,
		sdGS_dobs_lower=0,
		file="Pit_Tundra_DST_05611_B1_2011_02_25.xls"
		)
allsweobs["2011-057-B-TD",]=data.frame( #There is no pitfile as there is no snow!
		date="03/03/2011",
		Date=strptime("03/03//2011-1330","%d/%m/%Y-%H%M"),
		sweA=0,
		sweB=0,
		sweE=0,
		sweT=0,
		tgnd=-1.899,
		tair=-5.414,
		tsnow=-5.414,
		depthpit=0,
		depthmagna=0,
		pdsw=0,
		DHfrac=0,
		meanGS_dobs=0,
		meanGS_dobs_upper=0,
		meanGS_dobs_lower=0,
		sdGS_dobs_upper=0,
		sdGS_dobs_lower=0,
		file="Pit_Tundra_DST_05711_B1_2011_02_26.xls"
		)
allsweobs$nor=normalize(allsweobs$sweT,datarange=c(0,50))
#allsweobs$Date=strptime(gsub("(.*?)(201[01]_\\([[:digit:]]{2})\\){}(.*?[.]xls)",allsweobs$date),"%d/%m/%Y")

allsweobs=allsweobs[order(allsweobs$Date),]

write.csv(allsweobs,file="allsweobs.csv")

######################################
##	Adjust paramters for MEMLS
allsweobs$snowWet=0
allsweobs$salinity=0
allsweobs$corr_exp=0.5*allsweobs$meanGS_dobs*(1-(allsweobs$pdsw*1000/917.00))

#subset to working site
(sweobs=allsweobs[rownames(allsweobs) %in% c("2010-308-S-TS","2010-319-S-TS", "2010-327-S-TS", "2010-340-S-TS", "2010-345-S-TS", "2010-352-S-TS", "2011-004-S-TS", "2011-016-S-TS", "2011-023-S-TS", "2011-033-S-TS", "2011-044-S-TS", "2011-047-S-TS", "2011-057-S-TS", "2011-056-A-TD", "2011-056-B-TD", "2011-062-A-TD", "2011-062-A-TD"),])


######################################
##OneLayer
for(i in 1:nrow(sweobs)){
	pit=sweobs[i,]
	print(pit)
	filename=paste(c(row.names(pit),"_L1",".txt"),collapse="")
	#Layer#, layerTemp (K), Vol LWC (0-1), PDSW (kg/m3), thickness (cm), salinity (ppthousand), exp corr len (mm) 
	onelayer=c(1,pit$tsnow+273.13,pit$snowWet,pit$pdsw*1000,pit$depthpit,pit$salinity,pit$corr_exp)
	print(onelayer)
	write.table(t(onelayer),file=filename,sep=" ",row.names=F,col.names=F)
}



