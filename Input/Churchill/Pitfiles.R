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
setwd("~/Documents/School/Projects/MEMLS/MEMLS3&a/Input/Churchill/Pitfiles")

matlabcli=paste("/Applications/MATLAB_R2014a.app/bin/matlab -nodisplay -nosplash -nodesktop -r \"cd('",getwd(),"');run('run_MEMLS_Active_v2(|pitfile|,|soilfile|,|outfile|,|freq|,|theta|,|m|,|q|,|database|)');exit;\"",sep="")

########################################
########################################
## Sourcefiles
########################################
########################################

source('~/Documents/School/Projects/Rfunctions/R functions env.R')
source('~/Documents/School/Projects/Paper-Peatmodel/MODELS/MironovModel/mironov2010_newCoeff2012.R')
source('~/Documents/School/Projects/Paper-Peatmodel/PAPER/AGU2012/mironovSims/JRT_functions.R')
#Have run the following:
source('~/Documents/School/Projects/Visual-Scatterometry/pitfiles.R')
#To produce:
load('~/Documents/School/Projects/Paper-PeatModel/PAPER/AGU2012/mironovSims/pitfiles/pitfiles_tundra.rdata')

#Run inputbuilder portion of source('~/Documents/School/Projects/Paper-PeatModel/AGU2012/mironovSims/AGU2012-MironovSims.R')
#This builds the following files:
#allsweobs=read.csv("~/Documents/School/Projects/Paper-PeatModel/PAPER/AGU2012/mironovSims/pitfiles/allsweobs.csv")

########################################
########################################
## Build MEMLS3a inputfiles
########################################
########################################

allsweobs=data.frame(date="s",Date=strptime("19990101-0000","%Y%m%d-%H%M"),sweA=1,sweB=1,sweE=1,sweT=1,tgnd=1,tair=1,tsnow=1,depthpit=1,depthmagna=1,pdsw=1,DHfrac=1,meanGS_dobs=1,meanGS_dobs_upper=1,meanGS_dobs_lower=1,sdGS_dobs_upper=1,sdGS_dobs_lower=1,soil_rough=1,file="")

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
		soil_rough=NA,
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
		soil_rough=NA,
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
		soil_rough=NA,
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
		soil_rough=NA,
		file="Pit_Tundra_DST_05711_B1_2011_02_26.xls"
		)
allsweobs$nor=normalize(allsweobs$sweT,datarange=c(0,50))
#allsweobs$Date=strptime(gsub("(.*?)(201[01]_\\([[:digit:]]{2})\\){}(.*?[.]xls)",allsweobs$date),"%d/%m/%Y")

allsweobs=allsweobs[order(allsweobs$Date),]

#write.csv(allsweobs,file="allsweobs.csv")

######################################
##	Adjust paramters for MEMLS
allsweobs$snowWet=0
allsweobs$salinity=0
allsweobs$corr_exp=0.5*allsweobs$meanGS_dobs*(1-(allsweobs$pdsw*1000/917.00)) #TODO: better transformation
allsweobs$soilWVF=0.3 #(m3/m3) TODO: real value
allsweobs$soil_rough=0.4*100 #(m), loadsoil.m divides this by 100 TODO: real value
allsweobs$soil_sdh=0.4 #(cm) standard deviation of soil height (wegmuller and matzler 1999 's' variable)


#subset to working site
(sweobs=allsweobs[rownames(allsweobs) %in% c("2010-308-S-TS","2010-319-S-TS", "2010-327-S-TS", "2010-340-S-TS", "2010-345-S-TS", "2010-352-S-TS", "2011-004-S-TS", "2011-016-S-TS", "2011-023-S-TS", "2011-033-S-TS", "2011-044-S-TS", "2011-047-S-TS", "2011-057-S-TS", "2011-056-A-TD", "2011-056-B-TD", "2011-062-A-TD", "2011-062-A-TD"),])


sweobs[sweobs$date=="04/11/2010","Tgnd_K5lh"]=-0.72
sweobs[sweobs$date=="15/11/2010","Tgnd_K5lh"]=0.033
sweobs[sweobs$date=="23/11/2010","Tgnd_K5lh"]=-5.987
sweobs[sweobs$date=="06/12/2010","Tgnd_K5lh"]=-5.409
sweobs[sweobs$date=="11/12/2010","Tgnd_K5lh"]=-8.66
sweobs[sweobs$date=="18/12/2010","Tgnd_K5lh"]=-2.841
sweobs[sweobs$date=="04/01/2011","Tgnd_K5lh"]=-5.787
sweobs[sweobs$date=="16/01/2011","Tgnd_K5lh"]=-7.16
sweobs[sweobs$date=="23/01/2011","Tgnd_K5lh"]=-11.37
sweobs[sweobs$date=="02/02/2011","Tgnd_K5lh"]=-9.76
sweobs[sweobs$date=="13/02/2011","Tgnd_K5lh"]=-12.34
sweobs[sweobs$date=="16/02/2011","Tgnd_K5lh"]=-11.09
sweobs[sweobs$date=="25/02/2011","Tgnd_K5lh"]=-11.01
sweobs[sweobs$date=="26/02/2011","Tgnd_K5lh"]=-11.29
sweobs[sweobs$date=="03/03/2011","Tgnd_K5lh"]=-12.41

sweobs$epsg=0+0i
sweobs$epsg= sapply(sweobs$Tgnd_K5lh,function(x)mironov(9.6,x,75.4/1000,0.51))

sweobs=subset(sweobs, depthpit>0)

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
	write.table(t(c(pit$Tgnd_K5lh+273.13,pit$soilWVF,pit$soil_rough)),file=gsub("\\.txt","_soil.txt",filename),sep=" ",row.names=F,col.names=F)
}



######################################
##	Run MEMLS

freq=c(10,8,18)  #GHz
theta=c(21,3,81) #degrees
MEMLS_m=0.075  #amemlsmain%     m:     mean slope of surface undulations (typical 0.05 to 0.1)
MEMLS_q=0.3    #amemlsmain%     q:     cross pol fraction (typical 0.1 to 0.3)
MEMLS_sccho=c(1,2,4,5,7,8,9,10,11,12) #amemlsmain%     sccho: type of scattering coefficient (11 recommended)
MEMLS_echo_to_prompt=0

#build matlab command chain.
RUN_ID=round(as.numeric(Sys.time()))
matlabcli="/Applications/MATLAB_R2014a.app/bin/matlab -nodisplay -nosplash -nodesktop -r \"cd('~/Documents/School/Projects/MEMLS/MEMLS3&a');\""
cli_input=""
for(s in MEMLS_sccho[10]){
	for(m in seq(0.05,1,length=3)){
		for(q in seq(0.1,0.3,length=3)){
			for(i in 1:nrow(sweobs)){
				pit=sweobs[i,]
				#print(pit)
				pitfilename=paste(c(row.names(pit),"_L1",".txt"),collapse="")
				soilfilename=gsub("\\.txt","_soil.txt",pitfilename)
				#Layer#, layerTemp (K), Vol LWC (0-1), PDSW (kg/m3), thickness (cm), salinity (ppthousand), exp corr len (mm) 
				#pit=read.table(pitfilename,sep=" ")
				#soil=read.table(soilfilename,sep=" ")
				
				outfilename=paste( RUN_ID,"__",gsub("\\.txt",paste("_MEMLS.txt",sep=""),pitfilename),sep="")#_m",m,"_q",q,"
				
				command="run_MEMLS_Active_v2_External('|pitfile|','|soilfile|','|outfile|',|freq_start|,|freq_step|,|freq_stop|,|theta_start|,|theta_step|,|theta_stop|,|m|,|q|,|sccho|,|echo_to_prompt|,|database|);"
				command=gsub("|pitfile|",pitfilename, command,fixed=T)
				command=gsub("|soilfile|",soilfilename,command,fixed=T)
				command=gsub("|outfile|",outfilename,command,fixed=T)
				command=gsub("|freq_start|",paste(freq[1],collapse=","),command,fixed=T)
				command=gsub("|freq_step|",paste(freq[2],collapse=","),command,fixed=T)
				command=gsub("|freq_stop|",paste(freq[3],collapse=","),command,fixed=T)
				command=gsub("|theta_start|",paste(theta[1],collapse=","),command,fixed=T)
				command=gsub("|theta_step|",paste(theta[2],collapse=","),command,fixed=T)
				command=gsub("|theta_stop|",paste(theta[3],collapse=","),command,fixed=T)
				command=gsub("|m|",m,command,fixed=T)
				command=gsub("|q|",q,command,fixed=T)
				command=gsub("|sccho|", s,command,fixed=T)
				command=gsub("|echo_to_prompt|", MEMLS_echo_to_prompt,command,fixed=T)
				command=gsub("|database|",0,command,fixed=T)
				
				cli_input =c(cli_input, command)
				#cat(command)
				
				}
			}
		}
	}
#Run matlab on list of commands: (Bulk running nessisary due to extremly SLOW matlab loading.)
cli_input =c(cli_input,"exit;")[-1]
cat("Running MATLAB on ", length(cli_input)," runs of MEMLS. Standby. (ps -A |grep MATLAB to find process if it needs to be terminated...)")
envlist=list(model="MEMLS",modelver=4.0,modeldate=as.numeric(difftime(Sys.time(),strptime(1970,"%Y"))),freq,theta,MEMLS_m,MEMLS_q,MEMLS_sccho,command,cli_input,RUN_ID,sweobs, matlabcli,curdir=getwd())
save(envlist,file=paste("../Modelruns/",RUN_ID,"__RUN_INFO.RData",sep=""));rm(envlist)

###########
## SLOW
system(matlabcli, input= cli_input)
###########


#Read results files
results=data.frame()
setwd("../Modelruns")
files=dir(pattern=".*_MEMLS.txt")

for(outfilename in files){
	MEMLSresults=read.csv(outfilename, header=F)
	names(MEMLSresults)=c("ModelVersion","ModelVersionDate","nlayer","freq","theta","soil.meanslope","Xpol.frac","sccho","snow.T1","snow.density","soil.T","soil.mv","soil.rough","s0h", "s0v", "ss0h", "ss0v", "rv","rh","rdv","rdh","rsv","rsh","rs0","sigma0vv","sigma0hh","sigma0hv","Tbv","Tbh")
	MEMLSresults$run_id=as.numeric(gsub("^([[:digit:]]+)_.*","\\1", outfilename))
	MEMLSresults$m=m
	MEMLSresults$q=q
	
	pit=sweobs[rownames(sweobs)==gsub("(.*?__)(.*?)(_L._MEMLS.txt)","\\2",outfilename),]
	MEMLSresults=data.frame(pit,MEMLSresults)
	results=rbind(results,MEMLSresults)
	}

getRunInfo=function(resultsFileName){
	run_id=as.numeric(gsub("^([[:digit:]]+)_.*","\\1", resultsFileName))
	#load this run's envlist object
	load(paste(run_id,"__RUN_INFO.RData",sep=""))
	return(envlist)
}

results$modelver="VerGIT"
results$model="MEMLS"
results$modelRunDate=results$run_id
class(results$modelRunDate)="POSIXct"
results$simvar1="m"
results$simvar2="q"
results$inc=results$incidence=results$theta
save(results,file=paste("Results_",format(Sys.time(), "%Y_%m_%d__%H_%M"),".Rdata",sep=""))



#fold results by ScatFrequency
results$key=factor(paste(results$file,results$m,results$q,results$soil.rough,results$soil.mv,results$sccho,results$simvar1,results$simvar2,results$theta,sep="@"))
newdataX=subset(results,freq==unique(results$freq)[1])
newdataK=subset(results,freq==unique(results$freq)[2])
names(newdataX)[grepl("sigma0..",names(newdataX))]=c("Xvv","Xhh","Xvh")
names(newdataK)[grepl("sigma0..",names(newdataK))]=c("Kvv","Khh","Kvh")
stopifnot(all(newdataX$key==newdataK$key))#Assure both are ordered correctly
newdata=cbind(newdataX[,!grepl("X[vh]{2}",names(newdataX))],newdataX[,grepl("X[vh]{2}",names(newdataX))],newdataK[grepl("K[vh]{2}",names(newdataK))])
newdata$key=NULL



plotSim.NRCS.vs.Date(newdata,plotOBS=F,plotOBSavg=F,simModel="MEMLS",simModelVer="VerGIT",varname="m")
plotSim.NRCS.vs.Inc.byDate(newdata,plotOBS=F,plotOBSavg=F,simModel="MEMLS",simModelVer="VerGIT")





plotNRCS_date=function(data,dates=NA){
	
}


plotNRCS=function(data,mode=c("date","inc")){
	#
	
	#Assumptions: 
	
	###############
	##	SETUP WINDOW
	if(mode=="date"){
		MIN_DATE=round(min(results$Date),units="days")
		MAX_DATE=round(max(results$Date),units="days")
		DATE_RANGE=seq(from= MIN_DATE,to= MAX_DATE,"days")
		
		NWINDOW= length(unique(data$Date)) #Number of panels in window. Track so we can add to the plot
		NWINDOW_C=NWINDOW_R= floor(sqrt(NWINDOW))+1
		par(mfrow=c(NWINDOW_R,NWINDOW_C))
		
		#Draw empty plots
		NWINDOW_LIST=""
		for(i in 1:NWINDOW){
			
		}
		
		plot(dates,rep(10,length(dates)),ylim=ylim,type="l" ,col="white",main=paste(freqpol[i,]$varname,moretitle,sep=""),xlab="Date",ylab=paste("Response (dB)"))#,panel.first=grid(nx=5), xaxt="n")

		abline(h=seq(-30,0,by=5), v=seq(from=as.Date("2010-11-01"), to=as.Date("2011-03-03"),"months"), lty="dotted", col="lightgrey",lwd=0.7)
	}
	if(mode=="inc"){
		
	}
	
	
	invisible(par)
}






plotSim.NRCS.vs.Date=function(bigtable, writePDF=F,varname="none",plotSWE=F,plotDHF=F,plotRHO=F,plotGZ=F,plotTGND=F,plotEPSG=F,plotTSX=T,plotOBS=T,plotSIM=T,plotOBSavg=T, plotSIMavg=F,plotSIMdeco=F,plotALLextra=F, simModel="sRT",simModelVer="VerSVN", plotLegend=c(F,F,F,T),INC=39,moretitle="",simLayers=1,ylim=c(-30,5)){#makeplot
	if(plotALLextra){
		plotSWE=plotDHF=plotRHO=plotGZ=plotTGND=plotTSX=plotOBS=plotSIMdeco=plotEPSG=T
	}
	
	#TODO: when this is loadedd into an env, grabbing the global sweobs doesnt always work. (dataM vs dataC)

	style=data.frame(name=c("obs flood","obs narrow","SRT","TSX","SWE","DHF","RHOS","GZ","TGND","SRTg","SRTas","SRTv","EPSGr","EPSGi"), disp=c("Obs.","Obs.","sRT","TSX","SWE","Depth Hoar frac.", "Snow density","Grain size","Temperature:GND","s0:ground","s0:air-snow","s0:snow vol","epsg.r","epsg.i"), note=c("","", "",""," (mm)"," (%)", " (g/cm3)"," (mm)", "(C)","","","","",""), lty=c("dashed","12","solid","solid","dotdash","dotdash","dotdash","dotted","dotted","solid","solid","solid","dotted","dotted"), pch=c(7,0,22,1,1,9,8,10,3,2,6,5,7,9), lwd=c(1,1,1,1,.5,.5,.5,.5,.5,.8,.8,.8,.5,.5), col=c("black","black","black","grey","red","green","cyan","orange","darkgreen","black","black","black","darkgreen","darkgreen"), cex=c(0.5,0.5,0.8,0.8,.8,.8,.8,.8,.8,.8,.8,.8,.5,.5),stringsAsFactors=F) #For consistent plotting styles within this func. #var="";s=style[which(style$name==var),]; plot...(s$pch, etc.)
	style=style[-1,]#disable obs flood.
	freqpol=data.frame(var=c("Xvv","Xvh","Kvv","Kvh"),varname=c("X-VV","X-VH","Ku-VV","Ku-VH"),col=c("black","black","black","black"),stringsAsFactors=F)
	#dates=strptime(as.character(apply(rbind(expand.grid(a=c(2010),b=c(11,12),c=seq(1,31)),expand.grid(a=c(2011),b=c(1,2,3),c=seq(1,31))),1,function(x)paste(x,collapse="-"))),"%Y-%m-%d")
	#dates=sort(dates[!is.na(dates)])[13:124]
	#dates=dates[!is.na(dates)]
	dates=seq(from=as.Date("2010-11-04"),to=as.Date("2011-03-03"),"days")
	
	##Fill note for legend name incidence angle
	style[grepl("obs",style$name),]$note=ifelse(plotOBSavg," over 30-45º", paste(" at ",INC,"º",sep=""))
	style[grepl("SRT",style$name),]$note=ifelse(plotSIMavg," over 30-45º", paste(" at ",INC,"º",sep=""))

	
	##Subset data
		y=subset(bigtable,modelver=="obs" & model=="obs_narrow" & incidence==INC[1])
	if(plotOBSavg & any(is.na(INC)))
		y=subset(bigtable,modelver=="obs" & model=="obs_narrow" & inc=="AVG")
	
	if(plotSIM){
	sims=subset(bigtable,modelver==simModelVer & model==simModel & incidence==INC[1] & simvar1==varname & (nlayer==simLayers | nlayer==0))

	if(plotSIMavg)
		sims=subset(avgInc(subset(bigtable,modelver==simModelVer & model==simModel & simvar1==varname)),inc=="AVG")
	if((plotSIMavg&plotSIM)+(plotOBSavg&plotOBS)==1)
		warning("Only one of SIMs and OBS presented as an average.")
		
	if(nrow(sims)==0)stop("No sims!")
	}else{sims=NULL}
	
	par(mfrow=c(2,2))
	for(i in 1:4){#for each frequency:pol...
		
		##Setup plot
		plot(dates,rep(10,length(dates)),ylim=ylim,type="l" ,col="white",main=paste(freqpol[i,]$varname,moretitle,sep=""),xlab="Date",ylab=paste("Response (dB)"))#,panel.first=grid(nx=5), xaxt="n")

		abline(h=seq(-30,0,by=5),v=seq(from=as.Date("2010-11-01"),to=as.Date("2011-03-03"),"months"),lty="dotted",col="lightgrey",lwd=0.7)

		##Plot SWE
		if(plotSWE){
			var="SWE";s=style[which(style$name==var),]
				points(sweobs$Date,normalize(sweobs$sweT,datarange=c(0,50),ylim=ylim),type="b",pch=s$pch,cex=s$cex,lwd=s$lwd,lty=s$lty,col=s$col)
			at=axTicks(4)
			at.l=as.character(paste(round(seq(0,50,length=length(at)),0),"")) 
			axis(4,at,labels=at.l)
			}
		#stop()
		##Plot DHF
		if(plotDHF){
			var="DHF";s=style[which(style$name==var),]
				points(as.Date(sweobs$Date),normalize(sweobs$DHfrac,datarange=c(0,1),ylim=ylim),pch=s$pch,cex=s$cex,type="b",lwd=s$lwd,lty=s$lty,col=s$col)
		}
		
		##Plot Density
		if(plotRHO){
			var="RHOS";s=style[which(style$name==var),]
				points(as.Date(sweobs$Date),normalize(sweobs$pdsw,datarange=c(0.1,0.3),ylim=ylim),pch=s$pch,cex=s$cex,type="b",lwd=s$lwd,lty=s$lty,col=s$col)
		}
		
		##Plot Grain size
		if(plotGZ){
			var="GZ";s=style[which(style$name==var),]
				points(as.Date(sweobs$Date),normalize(sweobs$meanGS_dobs,datarange=c(0.2,4),ylim=ylim),pch=s$pch,cex=s$cex,type="b",lwd=s$lwd,lty=s$lty,col=s$col)
		}
		
		##Plot Temperature:GND
		if(plotTGND){
			var="TGND";s=style[which(style$name==var),]
				points(as.Date(sweobs$Date),normalize(abs(sweobs$tgnd),ylim=ylim,datarange=c(0,30),invert=T),pch=s$pch,cex=s$cex,type="b",lwd=s$lwd,lty=s$lty,col=s$col)
		}
		if(plotEPSG){
			var="EPSGr";s=style[which(style$name==var),]
				points(as.Date(sweobs$Date),normalize(abs(Re(sweobs$epsg)),ylim=ylim,datarange=c(0,6),invert=T),pch=s$pch,cex=s$cex,type="b",lwd=s$lwd,lty=s$lty,col=s$col)
			var="EPSGi";s=style[which(style$name==var),]
				points(as.Date(sweobs$Date),normalize(abs(Im(sweobs$epsg)),ylim=ylim,datarange=c(0,2),invert=T),pch=s$pch,cex=s$cex,type="b",lwd=s$lwd,lty=s$lty,col=s$col)
		}
		
		##Plot TSX
		if(exists("tsxtundra")& grepl("[xX]",freqpol[i,"var"]) & plotTSX){
			var="TSX";s=style[which(style$name==var),]
				points(as.Date(tsxtundra$Date), tsxtundra[,freqpol[i,]$var],pch=s$pch,type="p",ylim=c(-30,5),col=s$col,cex=s$cex,lwd=s$lwd) #hollow circle
				points(as.Date(tsxtundra$Date), tsxtundra[,freqpol[i,]$var],pch=s$pch,type="l",ylim=c(-30,5),col=s$col,cex=s$cex,lwd=s$lwd) #hollow circle
			}else{plotTSX=F}#if no TSX, don't include in legend.
		
		##Plot OBS
		if(plotOBS){
			var="obs narrow";s=style[which(style$name==var),]
				for(iinc in INC){
					if(length(INC)>1){
						y=subset(bigtable,modelver=="obs" & model=="obs_narrow" & incidence==iinc)
						s$col=rainbow(n=81,start=0.1)[iinc]
							}
					points(as.Date(y[!is.na(y[,freqpol[i,]$var]),]$Date),y[!is.na(y[,freqpol[i,]$var]),freqpol[i,]$var],pch=s$pch,type="b",ylim=c(-30,5),bg=freqpol[i,"col"],col=s$col,lty=s$lty,lwd=s$lwd,cex=s$cex)#ifelse(length(INC)==1,s$col,pal[iinc])
					}
			}
		
		##Plot Sim
		if(plotSIM){
		var="SRT";s=style[which(style$name==var),]
			lines(as.Date(sims$Date), sims[,freqpol[i,]$var],ylim=c(-30,5),col=s$col,lwd=s$lwd)
			points(as.Date(sims$Date), sims[,freqpol[i,]$var],pch=s$pch,type="p",ylim=c(-30,5),col=s$col,bg=freqpol[i,"col"],cex=s$cex)
		
		##Plot Decomposed Sim values
		if(plotSIMdeco){
			f=ifelse(grepl("X",freqpol[i,]$var),"X","K")
			p=ifelse(grepl("vv",freqpol[i,]$var),"vv","vh")
						
			var="SRTg";s=style[which(style$name==var),]
			points(as.Date(sims$Date),sims[,paste(f,p,"g",sep="")],pch=s$pch,type="b",ylim=c(-30,5),bg=freqpol[i,"col"],col=s$col,lty=s$lty,lwd=s$lwd,cex=s$cex)
			
			var="SRTas";s=style[which(style$name==var),]
			points(as.Date(sims$Date),sims[,paste(f,p,"as",sep="")],pch=s$pch,type="b",ylim=c(-30,5),bg=freqpol[i,"col"],col=s$col,lty=s$lty,lwd=s$lwd,cex=s$cex)
			#style=style[style$name!="SRTas",]#clip out as we are not using it.
			
			var="SRTv";s=style[which(style$name==var),]
			points(as.Date(sims$Date),sims[,paste(f,p,"v",sep="")],pch=s$pch,type="b",ylim=c(-30,5),bg=freqpol[i,"col"],col=s$col,lty=s$lty,lwd=s$lwd,cex=s$cex)
	
				
			}
			}
				
		
		##Plot legend
		
		if(!plotSWE)style=style[style$name!="SWE",]
		if(!plotTSX)style=style[style$name!="TSX",]
		if(!plotDHF)style=style[style$name!="DHF",]
		if(!plotRHO)style=style[style$name!="RHOS",]
		if(!plotGZ)style=style[style$name!="GZ",]
		if(!plotOBS)style=style[style$name!="obs narrow",]
		if(!plotSIM)style=style[style$name!="SRT",]
		if(!plotEPSG)style=style[style$name!="EPSGr" & style$name!="EPSGi",]
		if(!plotTGND)style=style[style$name!="TGND",]
		if(!plotSIMdeco){style=style[style$name!="SRTg",];style=style[style$name!="SRTv",];style=style[style$name!="SRTas",]}
		if(!plotOBS)style=style[!grepl("obs",style$name),]
		legendloc=c("bottomright","topleft","bottomright","topleft")
		if(plotLegend[i])
			legend(x=legendloc[i],legend=paste(style$disp,style$note,sep=""),lty=style$lty,pch=style$pch,pt.lwd=style$lwd,bg="white",cex=1,pt.bg="black",col=style$col)#draw legend
	
		}#end for freq:pol
	
	if(writePDF)dev.copy2pdf(file=paste("Qpol_obs_vs_Simulated_inc39_tsx","_",moretitle,".pdf",sep=""))
	invisible(list(y=y,z=sims,s=style))
	}
##############################
#temp=plotSim.NRCS.vs.Date(rbind(obs,paramSweeps))
#data=dataC;temp=plotSim.NRCS.vs.Date(data,writePDF=F,plotOBSavg=T,plotSIMavg=F,plotALLextra=T)
#data=dataM;temp=plotSim.NRCS.vs.Date(data,writePDF=F,plotOBSavg=T,plotSIMavg=F,plotALLextra=T)















