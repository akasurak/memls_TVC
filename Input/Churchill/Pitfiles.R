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
## Regenerate Observation files (TODO: fix inputbuilder sourcefiles)
########################################
########################################

obs=loadOBSTable("~/Documents/School/Projects/Paper-Peatmodel/DATA/SimOBS-ALL_FW__results_2012-12-02_2057=AGU2012_newpit_TGNDker_coreh2o.csv")#loadOBSTable("~/Documents/School/Projects/Paper-Peatmodel/PAPER/AGU2012/mironovSims/results_2012-12-02_2057=AGU2012_newpit_TGNDker_coreh2o/SimOBS-ALL_FW.csv")

obs=obs[!duplicated(obs),]#grab first instance of each



obs=avgInc(obs)#add a value for the average sigma0 over 30--45 degrees.

#For decompositions:
#temp=expand.grid(c("g","as","v","t"),c("vv","vh","hh"),c("X","K"))
#temp=paste(as.character(temp$Var3),as.character(temp$Var2),as.character(temp$Var1),sep="")
#for(i in temp)
#	obs[,i]=NA
	
allobs=obs

obsfeb25=obs[grepl("^(25022011-15[13]0)",obs$filename),]
obsfeb25SR=obs[grepl("^(25022011-1600)",obs$filename),]
obsSR=obs[grepl("^(03032011-1330)",obs$filename),]
obs=obs[!grepl("^(03032011-1330)|(25022011-)",obs$filename),]

OBS=obs

temp=plotSim.NRCS.vs.Date(obs,plotSIM=F)

###########################
###########################

tsx=read.csv("~/Documents/School/Projects/Paper-Peatmodel/DATA/TundraTSX_jul2012.csv")
names(tsx)=gsub("[.]","",tolower(names(tsx)))
tsx$Date=strptime(tsx$date,"%y-%m-%d")
TSXmean=data.frame(date=unique(tsx$date),Date=unique(tsx$Date),vv2x2=NA,vv3x3=NA,vh2x2=NA,vh3x3=NA)
for(i in 1:nrow(TSXmean)){
	for(j in names(TSXmean)[3:6])
	TSXmean[i,j]=mean(tsx[tsx$date==TSXmean[i,"date"],j],na.rm=T)
}


tsx3x3=function(){
	tsxtundra<<-TSXmean[,c("Date","date","vv3x3","vh3x3")]
	names(tsxtundra)<<-c("Date","date","Xvv","Xvh")
}

tsx2x2=function(){
	tsxtundra<<-TSXmean[,c("Date","date","vv2x2","vh2x2")]
	names(tsxtundra)<<-c("Date","date","Xvv","Xvh")
}
tsx2x2()
tsx3x3()


###########################
###########################


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





########################################
########################################
## Build MEMLS3a inputfiles
########################################
########################################



######################################
##OneLayer
for(i in 1:nrow(sweobs)){
	pit=sweobs[i,]
	print(pit)
	if(pit$depthpit<=0){
		print("!!!!!!!!!!!!!!!!!!	NO SNOW, NOT GENERATING PITFILE		!!!!!!!!!!!!!!!!!!")
		#this is because MEMLS currently requires snow (go fig.)
		next
	}
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
	for(m in seq(0.05,0.1,length=3)){
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
if(F)system(matlabcli, input= cli_input)
###########


#Read results files
results=data.frame()
setwd("../Modelruns")
files=dir(pattern=".*_MEMLS.txt")

for(outfilename in files){
	MEMLSresults=read.csv(outfilename, header=F)
	names(MEMLSresults)=c("ModelVersion","ModelVersionDate","nlayer","freq","theta","soil.meanslope","Xpol.frac","sccho","snow.T1","snow.density","soil.T","soil.mv","soil.rough","s0h", "s0v", "ss0h", "ss0v", "rv","rh","rdv","rdh","rsv","rsh","rs0","sigma0vv","sigma0hh","sigma0hv","Tbv","Tbh")
	MEMLSresults$run_id=as.numeric(gsub("^([[:digit:]]+)_.*","\\1", outfilename))
	
	pit=sweobs[rownames(sweobs)==gsub("(.*?__)(.*?)(_L._MEMLS.txt)","\\2",outfilename),]
	pit$name=row.names(pit)
	row.names(pit)=NULL
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
save(results,file=paste(results$run_id[1],"_Results_",format(Sys.time(), "%Y_%m_%d__%H_%M"),".Rdata",sep=""))



#fold results by ScatFrequency
results$key=factor(paste(results$file,results$m,results$q,results$soil.rough,results$soil.mv,results$sccho,results$simvar1,results$simvar2,results$theta,sep="@"))
newdataX=subset(results,freq==unique(results$freq)[1])
newdataK=subset(results,freq==unique(results$freq)[2])
names(newdataX)[grepl("sigma0..",names(newdataX))]=c("Xvv","Xhh","Xvh")
names(newdataK)[grepl("sigma0..",names(newdataK))]=c("Kvv","Khh","Kvh")
stopifnot(all(newdataX$key==newdataK$key))#Assure both are ordered correctly
newdata=cbind(newdataX[,!grepl("X[vh]{2}",names(newdataX))],newdataX[,grepl("X[vh]{2}",names(newdataX))],newdataK[grepl("K[vh]{2}",names(newdataK))])
newdata$key=NULL


#Add in Obs
obs$file=obs$filename
cat("Convert from factor:\n")
for(i in names(obs))if(any(class(obs[,i])=="factor")){print(i);obs[,i]=as.character(obs[,i])} 

for(i in names(sweobs)[!names(sweobs)%in%names(obs)]) obs[,i]=NA
for(i in 1:nrow(sweobs)){
	pit=sweobs[i,]
	obs[format(obs$Date,"%Y-%j")==substr(rownames(pit),1,8),names(sweobs)]=pit
}
obs$soil.T=obs$Tgnd_K5lh #or obs$tgnd
#obs$soil.mv =obs$soilWVF
#obs$epsg =complex(obs$)
obs$soil.rough =obs$soil_rough
obs$snow.T1 =obs$tsnow
obs$snow.density =obs$pdsw
obs$soil.meanslope =obs$soil_rough



for(i in names(newdata)[!names(newdata)%in%names(obs)]) obs[,i]=NA
obs=obs[,names(newdata)]
newdata=rbind(obs,newdata)

#Convert appropriate vars to factors:
cat("Convert to factor:\n")
for(i in c("file","inc")){print(i);newdata[,i]=factor(newdata[,i])}






















plotNRCS=function(data,X=c("date","inc"),Y=c("dB"),Z=c("FreqPol","date"),ylim=c(-30,5),plotGroup=c("OBS","ALLSIM"), plotLegend=T,legendOnOwnPanel=F, writePDF=F){
	def.par=par(no.readonly = TRUE) # save default, for resetting...
	#Assumptions: 
	stopifnot(Y=="dB")
	
	###############
	##	SETUP WINDOW
	par(mar=c(3,4,3,4)+0.1)
	par(mgp=c(2,1,0))
	if(any(grepl("FIELD",plotGroup)))
		par(par("mar")+c(0,0,0,1))
	
				par(mar=par("mar")+c(0,0,0,1))
	
	if(Z=="FreqPol"){
		FREQPOL=grep("[XK][vh]{2}",names(data),value=T)
		NWINDOW= length(FREQPOL)+ifelse(legendOnOwnPanel==2,1,0) #Number of panels in window. Track so we can add to the plot
		NWINDOW_R= floor(sqrt(NWINDOW))
		NWINDOW_C=ceiling(NWINDOW/NWINDOW_R) #prefer extra to be made up in columns.
		
		par(mfrow=c(NWINDOW_R,NWINDOW_C)) #draw by rows, nr * nc
		
	}
	if(Z=="date"){
		stop()
		title=paste("NRCS ", format(DATES[i],"%Y-%m-%d"),sep="")
	}
	
	plotlist=list(X=X,Y=Y,Z=Z, NWINDOW= NWINDOW, NWINDOW_R= NWINDOW_R, NWINDOW_C= NWINDOW_C)
	
	NWINDOW_LIST=data.frame(title="",r=NA,c=NA)
	
	for(i in 1:(NWINDOW-ifelse(legendOnOwnPanel==2,1,0))){
			title=paste("NRCS (",FREQPOL[i],")",sep="")
			
			
			#Setup plot frame
			plot.new()
			
			

		if(X=="date"){
			data$date=as.Date(round(data$Date,units="days"))
			plotlist$MIN_X=MIN_DATE=round(min(data $Date),units="days")
			plotlist$MAX_X=MAX_DATE=round(max(data $Date),units="days")
			plotlist$X_RANGE=DATE_RANGE=as.Date(seq(from= MIN_DATE,to= MAX_DATE,"days"))
			plotlist$X_VALS=DATES=as.Date(sort(unique(round(data$Date,units="days"))))
			plotlist$xlim=xlim=as.Date(c(MIN_DATE,MAX_DATE))
			
			plotlist$LEGEND_LOC= LEGEND_LOC="bottomright"
					
			
			plot.window(ylim=ylim,xlim= xlim)
			axis.Date(1, at = as.Date(DATE_RANGE[which(unclass(as.POSIXlt(DATE_RANGE))$mday==15)]), format = "%b-'%y",line=1,lwd=0)#bottom axis
			axis.Date(1, at =as.Date(DATE_RANGE[which(unclass(as.POSIXlt(DATE_RANGE))$wday==6)]), format = "%V")#bottom axis
			mtext("Date",1,cex=par("cex.axis")*par("cex")*1.1,line=3)
			
			#add extra axis for FIELD: measurements
			if(any(grepl("FIELD",plotGroup))){
				axis(4,at=seq(ylim[2],ylim[1],length=11),labels=round(seq(1,0,length=11),1),line=0.3,lty="dotted")
				mtext("Scaled Field Measurements",4,line=2.2,cex=par("cex.axis")*par("cex"))
				}
			
			
			#add grid lines
			abline(v=as.Date(DATE_RANGE[which(unclass(as.POSIXlt(DATE_RANGE))$mday==1)]), lty="dotted", col="lightgrey",lwd=0.7)
			}
	
		if(X=="inc"){
			data=data[data$incidence<=90,]
			plotlist$MIN_X=MIN_INC=min(data$incidence)
			plotlist$MAX_X=MAX_INC=max(data$incidence)
			plotlist$X_RANGE=INC_RANGE=seq(from= MIN_INC,to= MAX_INC)
			plotlist$X_VALS=INCS=sort(unique(data$incidence))
			plotlist$xlim=xlim=c(MIN_INC,MAX_INC)
			
			plotlist$LEGEND_LOC= LEGEND_LOC="topright"
			
			plot.window(ylim=ylim,xlim= xlim)
			axis(1,at=INCS)
			mtext("Incidence Angle (deg)",1,cex=par("cex.axis")*par("cex")*1.1,line=2)
			
			abline(v=seq(from=20,to=80,by=10), lty="dotted", col="lightgrey",lwd=0.7)
			
			if(length(unique(data$date))==1)title=paste(title,unique(data$date))
			}
		
		title(main=title,line=0.5)
		box()
		NWINDOW_LIST[i,]=data.frame(title,par("mfg")[1],par("mfg")[2])#must be after plot so mfg is current.
		
		axis(2,seq(ylim[2],ylim[1],by=-5)) #left axis
		mtext(paste("Response (dB)"),side=2,cex=par("cex.axis")*par("cex"),line=2)

		
	}
	
	
	
	
	plotlist$NWINDOW_LIST= NWINDOW_LIST
	
	if(Y=="dB"){
		for(xi in 1:NWINDOW_C)
			for(yi in 1:NWINDOW_R){
				if(xi* NWINDOW_R+yi > NWINDOW)next #do not draw grid on unused panels
				par(mfg=c(yi,xi))
				#add grid lines
				abline(h=seq(min(ylim),max(ylim),by=5), lty="dotted", col="lightgrey",lwd=0.7)
				}	
	MIN_Y =plotlist$MIN_Y=min(ylim)
	MAX_Y =plotlist$MAX_Y=max(ylim)
	Y_VALS =plotlist$Y_VALS=Y_RANGE =plotlist$Y_RANGE=seq(MIN_Y,MAX_Y,by=1)
	
	plotlist$ylim=ylim
	}
	
	
	
	
	#########################
	## Plot OBS/SIM
	if(any(plotGroup=="ALLSIM")){
		plotGroup=unique(c(plotGroup,unique(data$model)[unique(data$model)!="OBS"]))
		plotGroup=sort(plotGroup[plotGroup!="ALLSIM"])
		}
	
	if(any(grepl("OBSALL:",plotGroup))){
				plotGroup=gsub("OBSALL","obs_narrow",plotGroup)
				plotGroup=c(plotGroup[grepl("obs_narrow",plotGroup)],gsub("narrow","flood", plotGroup[grepl("obs_narrow",plotGroup)]),plotGroup[!grepl("obs",plotGroup)])
				}
	
	legend=NULL
	for(datagroup in plotGroup){
		datagroupRaw=datagroup
		plotdata=NA
		
		### Attend to slicing plotGroup's
		if(grepl("DB",datagroup)){
			if(X!="inc")
				stop("not coded")
			
			scanID=gsub("(DB:)([[:digit:]]+)","\\2")
			stop("not coded")
		}
		if(grepl("OBS",datagroup)){
			if(grepl("OBS:N",datagroup)){
				datagroup=gsub("OBS","obs_narrow",datagroup)
				}else if(grepl("OBS:F",datagroup)){
				datagroup=gsub("OBS","obs_flood",datagroup)
				}else{
				datagroup=gsub("OBS","obs_narrow",datagroup)
			}
		}
		if(grepl(".*VAR:.*", datagroup)){ #eg. "MEMLS:VAR:incidence@39"
			#Has @:
			if(grepl("@", datagroup)){
				#Draw the trace at this value of VAR
				
				#Has @AVG<#-#>
				if(grepl("VAR:.*?@AVG", datagroup)){
					
					#draw the trace at this average value
					#eg datagroup="MEMLS:VAR:noscale:incidence:AVG<33-45>"
					avgstart=as.numeric(gsub(".*?@AVG<([[:digit:]]+)-([[:digit:]]+)>","\\1", datagroup))
					avgstop=as.numeric(gsub(".*?@AVG<([[:digit:]]+)-([[:digit:]]+)>","\\2", datagroup))
					
					modelgroup=gsub("(.*?)(:VAR:)(.+)","\\1", datagroup)
					datagroup=gsub("(.*?)(:.*?)@AVG.*","\\1_AVG\\2",datagroup)
					trace=gsub(".*?VAR:(.*?scale:)?(.+?)(:.*)?","\\2", datagroup)
					
	
					plotdata =data[data[, trace]>=avgstart & data[, trace]<=avgstop & data$model== modelgroup,]
					if(plotlist$X=="date"){
						temp =nlme::gsummary(plotdata,form=~model/as.character(Date),FUN=list(character=function(x)names(summary(as.factor(x)))[1]))
						temp2=nlme::gsummary(plotdata,form=~model/as.character(Date), FUN=list(numeric=function(x)dbfun(x,mean, na.rm=T),character=function(x)names(summary(as.factor(x)))[1]))	
						}else{
							stop()
							temp =nlme::gsummary(plotdata,form=~model/as.character(Date),FUN=list(character=function(x)names(summary(as.factor(x)))[1]))
							temp2=nlme::gsummary(plotdata,form=~model/as.character(Date), FUN=list(numeric=function(x)dbfun(x,mean, na.rm=T),character=function(x)names(summary(as.factor(x)))[1]))
							}
					
					
					plotdata=data.frame(temp[,!grepl("[XK][vh]{2}",names(temp))],temp2[,grepl("[XK][vh]{2}",names(temp))])
					plotdata[,trace]=999
					plotdata$model=paste(modelgroup,"_AVG",sep="")
					trace=999
					data=rbind(data,plotdata)
					
				}else{
					trace=gsub("^(.*?@)(.+)$","\\2", datagroup)
					datagroup =gsub("(.*?)(@.+)$","\\1", datagroup)
					}
			}#end if("@")
			else{
				trace=NA
			}
			
			#Has :*scale:
			if(grepl(".*?VAR:[[:alpha:]]+scale+:.*", datagroup)){
				varscale=gsub(".*?VAR:(.+?):.+","\\1", datagroup)
			}else{
				varscale=F
				datagroup=gsub("VAR:","VAR:noscale:",datagroup)#so var=.. can be consistant
			}
			
				
			
			var=gsub("(.*?VAR:)(.+?):(.*)","\\3", datagroup)
			
			datagroup=gsub("(.*?)(:VAR:)(.+)","\\1", datagroup)
			plotNRCS.addGroup(subset(data,model==datagroup),plotlist,drawRange=T) #plot ranges, no line, no legend
			if(!is.na(trace))plotdata=subset(data,as.character(data[,var])==trace & model==datagroup)
			}#end if("VAR")
		
		
		if(grepl("FIELD:.+",datagroup)){
			if(X!="date"){warning("Asked to plot ",datagroup," but X axis is not date. Skipping.");next}
			var=gsub("(FIELD:)(.+)(:.*)","\\2", datagroup)
			

			plotdata=subset(data,model=="obs_narrow")
			plotdata$model=var
		}
		
		
		if(is.null(nrow(plotdata)))plotdata=subset(data,model==datagroup)
		
		
		if(any(duplicated(plotdata[,X])) & !grepl("FIELD:.+",datagroupRaw)){
			warning("Multiple Yvals for each Xval, while plotting as/with a line. Plotting once as line at ???SEE FUNC??? value, then as range.")
			if(grepl("VAR:.*?incidence",datagroupRaw)|!grepl("VAR",datagroupRaw)  ){
				if(X=="date"){dedupe="inc"} else{ dedupe="date"}
			}else{dedupe=gsub(".*?VAR:(.*)","\\1", datagroupRaw) }#which is the plotvar
				
			
			plotNRCS.addGroup(plotdata,plotlist,drawRange=T)
			
			if(dedupe=="GROUPED:SOMEFUNC"){
				#yvals=tapply(yvals,INDEX=factor(xvals), FUN=function(x)max(x))
				plotdata=nlme::gsummary(plotdata,form=formula(paste("~",X)),FUN=list(numeric=function(x)dbfun(x,mean)))
				plotdata=plotdata[order(plotdata[,X]),]
				}
			else if(dedupe=="FIRST"){
				plotdata=plotdata[match(as.Date(plotlist$X_VALS), plotdata[,X]),]
				}
			
			else {#if(dedupe=="inc"|dedupe=="date"){ #recall dedupe is opposite of X
				vals=unique(plotdata[,dedupe])
				for(i in vals[-1]){
					temp=plotdata[plotdata[,dedupe]==i,]
					plotNRCS.addGroup(temp,plotlist)
					
					}
				plotdata=plotdata[plotdata[,dedupe]==vals[1],]
				}
			}
		
		## Plot this group
		sty=plotNRCS.addGroup(plotdata,plotlist)
		legend =rbind(legend,sty)
	}

	
	#########################
	##	LEGEND
	
	if(plotLegend){
		Horiz=ifelse(nrow(legend)<=2,T,F)#if the legend is small, go horizontal. Where 3 is arbitrary and depends on screen size, etc...
		if(legendOnOwnPanel){
			Horiz=F
			plot.new()
			par(mfg=c(NWINDOW_R,NWINDOW_C))
			
		}
		legend(x=LEGEND_LOC,legend=paste(legend$disp, legend $note,sep=" "),lty= legend $lty,pch= legend $pch,pt.lwd= legend $lwd,bg="white",cex=1,pt.bg="black",col= legend $col, horiz= Horiz)#draw legend
	}
	
	#########################
	##	Cleanup
	
	if(writePDF)dev.copy2pdf(file=paste("NRCS_by_",X,"_vs_",paste(plotGroup,collapse="-"),".pdf",sep=""))
	
	par(def.par)#reset graphics window parameters..
	
	#########################	
	##	Return plotting structure so other plots can be overlaid.
	invisible(plotlist)
}

 #temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05  |grepl("obs",model))), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS","MEMLS:VAR:incidence@AVG<33-45>"))




temp=plotNRCS(subset(newdata, incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("MEMLS:VAR:Xpol.frac","OBS","MEMLS:VAR: soil.meanslope" ))





















plotNRCS.addGroup=function(data,plotlist, drawRange=list(F,T,"logscale","linscale" )  ){
	stopifnot(plotlist$Y=="dB")
	
	rangeFactor=function(x)1 #default value, scale CEX to be unscaled.
	
	switch(class(drawRange), list={drawRange=F}, logical=drawRange, character={rangeFactor=switch(drawRange, logscale=function(x){normalize(data=log(1+x/min(x+1))/(2*max((log(1+x/min(x+1))))),ylim=c(0.1,1.3),datarange=c(0.1,1.3))}, linscale=function(x){normalize(data=x/max(2*x,na.rm=T),ylim=c(0.1,1.3),datarange=c(0.1,1.3))}  );drawRange=T}, `function`={rangeFactor=drawRange;		drawRange=T},stop())

	print(drawRange)
	
	

	#For consistent plotting styles within this func. #var=""; sty=style[which(style$name==var),]; plot...(sty$pch, etc.) ;; style=edit(style);dput(style)
	style= structure(list(name = c("obs_flood", "obs_narrow", "SRT", "TSX", 
"SWE", "DHF", "RHOS", "GZ", "TGND", "SRTg", "SRTas", "SRTv", 
"EPSG.Re", "EPSG.Im", "MEMLS", "TSANG", "DMRTML", "TSNOW"), disp = c("Obs. Flood", 
"Obs. Narrow", "sRT", "TSX", "SWE", "Depth Hoar frac.", "Snow density", 
"Grain size", "Temperature:GND", "s0:ground", "s0:air-snow", 
"s0:snow vol", "Soil.permittivity (Real)", "Soil.permittivity (Imag.)", 
"MEMLS", "DMRT-QMS-ML (LACEO)", "DMRT-ML (LGGE)", "Temperature:Snow"
), note = c("(dB)", "(dB)", "(dB)", "(dB)", " (mm)", " (%)", 
" (g/cm3)", " (mm)", "(C)", "(dB)", "(dB)", "(dB)", "", "", "(dB)", 
"(dB)", "(dB)", "(C)"), lty = c("dotted", "dotted", "solid", "solid", 
"dotdash", "dotdash", "dotdash", "dotted", "dotted", "solid", 
"solid", "solid", "dotted", "dotted", "solid", "solid", "solid", 
"dotted"), pch = c(0, 7, 22, 1, 1, 9, 8, 10, 3, 2, 6, 5, 10, 
9, 22, 22, 22, 3), lwd = c(2, 2, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 
0.8, 0.8, 0.8, 0.5, 0.5, 1, 1, 1, 0.5), col = c("grey", "grey", 
"black", "grey", "red", "green", "cyan", "orange", "darkgreen", 
"black", "black", "black", "darkgreen", "darkgreen", "black", 
"black", "black", "lightgreen"), cex = c(0.7, 0.7, 0.8, 0.8, 
0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 
0.8), type = c("o", "o", "o", "o", "o", "o", "o", "o", "o", "o", 
"o", "o", "o", "o", "o", "o", "o", "o"), norm_min = c(NA, NA, 
NA, NA, 0, 0, 0.1, 0.2, 0, NA, NA, NA, 0, 0, NA, NA, NA, 0), 
    norm_max = c(NA, NA, NA, NA, 50, 1, 0.3, 4, 50, NA, NA, NA, 
    6, 3, NA, NA, NA, 50), bg = c("", "", "green", "", "", "", 
    "", "", "", "", "", "", "", "", "red", "cyan", "yellow", 
    ""), var = c("", "", "", "", "sweT", "DHfrac", "pdsw", "meanGS_dobs", 
    "Tgnd_K5lh", "", "", "", "epsg", "epsg", "", "", "", "tsnow"
    )), .Names = c("name", "disp", "note", "lty", "pch", "lwd", 
"col", "cex", "type", "norm_min", "norm_max", "bg", "var"), row.names = c(NA, 
18L), class = "data.frame")
style[style$bg=="","bg"]=NA#style[style$bg=="","col"]

	if(grepl("_AVG",unique(data$model))){
		warning("add style for average")
		data$model=gsub("_AVG","",unique(data$model))
		style$note=paste(" Avg.",style$note,sep="")
	}

	style<<-style;warning("styledump",call.=F)
	
	sty=style[style$name==unique(data$model),]
	
	if(drawRange & sty$type!="p"){
		sty$type="p"
	}else if (drawRange==F & sty$type=="p"){
		print("force style:type:b")
		sty$type="b"
	}


	for(zi in 1:nrow(plotlist$NWINDOW_LIST)){
		title=plotlist$NWINDOW_LIST[zi,'title']
		par(mfg=c(plotlist$NWINDOW_LIST[zi,'r'],plotlist$NWINDOW_LIST[zi,'c'])) #Set plot sub-window

		if(plotlist$Y=="dB"){
			if(grepl("dB",sty$note)){#plotting obs/sim
				yvar=gsub("(NRCS .)([XK][vh]{2})(.+)","\\2",title)
				}else{#plotting field measurements
					data=data[!duplicated(data[,plotlist$X,]),]#grab one copy of data[] for each Xval.
					yvar=sty$var
					}
			}
		
		
		if(plotlist$X=="date"){
			xvar="Date"
			data$Date=as.Date(data$Date)
			}else if (plotlist$X=="inc"){
				xvar="incidence"
				}
		xvals=data[,xvar]
		
		if(!drawRange & !is.na(sty$norm_min)){
			yvals=normalize(data[,yvar],datarange=c(sty$norm_min,sty$norm_max),ylim=plotlist$ylim)
			}
		else{
			yvals=data[,yvar]
			}
		
		if(any(is.complex(yvals))){
			if(grepl("Im",yvar)){
				yvals=Im(yvals)
			}else{
				yvals=Re(yvals)
			}
			
		}
		if(all(is.na(yvals))){warning("No Yvals:",title);next;}
		
		

		points(yvals~xvals,lty=sty$lty,lwd=sty$lwd,col=sty$col,cex=ifelse(rep(drawRange,length(yvals)),rangeFactor(yvals),1)*sty$cex,pch=sty$pch,type=sty$type,bg=sty$bg)
		#title(sub=paste("subplot: ",paste(plotlist$NWINDOW_LIST[zi,],collapse=" ")))
				
			#if(title=="NRCS (Kvv)")stop()
			
			}
sty<<-sty;	invisible(sty) #for legend plotting.

}

#temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05  |grepl("obs",model))& date=="16/02/2011" ), X="inc", Y="dB", Z="FreqPol", plotGroup=c("MEMLS","OBS"))






















#testing:
if(F){#By DATE
	####all these should plot successfully...
	
	#Test plotting Sims: MEMLS
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 & incidence ==39), X="date", Y="dB", Z="FreqPol", plotGroup="MEMLS")
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 ), X="date", Y="dB", Z="FreqPol", plotGroup="MEMLS")
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 ), X="date", Y="dB", Z="FreqPol", plotGroup="MEMLS:VAR:incidence@39")
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 ), X="date", Y="dB", Z="FreqPol", plotGroup="MEMLS:VAR:logscale:incidence@39")
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 ), X="date", Y="dB", Z="FreqPol", plotGroup="MEMLS:VAR:linscale:incidence@39")
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05  |grepl("obs",model))), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS","MEMLS:VAR:incidence@AVG<33-45>"))
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05  |grepl("obs",model))), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS:VAR:incidence@AVG<33-45>","MEMLS:VAR:incidence@AVG<33-45>"))

	
	#test plotting obs
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup="OBS")
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model))), X="date", Y="dB", Z="FreqPol", plotGroup="OBS")
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("MEMLS","OBS"))
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model))), X="date", Y="dB", Z="FreqPol", plotGroup=c("MEMLS","OBS"))
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05| grepl("obs",model)), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS:VAR:incidence@39","MEMLS:VAR:incidence@39"))
	
	#test plotting overlays
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model) ), X="date", Y="dB", Z="FreqPol", plotGroup="FIELD:SWE:")
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS","MEMLS","FIELD:SWE:","FIELD:DHF:", "FIELD:RHOS:", "FIELD:GZ:", "FIELD:TGND:", "FIELD:EPSG.Re:", "FIELD:EPSG.Im:"), legendOnOwnPanel=F)
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS","MEMLS","FIELD:SWE:","FIELD:DHF:", "FIELD:RHOS:", "FIELD:GZ:", "FIELD:TGND:", "FIELD:EPSG.Re:", "FIELD:EPSG.Im:"), legendOnOwnPanel=T)
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS","MEMLS","FIELD:SWE:","FIELD:DHF:", "FIELD:RHOS:", "FIELD:GZ:", "FIELD:TGND:", "FIELD:EPSG.Re:", "FIELD:EPSG.Im:"), legendOnOwnPanel=2)
	
	#test plotting by sweepvar
	temp=plotNRCS(subset(newdata, (soil.meanslope==0.05 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("MEMLS","OBS"))#TODO: fix dedupe selection here.
	temp=plotNRCS(subset(newdata, (soil.meanslope==0.05 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("MEMLS:VAR:Xpol.frac","OBS" ))
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 | grepl("obs",model))&incidence==39), X="date", Y="dB", Z="FreqPol", plotGroup=c("OBS","MEMLS:VAR:soil.meanslope" ))
	}



if(F){ #By INC
	temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 &date=="16/02/2011"), X="inc", Y="dB", Z="FreqPol", plotGroup="MEMLS")
	temp=plotNRCS(subset(newdata, (Xpol.frac==0.2 & soil.meanslope==0.05 |grepl("obs",model))&date=="16/02/2011"), X="inc", Y="dB", Z="FreqPol", plotGroup=c("OBS","MEMLS"))

temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 ), X="inc", Y="dB", Z="FreqPol", plotGroup="MEMLS")
temp=plotNRCS(subset(newdata, Xpol.frac==0.2 & soil.meanslope==0.05 |grepl("obs",model)), X="inc", Y="dB", Z="FreqPol", plotGroup=c("MEMLS","OBS"))
	
	}





