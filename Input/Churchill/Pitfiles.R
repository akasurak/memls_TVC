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

matlabcli=paste("/Applications/MATLAB_R2014a.app/bin/matlab -nodisplay -nosplash -nodesktop -r \"cd('",getwd(),"');run('run_MEMLS_Active_v2(|pitfile|,|soilfile|,|outfile|,|freq|,|theta|,|m|,|q|,|database|)');exit;\"",sep="")

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

write.csv(allsweobs,file="allsweobs.csv")

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
m=0.075  #amemlsmain%     m:     mean slope of surface undulations (typical 0.05 to 0.1)
q=0.3    #amemlsmain%     q:     cross pol fraction (typical 0.1 to 0.3)

results=data.frame()
for(i in 1:nrow(sweobs)){
	pit=sweobs[i,]
	print(pit)
	pitfilename=paste(c(row.names(pit),"_L1",".txt"),collapse="")
	soilfilename=gsub("\\.txt","_soil.txt",pitfilename)
	#Layer#, layerTemp (K), Vol LWC (0-1), PDSW (kg/m3), thickness (cm), salinity (ppthousand), exp corr len (mm) 
	#pit=read.table(pitfilename,sep=" ")
	#soil=read.table(soilfilename,sep=" ")
	
	outfilename=gsub("\\.txt","_MEMLS3a_VV.txt",pitfilename)
	
	matlabcli=paste("/Applications/MATLAB_R2014a.app/bin/matlab -nodisplay -nosplash -nodesktop -r \"cd('~/Documents/School/Projects/MEMLS/MEMLS3&a');run_MEMLS_Active_v2_External('|pitfile|','|soilfile|','|outfile|',|freq|,|theta|,|m|,|q|,|database|);exit;\"",sep="")
	command=gsub("|pitfile|",pitfilename,matlabcli,fixed=T)
	command=gsub("|soilfile|",soilfilename,command,fixed=T)
	command=gsub("|outfile|",outfilename,command,fixed=T)
	command=gsub("|freq|",paste(freq,collapse=","),command,fixed=T)
	command=gsub("|theta|",paste(theta,collapse=","),command,fixed=T)
	command=gsub("|m|",m,command,fixed=T)
	command=gsub("|q|",q,command,fixed=T)
	command=gsub("|database|",0,command,fixed=T)
	
	cat(command)
	
	system(command)
	
	MEMLSresults=read.csv(outfilename, header=F)
	names(MEMLSresults)=c("freq","theta","s0h", "s0v", "ss0h", "ss0v", "rv","rh","rdv","rdh","rsv","rsh","rs0","sigma0vv","sigma0hh","sigma0hv","Tbv","Tbh")
	MEMLSresults=data.frame(pit,MEMLSresults)
	results=rbind(results,MEMLSresults)
}


results$modelver="VerGIT"
results$model="MEMLS3a"
results$simvar1=results$simvar2="none"
results$nlayer=1


names(simtable)[1:17]=c( "filename","datenum", "modelnum", "modelver", "incidence", "nlayer", "varnum1", "varnum2", "simvalue1", "simvalue2", "Xvv", "Xvh", "Kvv", "Kvh", "epsgRE", "epsgIM", "experiment")

results[results$freq==10,"Xvv"]=results[results$freq==10,"sigma0vv"]
results[results$freq==10,"Xvh"]=results[results$freq==10,"sigma0hv"]
results[results$freq==18,"Kvv"]=results[results$freq==18,"sigma0vv"]
results[results$freq==18,"Kvh"]=results[results$freq==18,"sigma0hv"]






plotSim.NRCS.vs.Inc.byDate=function(bigtable,dates=as.character(unique(bigtable$Date)), writePDF=F,plotSWE=F,plotDHF=F,plotRHO=F,plotGZ=F,plotTGND=F,plotEPSG=F,plotTSX=F,plotOBS=T,plotSIMdeco=T,plotALLextra=F, plotLegend=c(F,F,F,T),moretitle=""){ #makeplot2
	
	style_pol = data.frame(var = c("Xvv", "Xvh", "Kvv", "Kvh"), varname = c("X-VV", "X-VH", "Ku-VV", "Ku-VH"), col = c("black", "black", "black", "black"), stringsAsFactors = F)

style = data.frame(name = c("obs flood", "obs narrow", "SRT", "TSX", "SWE", "DHF", "RHOS", "GZ", "TGND", "SRTg", "SRTas", "SRTv", "EPSGr", "EPSGi"), 		disp = c("Obs.", "Obs.", "sRT", "TSX", "SWE", "Depth Hoar frac.", "Snow density", "Grain size", "Temperature:GND", "s0:ground", "s0:air-snow", "s0:snow vol", "epsg.r", "epsg.i"), note = c("", "", "", "", " (mm)", " (%)", " (g/cm3)", " (mm)", "(C)", "", "", "", "", ""), lty = c("dashed", "12", "solid", "solid", "dotdash", "dotdash", "dotdash", "dotted", "dotted", "solid", "solid", "solid", "dotted", "dotted"), pch = c(7, 			0, 22, 1, 1, 9, 8, 10, 3, 2, 6, 5, 7, 9), lwd = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.8, 0.8, 0.8, 0.5, 0.5), col = c("black", "black", 			"black", "grey", "red", "green", "cyan", "orange", "darkgreen", "black", "black", "black", "darkgreen", "darkgreen"), cex = c(0.5, 0.5, 			0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.5, 0.5), stringsAsFactors = F) #For consistent plotting styles within this func. #var=\\";s=style[which(style$name==var),]; plot...(s$pch, etc.)"
	style = style[-1, ] #disable obs flood.

	
	incs=seq(min(c(21,bigtable$theta),na.rm=T),81,by=3)
	ylims=c(min(c(-30,bigtable$Xvv,bigtable$Kvv),na.rm=T),max(c(5,unlist(bigtable[,c("Xvv","Xvh","Kvv","Kvh")])),na.rm=T))

	
		
	for(ii in dates){ #for each date supplied, or all dates if none specified...
		par(mfrow=c(2,2))
		for( i in 1:4){ #for each of Xvv, Xvh, Kvv, Kvh...

			##subset data
			obsF=subset(bigtable,modelver=="obs" & model=="obs_flood" &date==ii & theta!="AVG")
			obsN=subset(bigtable,modelver=="obs" & model=="obs_narrow" &date==ii & theta!="AVG")
	
			sims2=sims1=subset(bigtable,modelver=="VerGIT" & model=="MEMLS3a"  &(date==ii |is.na(date)) &simvar1=="none" & nlayer==1 & theta!="AVG")

			
			
			if(  (nrow(obsN)==0 &plotOBS)  | ((nrow(sims1)==0 & nrow(sims2)==0)  )   ){
				stop("No rows in obs or sims (check) when expected.")
			}
			


			##make plot
			plotSim.NRCS.vs.Inc.singleplot(obsF,obsN,sims1,sims2,style,style_pol[i,],incs,ylims,plotSWE,plotDHF,plotRHO,plotGZ,plotTGND,plotEPSG,plotTSX,plotOBS,plotSIMdeco,plotALLextra, plotLegend[i],title=paste(ii," ", style_pol[i,"varname"],moretitle,sep=""),legendloc=c("bottomright","topleft","bottomright","topleft")[i])
		}#end for freq:pol

		if(writePDF)dev.copy2pdf(file=paste("Qpol_vs_inc__",ii,"_",moretitle,".pdf",sep=""))
		}#end for freq:pol
	
	invisible(list(obsf=obsF,obsn=obsN,sims_n1=sims1,sims_n2=sims2))
	}#end for dates
##############################

plotSim.NRCS.vs.Inc.byDate(results,unique(results$Date))
