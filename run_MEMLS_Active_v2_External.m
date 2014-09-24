function  MEMLS_result = run_MEMLS_Active_v2_External(pitfile,soilfile,outfile,freq_start,freq_step,freq_stop,theta_start,theta_step,theta_stop,m,q,sccho,echo_to_prompt,database)
%     Andrew's MEMLS_Active Framework
%     Date: 2014/09/24
%     Authors: A. Kasurak, modifying run_MEMLS_Active_v2 by J. King, B. Montpetit (see sub-scripts for respective authorship and acknolwedgement)

%     Input:
%	  pitfile: input snow layering description filename.
%	  soilfile: input soil description filename.
%	  outfile: filename to APPEND results to, .CSV format.
%     freq:     frequency [in GHz]
%     theta:    incident angle [in Deg]
	% ^^_(start|step|stop):		To sweep over the above in one run, give MATLAB style start:step:stop values.
%     m:     mean slope of surface undulations (typical 0.05 to 0.1) 
%     q:     cross pol fraction (typical 0.1 to 0.3)
%	  sccho: scattering coefficient selector for sccoeff.m
%     database: save records to mysql database (0 or 1)
%     echo_to_prompt (1): 1:Echo output to interface, 0: Print nothing

%     Output:
%     result.reflec: layer reflectivities from bottom to top
%     result.sigSoilma0: sigSoilma naught in vv, hh, and hv polarizations
%	  CSV file in 'outfile' containing all produced and consumed values.
%     result

%     Uses: loadsnowpit, loadsoil

if nargin<14
	echo_to_prompt=1 %write output to prompt
	database=0 %disable
	end

MEMLSversion=4.0; %as amemlsmain.m
[temp,MEMLSdate]=system('R --slave -e "cat(difftime(Sys.time(),strptime(\"1970\",\"%Y\")))"');%floor(now-datenum([1970 0 0 0 0 0])); %number of days since 1970/01/01, as an indication of which sub-version of MEMLS we are using (crossref with git). Done via R because matlab does dates funky, and we want to use this to match with R.

MEMLSdate=str2num(MEMLSdate);

% BEGIN Database connection (see tech note to replicate expected schema)
% TODO: Finish database output stucture, non-functional as is.
if database == 1
    dbname = '';
    username = '';
    password = '';
    driver = '';
    dburl = ['' dbname];
    javaclasspath('');
    conn = database(dbname, username, password, driver, dburl);
end 

fIndex=1;
snowInput = loadsnowpit(fullfile(pwd,'Input','Churchill', 'Pitfiles',pitfile));
soilInput = loadsoil(fullfile(pwd,'Input','Churchill', 'Pitfiles',soilfile)); 
%soilInput = loadsoil(['Soil_' fName(a).name]);

nLayers=size([snowInput]);
nLayers=nLayers(1);

freqList=freq_start:freq_step:freq_stop;
thetaList=theta_start:theta_step:theta_stop;



   for freqIndex = 1:length(freqList)
       for thetaIndex = 1:length(thetaList)
       
		[s0h, s0v, ss0h, ss0v] = snowsoilreflectivity(freqList(freqIndex), thetaList(thetaIndex), snowInput.Ti(1), snowInput.roi(1), soilInput.t, soilInput.mv,  soilInput.sig);

				 %Note: input for soil reflectivity needs to be addressed. 
		MEMLS_result(fIndex,freqIndex,thetaIndex) = amemlsmain(freqList(freqIndex),thetaList(thetaIndex),s0h,s0v,s0h/1.1,s0v/1.1,fullfile(pwd,'Input','Churchill', 'Pitfiles',pitfile),13,soilInput.t,sccho,m,q); %amemlsmain result.sigma0 = [sigma0vv,sigma0hh,sigma0hv];
		
		MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0=10*log10(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0); %Convert to dB
		
		
		if echo_to_prompt == 1
		disp(['amemlsmain results (dB):   sigma0vv= ' num2str(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0(1)) ' sigma0hh = ' num2str(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0(2)) ' sigma0hv = ' num2str(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0(3))]); 
		end
		
	%[commandStatus,unixtime]=system('date +"%s"');
		
		
	dlmwrite(fullfile(pwd,'Input','Churchill','Modelruns',outfile),[MEMLSversion, MEMLSdate, nLayers, freqList(freqIndex), thetaList(thetaIndex), m, q, sccho, snowInput.Ti(1), snowInput.roi(1), soilInput.t, soilInput.mv, soilInput.sig, s0h, s0v, ss0h, ss0v, MEMLS_result(fIndex,freqIndex,thetaIndex).reflec, MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0, MEMLS_result(fIndex,freqIndex,thetaIndex).Tb], 'delimiter',',','-append');
	
		end
	end
	
%MEMLS_result(1,1,1).Snow=[snowInput];
%MEMLS_result(1,1,1).Soil=[soilInput];
%MEMLS_result(1,1,1).m=m;
%MEMLS_result(1,1,1).q=q;
%MEMLS_result(1,1,1).sccho=sccho;

end