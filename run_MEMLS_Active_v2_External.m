function  MEMLS_result = run_MEMLS_Active_v2(pitfile,soilfile,outfile,freq_start,freq_step,freq_stop,theta_start,theta_step,theta_stop,m,q,database)
%     Environment Canada MEMLS_Active Framework
%     Date: 15/05/14
%     Authors: J. King, B. Montpetit (see sub-scripts for respective
%     authorship and acknolwedgement)

%     Input:
%     freq:     frequuency [in GHz]
%     theta:    incident angle [in Deg]
%     m:     mean slope of surface undulations (typical 0.05 to 0.1) 
%     q:     cross pol fraction (typical 0.1 to 0.3)
%     database: save records to mysql database (0 or 1)

%     Output:
%     result.reflec: layer reflectivities from bottom to top
%     result.sigSoilma0: sigSoilma naught in vv, hh, and hv polarizations
%     result
%
%     Uses: loadsnowpit, loadsoil


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
snowInput = loadsnowpit(fullfile(pwd,'Input','Churchill',pitfile));
soilInput = loadsoil(fullfile(pwd,'Input','Churchill',soilfile)); 
%soilInput = loadsoil(['Soil_' fName(a).name]);

freqList=freq_start:freq_step:freq_stop;
thetaList=theta_start:theta_step:theta_stop;

%delete(fullfile(pwd,'Input','Churchill',outfile));

   for freqIndex = 1:length(freqList)
       for thetaIndex = 1:length(thetaList)
       
		[s0h, s0v, ss0h, ss0v] = snowsoilreflectivity(freqList(freqIndex), thetaList(thetaIndex), snowInput.Ti(1), snowInput.roi(1), soilInput.t, soilInput.mv,  soilInput.sig);

				 %Note: input for soil reflectivity needs to be addressed. 
		MEMLS_result(fIndex,freqIndex,thetaIndex) = amemlsmain(freqList(freqIndex),thetaList(thetaIndex),s0h,s0v,s0h/1.1,s0v/1.1,fullfile(pwd,'Input','Churchill',pitfile),13,soilInput.t,12,m,q); %amemlsmain result.sigma0 = [sigma0vv,sigma0hh,sigma0hv];
		
		MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0=10*log10(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0) %Convert to dB
		
		
		disp(['amemlsmain results (dB):   sigma0vv= ' num2str(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0(1)) ' sigma0hh = ' num2str(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0(2)) ' sigma0hv = ' num2str(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0(3))]); 
		
	dlmwrite(fullfile(pwd,'Input','Churchill',outfile),[freqList(freqIndex), thetaList(thetaIndex),  s0h, s0v, ss0h, ss0v, MEMLS_result(fIndex,freqIndex,thetaIndex).reflec, MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0, MEMLS_result(fIndex,freqIndex,thetaIndex).Tb],'delimiter',',','-append');
	
		end
	end

%results=struct2table(MEMLS_result)
%csvwrite(fullfile(pwd,'Input','Churchill',outfile), results )

%M






%hist(10*log10(MEMLS_Sig0VV))
MEMLS_result
end