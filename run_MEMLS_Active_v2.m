function  result = run_MEMLS_Active_v2(freq,theta,m,q,database)
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

%Gather all filenames with a single prefix within the input folder
%TODO Allow prefix input from commandline
fName=dir(fullfile(pwd,'Input','TVC*.txt'))
fNum=size(fName);


for fIndex=1:fNum(1,1)
   
    snowInput = loadsnowpit(fullfile(pwd,'Input',fName(fIndex).name));
    soilInput = loadsoil(fullfile(pwd,'Input','Soil_TVC.txt')); 
    %soilInput = loadsoil(['Soil_' fName(a).name]);
  
   for freqIndex = 1:length(freq)
       for thetaIndex = 1:length(theta)
         [s0h, s0v, ss0h, ss0v] = snowsoilreflectivity(freq(freqIndex), theta(thetaIndex), snowInput.Ti(1), snowInput.roi(1), soilInput.t, soilInput.mv,  soilInput.sig);
    
           
         %Note: input for soil reflectivity needs to be addressed. 
           MEMLS_result(fIndex,freqIndex,thetaIndex) = amemlsmain(freq(freqIndex),theta(thetaIndex),s0h,s0v,s0h/1.1,s0v/1.1,fullfile(pwd,'Input',fName(fIndex).name),13,soilInput.t,11,m,q);
           MEMLS_Sig0VV(fIndex) = MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0(1);
           
           10*log10(MEMLS_result(fIndex,freqIndex,thetaIndex).sigma0)
       end
   end
end

hist(10*log10(MEMLS_Sig0VV))
