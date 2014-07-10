function [r_h_mod, r_v_mod , fresnel_h, fresnel_v] = snowsoilreflectivity(freq, theta, tSnowK, rhoSnow, tSoilK, mvSoil, sigSoil)
%   Reflectivity estimate for rough bare soil covered with snow
%   Date: 13/05/14
%   Author: J. King
%   Based on code by B. Montpetit (ss_ref_grass.m, ruffsoil.m) and
%   reflectivity approach detailed in Wagmuller & Matzler 2009
%   Uses: ruffsoil.m, huteps.m, tdgrmdm.m. epss.m, epsr.m, gammah.m
%   Date: 13/05/14
%
%   Input
%   freq: Centre frequency [in GHz]
%   theta: Incident angle at snow surface [in Deg]
%   tSnow: Temperature of snow layer above the soil (in K)
%   rhoSnow: Density of the snow layer above the soil (in Kg m^-3) 
%   tSoil:  Temperature of the soil surface (in K)
%   mvSoil: Volumetric soil moisture (m3/m3)
%   sigSoil: Soil roughness (m);
%   s0: Vert and horizontal reflectivity at the snow-soil interface


%Constants
epiSnowMethod = 0; % Default = MEMLS, 1 = HUT 
epiSoilMethod = 1; % Default = HUT, 1 = Mironov 2010
sihMethod = 0;

mgSoil = 0.6; %Gravametric soil density for Mironov 2010
tSoil = tSoilK - 273.15;
tSnow = tSnowK - 273.15; %Temp of base snow layer
rhoSnow = rhoSnow./1000; %Density of base snow layer 

  if epiSoilMethod == 1
        epiSoil=real(tdgrmdm(tSoilK,mgSoil,mvSoil,freq)); %Method frmo Mironov 2010
   else
        %epiSoil=real(epss(mvSoil,tSoil,freq*1e9)) %Method from MEMLS (Dobson)
        %higheps(theta,f,sand,clay,rhob,rhos)
        [er,ei]=higheps(mvSoil,freq*1e9,0.4,0.3,1.6,2.6);
        epiSoil = er;
  end
    
  if epiSnowMethod == 1
        epiSnow = real(huteps(freq,tSnowK,rhoSnow,0,0)); % HUT method
  else
        epiSnow = real(epsr(rhoSnow)); % MEMLS method
  end
   
  theta_r = (theta * pi) / 180;
  
  
  if sihMethod == 1   
       nsSnow  = sqrt(epiSnow);
         nsSoil = sqrt(epiSoil);
         tei = [asin(sin(theta_r)./nsSoil);asin(sin(theta_r)./nsSnow)];
         epsEff = epiSoil./epiSnow;
         fresnel_h = gammah(epsEff,theta);
         fresnel_v = gammah(epsEff,theta); %Need to replace with HUT gammav
  else
         nsSnow  = sqrt(epiSnow);
         nsSoil = sqrt(epiSoil);
         tei = [asin(sin(theta_r)./nsSoil);asin(sin(theta_r)./nsSnow)];
         [sih,siv] = fresnelyc(tei,[epiSoil;epiSnow]);
         fresnel_h = sih;
         fresnel_v = siv;
         
       
  end
  
  theta_r_local = (tei(2) * pi) / 180; %Not sure if this should be the snow or ground theta.
  
  kSigma = real(2*pi*freq*1e9*sqrt(pi*4e-7*8.8542e-12*epiSnow))*sigSoil; % Soil roughness paramterization
   

  r_h_mod = fresnel_h.*exp(-kSigma.^(sqrt(0.1.*cos(theta_r))));
  r_v_mod = r_h_mod.*cos(theta_r).^0.655;%1.72;%0.655;%1.38;
  %r_v_mod = r_h_mod.*(0.635-(0.0014*(theta-60)))%1.72;%0.655;%1.38;