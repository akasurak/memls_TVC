function [r_h_mod, r_v_mod , fresnel_h, fresnel_v] = snowsoilreflectivity(freq, theta, tSnowK, rhoSnow, tSoilK, mvSoil, sigSoil)
%	 Reflectivity estimate for rough bare soil covered with snow
%	 Date: 13/05/14
%	 Author: J. King
%	 Based on code by B. Montpetit (ss_ref_grass.m, ruffsoil.m) and
%	 reflectivity approach detailed in Wegmuller & Matzler 2009
%	 Uses: ruffsoil.m, huteps.m, tdgrmdm.m. epss.m, epsr.m, gammah.m
%	 Date: 13/05/14
%
%	 Input
%	 freq: Centre frequency [in GHz]
%	 theta: Incident angle at snow surface [in Deg]
%	 tSnow: Temperature of snow layer above the soil (in K)
%	 rhoSnow: Density of the snow layer above the soil (in Kg m^-3) 
%	 tSoil: Temperature of the soil surface (in K)
%	 mvSoil: Volumetric soil moisture (m3/m3)
%	 sigSoil: Soil roughness (cm) = standard deviation of surface height (W&M 1999);



%Constants
epiSnowMethod = 0; % Default = MEMLS, 1 = HUT 
epiSoilMethod = 0; % Default = Dobson, 1 = Mironov 2010, 2 = Constant

mgSoil = 0.6; %Gravametric soil density for Mironov model
tSoil = tSoilK - 273.15;
tSnow = tSnowK - 273.15; %Temp of base snow layer
rhoSnow = rhoSnow./1000; %Density of base snow layer 
theta_r = (theta * pi) / 180;

C = 299792458;


if epiSoilMethod == 1
		epiSoil=real(tdgrmdm(tSoil,mgSoil,mvSoil,freq))^2; %Method frmo Mironov 2010
elseif epiSoilMethod == 2
		epiSoil = 6; 
else
		[er,ei]=higheps(mvSoil,freq*1e9,0.4,0.3,1.6,2.6);
		epiSoil = er;
end
 
if epiSnowMethod == 1
		epiSnow = real(huteps(freq,tSnowK,rhoSnow,0,0)); % HUT method %AK: find multilayer HUT eps in coreh2o's snowemis_nlayer.m
else
		epiSnow = real(epsr(rhoSnow)); % MEMLS method
end


nsSnow	= sqrt(epiSnow);
nsSoil = sqrt(epiSoil);
tei = [asin(sin(theta_r)./nsSoil);asin(sin(theta_r)./nsSnow)];
[sih,siv] = fresnelyc(tei,[epiSoil;epiSnow]);
fresnel_h = sih;
fresnel_v = siv;		




theta_r_local = (tei(2) * pi) / 180; %Not sure if this should be the snow or ground theta.

effWave = (C/nsSnow)/freq*1e-7; %Effective wavelength (cm) in snow
effK = 2*pi/effWave; %Effective wavenumber (rad/cm) in snow

kSigma = effK * sigSoil;


if kSigma < 0.07 || kSigma > 27.5
	disp(['snowsoilreflectivity>>	kSigma Out of range! (27.5 < [kSigma= ',num2str(kSigma),'] < 0.07); k = ' num2str(wavenumber) ' sigma = '	 num2str(sigSoil)]); 
end

%kSigma = real(2*pi*freq*1e9*sqrt(pi*4e-7*8.8542e-12*epiSnow))*sigSoil; %
%JK: The above is ben's implimentation of kSigma, according to the paper the valid
%kSigma range is 0.07 to 27.5 so I'm not sure its valid. 

r_h_mod = fresnel_h.*exp(-kSigma.^(sqrt(0.1.*cos(theta_r))));
r_v_mod = r_h_mod.*cos(theta_r).^0.655; %TODO Build in beta tuning function (0.655 currently)
