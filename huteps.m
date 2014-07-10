% snow_eps.m
% This function calculates the complex permittivity of snow cover (single layer only)
% Function is used for single layer HUT snow emission model
% JP/5.12.1998
%
% Input parameters:
% f = frequency (GHz)
% Tlumi = snow temperature (degrees C)
% pdsw = snow density (g/cm^3)
% Mv = snow moisture (%)
% SS = salinity (ppm) (Initial to be used typically: 0 ppm)

function [epsr] = huteps(f,Tlumi,pdsw,Mv,SS)

T0 = 273.15;
pice = 0.916; % Density of ice
im = sqrt(-1);
mv = Mv/100;

% Density of dry snow
pds = (pdsw - mv)/(1.0-mv);

% Real part of epsilon for dry snow
reds = 1+1.58*pds/(1-0.365*pds);

% Real part of epsilon for ice
rei = 3.1884 + 9.1e-4*(Tlumi+T0-273);  % Matzler & Wegmuller 1987

% Calculation of imaginary part of epsilon for ice
A=0.0026;               % impure ice -5 C (Matzler)
B=0.00023;
C=0.87;
iei_S=A/f+B*f^C;

A_p=6e-4;               % pure ice -5 C (Matzler)
B_p=6.5e-5;
C_p=1.07;
iei_P=A_p/f+B_p*f^C_p;
delta_iei = iei_S - iei_P; %difference between pure and impure ice at -5 C
% Renewed formula for pure ice (Hufford 1991):
invT = (300/(Tlumi+T0))-1;
alf = (0.00504 + 0.0062*invT) * exp(-22.1*invT);
% Modification (Mishima, Matzler):
B_1 = 0.0207;
B_2 = 1.16e-11;
bb = 335;
bet_M = (B_1/(Tlumi+T0)) * (exp(bb/(Tlumi+T0)))/((exp(bb/(Tlumi+T0))-1)^2) + B_2*f^2;
bet_delta = exp(-10.02+0.0364*(Tlumi+T0-273));
bet = bet_M + bet_delta;
iei = alf/f + bet*f; % improved imaginary part for pure ice 
% Effect of salinity:
iei = iei + delta_iei*SS/13; %imaginary part for impure ice

% Imaginary part of epsilon for dry snow (PVS-model)
vi=pds/pice;
ieds = 3*vi*iei*reds^2*(2*reds+1)/((rei+2*reds)*(rei+2*reds^2));

% Epsilon of wet snow
% (according to Matzler)
if mv>0
  Aa = 0.005;
  Ab = 0.4975;
  Ac = 0.4975;
  epssw = 88;
  epsinf = 4.9;
  ff_0= 9;
  f_0a = ff_0 * (1 + Aa*(epssw-epsinf)/(reds+Aa*(epsinf-reds)) );
  f_0b = ff_0 * (1 + Ab*(epssw-epsinf)/(reds+Ab*(epsinf-reds)) );
  f_0c = ff_0 * (1 + Ac*(epssw-epsinf)/(reds+Ac*(epsinf-reds)) );
  epsinf_a = (mv/3) * (epsinf-reds)/(1+Aa*((epsinf/reds)-1));
  epsinf_b = (mv/3) * (epsinf-reds)/(1+Ab*((epsinf/reds)-1));
  epsinf_c = (mv/3) * (epsinf-reds)/(1+Ac*((epsinf/reds)-1));
  epss_a = (mv/3) * (epssw-reds)/(1+Aa*((epssw/reds)-1));
  epss_b = (mv/3) * (epssw-reds)/(1+Ab*((epssw/reds)-1));
  epss_c = (mv/3) * (epssw-reds)/(1+Ac*((epssw/reds)-1));
  eps_a = epsinf_a + (epss_a-epsinf_a)/(1+im*f/f_0a);
  eps_b = epsinf_b + (epss_b-epsinf_b)/(1+im*f/f_0b);
  eps_c = epsinf_c + (epss_c-epsinf_c)/(1+im*f/f_0c);
  epsmarka = eps_a + eps_b + eps_c + (reds-im*ieds);
  rews = real(epsmarka);
  iews = -1*imag(epsmarka);
else
  rews = reds;
  iews = ieds;
end
epsr = rews-im*iews;
return
