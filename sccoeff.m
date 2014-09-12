function [gbih,gbiv,gs6,ga2i] = sccoeff(roi,Ti,pci,freq,epsi,gai,sccho,kq)
%   calculates the scattering coefficient from structural parameters
%   different algorithms can be chosen, by changing "sccho"
%       gbih:  2-flux scattering coefficient at h pol
%       gbiv:  2-flux scattering coefficient at v pol
%       gs6:   6-flux scattering coefficient
%       ga2i:  2-flux absorption coefficient
%       roi:   density
%       Ti:    physical temperature
%       pci:   correlation length
%       freq:  frequency
%       epsi:  real part of relative dielectric constant of dry snow
%       gai:   absorption coefficient
%       sccho: scattering coefficient algorithm chosen
%       kq:    squared field ratio K^2 (for sccho=12)
%
%   Version history:
%   1.0b    wi 15.7.95
%   1.0     wi 23.9.97 bug fixed
%   1.1     wi 26.9.97 latest fit on experimental data was added (option 7)
%   1.2     wi 13.10.97 option 8 added, adapted scattering of a shell/sphere 
%   1.3     wi  4.11.97 option 9, 10 and 11 added 
%   1.4     wi 27.05.98 born approximation added (borna.m)
%   3.0     ma 03.04.2007 adapted to Version 3
%   4.0     ma May 2013   avoid con<0 for wahl=11
%   Uses:
%   borna, epsicereal (both for sccho = 12 only)
%
%   Copyright (c) 2007 by the Institute of Applied Physics, 
%   University of Bern, Switzerland
%   constants
c = 2.99793; 
roice = 0.917;
% specular components of vol. scattering coefficient set to 0: 
dgb0h = 0; dgb0v = 0;
% vacuum wave number
k0 = freq*(2*pi/0.299793);
eice = 3.18;
vfi = roi./roice;
wahl = sccho;             % selection of scattering algorithm
% 6-flux scattering coefficient
if wahl == 1
   gs6 = ((130 * ((freq/50)^2.7)) .* pci.^3) ./ (roi.^1.3 + 0.001);
end
% fit from 26.8.97 auf alle Daten v-pol, > 11 GHz
if wahl == 2
   gs6 = 0.0704 * (pci.^2.32).*(freq.^1.68).*roi.^(-0.63);
end
% for spheres
if wahl == 4
   epseff = (2-eice+3.*vfi.*(eice-1)+ sqrt((2-eice+3.*vfi.*(eice-1)).^2+8.*eice))./4;
   sphe = (3/32).*(0.001.*pci).^3.*k0.^4.*vfi.*(1-vfi).*abs((2.*epseff+1).*(eice-1)./(2.*epseff+eice)).^2;
   gs6 = sphe;
end
% for shells
if wahl == 5
   epseff = 1+(vfi.*(eice-1).*(2+1/eice))./(3-vfi.*(1-1/eice));
   shel = abs(2/3 + 1./(3.*eice.^2)).*(0.001.*pci).*k0.^2.*vfi.*(1-vfi).*(eice-1).^2./(16.*epseff);
   gs6 = shel;
end
% fit from 26.9.97
if wahl == 7
   gs6 = 73.21 * (pci.^3).*((freq./50).^2.68).*roi.^(-1);
end
% fit from 13.10.97 
if wahl == 8 %eqn 78 in memls3 doc
   gs6 = 136 .* (pci.^2.85) .* ((freq./50).^2.5) ./ (roi + 0.001);
end
% fit from 4.11.97 (without density) 
if wahl == 9
   gs6 = 564 .* (pci.^3.0) .* ((freq./50).^2.5);
end
% fit from 4.11.97 (without density, uses corr. length from exp. fit!)
if wahl == 10 %eqn 79 in memls3 doc
   gs6 = (3.16 .* pci + 295 .* (pci.^2.5)).* ((freq./50).^2.5);
end
% fit from 4.11.97 (with density, uses corr. length from exp. fit!)
if wahl == 11 %eqn 80 in memls3 doc
   con=9.20 .* pci - 1.23 .* roi + 0.54;
   gs6 = max(0,con).^2.5 .* ((freq./50).^2.5);
end

% Improved Born Approximation
if wahl == 12
  eice=epsicereal(Ti);
  [gb6,gc6,gf6,gs6] = borna(k0,vfi,pci,epsi,eice,kq); 
else
  omega = sqrt((epsi - 1)./epsi);
  gb6 = 0.5 .* gs6 .* (1-omega);
  gc6 = 0.25 .* gs6 .* omega;
end
gbiv = zeros(size(gs6));
gbih = zeros(size(gs6));
gtr = zeros(size(gs6));
ga2 = zeros(size(gai));

% -> 2 Flux Coefficients
gtr = (4 .* gc6) ./ (gai + 2 .* gc6);
ga2i = gai .* (1 + gtr);
gbih = (gb6 + dgb0h) + gtr .* gc6;
gbiv = (gb6 + dgb0v) + gtr .* gc6;
% result=[gbih,gbiv,gs6,ga2i];