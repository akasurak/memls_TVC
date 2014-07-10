function [gb6,gc6,gf6,gs6] = borna(k0,vfi,pci,epsi,eice,kq)

%   calculates the scattering coefficient using Born Approximation
%   [gb6,gc6,gf6,gs6] = borna(k0,vfi,pci,epsi,eice,kq)
%       gb6: 6-flux back scattering coefficient 
%       gc6: 6-flux cross scattering coefficient 
%       gf6: 6-flux forward scattering coefficient 
%       gs6: 6-flux scattering coefficient 
%       k0:  vacuum wave number (1/m)
%       vfi: volume fraction of ice
%       pci:  correlation length  (mm)
%       epsi: real dielectric constant of snow
%       eice: real dielectric constant of ice
%       kq :  squared E-field ratio (from bornsnk)
%
%   Version history:
%      1.0     wi 27.05.98, 
%      3.0     lb 2.3.2007, 3.1 cm 2.4.2007 
%
%   Uses: integrmui
%
%   Copyright (c) 2007 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

steps = 11;                 % number if integration steps:
pci = pci .* 0.001;         % conversion from mm to m
muc = sqrt((epsi-1)./epsi); % critical cos 
% Eq. (1) of Ma+Wi 1999, but without denominator of (4):
aa = 2.*(pci.*k0).^3.*k0.*vfi.*(1-vfi).*(eice-1).^2.*kq;
xx = pci .* k0 .* sqrt(epsi);
% triple integration
% backward scattering
bb = integrmui(xx,muc,1,-1,-1.*muc,steps);
% transverse scattering
bt = integrmui(xx,muc,1,-1.*muc,muc,steps);
% forward scattering
bf = integrmui(xx,muc,1,muc,1,steps);
btot = bb + bt + bf;
gb6 = aa .* bb;
gc6 = 0.25 .* aa .* bt;
gf6 = aa .* bf;
gs6 = aa .* btot;