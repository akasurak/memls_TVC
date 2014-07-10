function [dei,tei,tscat] = pfadc(teta,di,epsi,gs6)

%   calculates the effective path length in a layer
%
%   [dei,tei,tscat] = pfadc(teta,di,epsi,gs6)
%       dei:  effective path length [m]
%       tei:  local incidence angle
%       tscat: scattering 
%       teta: incidence angle at snow air interface
%       di:   thickness [m]
%       epsi: dielectric permittivity
%       gs6:  6-flux scattering coefficient
%
%   Version history:
%      1.0    wi 15.10.97
%
%   
%   Uses: -
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

N = length(di);
ns = sqrt(epsi);
costetasn = sqrt(1-(sin(teta)./ns).^2);
cosc = sqrt(1-(1./ns).^2);
costetasc = 0.5 .* (1 + cosc);
dei = di./costetasn;

tauscat = zeros(size(epsi)+[1,0]);
tscat = zeros(size(epsi));
costeta = zeros(size(epsi));

for m = N:-1:1
   tauscat(m) = tauscat(m+1) + dei(m) * gs6(m)/2;
   tscat(m) = exp(-1 * tauscat(m));
   costeta(m) = tscat(m) * costetasn(m) + (1-tscat(m)) * costetasc(m);
end

tei = acos(costeta);
tei*180/pi;
costeta;
tauscat;
tscat;
