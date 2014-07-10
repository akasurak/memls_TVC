function [FH,FV] = fresnelrc(tei,epsr)

%   fresnel reflection coefficients (assuming eps'' = 0)
%     (layer n+1 is the air above the snowpack)
%
%       FH:   Fresnel reflection coefficient at h pol
%       FV:   Fresnel reflection coefficient at v pol
%       tei:  local incidence angle
%       epsr: (real part) dielectric permittivity
%
%   Version history:
%      1.0    wi 15.7.95
%   
%   Uses:
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

N = length(epsr);

FH = zeros(size(epsr(1:N-1)));
FV = zeros(size(epsr(1:N-1)));

% epsi = [epsi;0];


for n = 1:(N-1)
  epsn = epsr(n)/epsr(n+1);
  tein = tei(n+1);
  sinq = sin(tein)^2;
  qeps = sinq/epsn;
  wurz = sqrt(1-qeps);
  wsub = epsn-sinq;
  nd = sqrt(epsn);
  
  FH(n) = ((nd*wurz-cos(tein))/(nd*wurz+cos(tein)));
  FV(n) = ((wurz-nd*cos(tein))/(wurz+nd*cos(tein)));
   
end
