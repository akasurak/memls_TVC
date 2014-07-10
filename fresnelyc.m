function [sih,siv] = fresnelyc(tei,epsi)

%   fresnel reflection coefficients (assuming eps'' = 0)
%     (layer n+1 is the air above the snowpack)
%
%       sih:  interface reflectivity at h pol
%       siv:  interface reflectivity at v pol
%       tei:  local incidence angle
%       epsi: real part of dielectric permittivity
%
%   Version history:
%      1.0    wi 15.7.97
%   
%   Uses:
%       epsr
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

N = length(epsi)-1;
siv = zeros(size(epsi(1:N)));
sih = zeros(size(epsi(1:N)));

for n = 1:N
  epso = epsi(n+1);
  epsu = epsi(n);
  tein = tei(n+1);
  sih(n) = ((sqrt(epso)*cos(tein) - sqrt(epsu - epso * sin(tein)^2))/(sqrt(epso)*cos(tein) + sqrt(epsu - epso * sin(tein)^2)))^2;
  siv(n) = ((epsu*cos(tein) - sqrt(epso)*sqrt(epsu - epso*sin(tein)^2))/(epsu*cos(tein) + sqrt(epso)*sqrt(epsu - epso*sin(tein)^2)))^2;
end

