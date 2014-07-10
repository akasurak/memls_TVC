function [ri,ti] = rt(gai,gbi,dei)

%   calculates the layer reflectivity and transmissivity   
%
%
%   [ri,ti] = rt(gai,gbi,dei)
%       ri:   layer reflectivity
%       ti:   layer transmissivity
%       gai:  absorption coefficient
%       gbi:  scattering coefficient
%       dei:  path length
%
%   Version history:
%      1.0    wi 15.7.95
%   
%   Uses: -
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

gamma = sqrt(gai .* (gai + 2 .* gbi));
t0i = exp(gamma .* dei .* (-1));
r0i = zeros(size(t0i));
i = find(gbi > 0.00001);
r0i(i) = gbi(i) ./ (gai(i) + gbi(i) + gamma(i));
t02 = t0i.^2;
r02 = r0i.^2;
ri = r0i .* (1 - t02) ./ (1 - t02 .* r02);
ti = t0i .* (1 - r02) ./ (1 - t02 .* r02);
