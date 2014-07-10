function epsi = epsr(roi)

%   calculates the dielectric permittivity for dry snow from 
%   density.
%
%   epsi = epsr(roi)
%       epsi:  real part of dielectric permittivity
%       roi:   density g/cm^3
%
%   Version history:
%      1.0    wi 15.7.95
%      1.1    wi 23.9.97 added Looyenga for snow denser than 0.4 g/cm^3
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

epsi = zeros(size(roi));
i = find(roi<=0.4);
j = find(roi>0.4);
vfi = roi./0.917;
ehb = 0.99913;
esb = 1.4759;
epsi(i) = 1 + 1.5995 .* roi(i) + 1.861 .* roi(i).^3;
epsi(j) = ((1 - vfi(j)) .* ehb + vfi(j) .* esb).^3; 
