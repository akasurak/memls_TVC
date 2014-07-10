function dei = pfadi(tei,di)

%   calculates the effective path length in a layer
%
%   dei = pfadi(tei,di)
%       dei:  effective path length [m]
%       tei:  local incidence angle
%       di:   thickness [m]
%
%   Version history:
%      1.0    wi 15.7.95
%   
%   Uses:
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

N = length(di);
dei = di./cos(tei(1:N));

