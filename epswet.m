function result = epswet(f,Ti,Wi,epsd);

%   calculates complex dielectric constant of wet snow
%   using Maxwell-Garnett Mixing rule of water in dry snow
%   for prolate spheroidal water with experimentally determined
%   depolarisation factors.
%   Water temperature is at 273.15 K, with epsilon
%   of water from Liebe et al. 1991.
%       epsd:  complex epsilon of dry snow
%       f:   frequency [GHz]
%       Ti:  physical snow temperature [K]
%       Wi:  wetness [volume fraction]
%
%   Version history:
%      1.0    wi 15.7.95
%      2.0    ma 31.5.2005: Wi is volume fraction (not %) 
%      3.0    ma 2.4.2007 : adjustments, new function name 
%   Uses: epswater (since Version 3)
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland
%
Aa = 0.005;    % depolarisation factors of prolate 
Ab = 0.4975;   % water inclusion (Matzler 1987)
Ac = Ab;
ew=epswater(f, 273.15); 
Ka=epsd./(epsd+Aa*(ew-epsd));
Kb=epsd./(epsd+Ab*(ew-epsd));
K =(Ka+2*Kb)/3;
epsz=(1-Wi).*epsd+Wi.*ew.*K;
epsn=1-Wi.*(1-K);
eps=epsz./epsn;  % Maxwell-Garnett Mixing of water in dry snow
result=eps;