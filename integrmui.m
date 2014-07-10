function integr = integrmui(xx,mini,maxi,mino,maxo,steps)

%   calculates integration  over incident directions, mui
%       from mini to maxi in steps intervals of the fi and
%       muo integrated phase function.
%
%   integr = integrmui(xx,mini,maxi,mino,maxo,steps)
%
%   Version history:
%      1.0     wi 27.05.98 
%
%   Uses: integrmuo
%
%   Copyright (c) 1998 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

% constants
steps = 11;

dmu = maxi - mini;
delta = dmu ./ steps;
f0 = 0.5 .* delta;
integr = 0;

for imu = 1:steps
  mui = mini + f0 + (imu -1) .* delta;
  func = integrmuo(xx,mui,mino,maxo,steps);
  integr = integr + func .* delta;
end

integr = integr ./ dmu;
 