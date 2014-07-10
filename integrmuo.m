function integr = integrmuo(xx,mui,mino,maxo,steps)

%   calculates integration  over incident directions, mui
%       from mino to maxo in steps intervals of the fi
%       phase function.
%   
%
%   integr = integrmuo(xx,mui,mino,maxo,steps)
%
%   Version history:
%      1.0     wi 27.05.98 
%
%   Uses: integrfi
%     
%
%   Copyright (c) 1998 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

dmu = maxo - mino;
delta = dmu ./ steps;
f0 = 0.5 .* delta;
integr = 0;

for imu = 1:steps
  muo = mino + f0 + (imu -1) .* delta;
  func = integrfi(xx,mui,muo,steps);
  integr = integr + func .* delta;
end

integr = integr ./ 2;
 