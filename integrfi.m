function integr = integrfi(xx,mui,muo,steps)

%   calculates integration  over incident directions, mui
%       from mini to maxi in 'steps' intervals of the fi
%       phase function.
%
%   integr = integrfi(xx,mui,muo,steps)
%
%   Version history:
%      1.0     wi 27.05.98 
%
%   Uses: -
%
%   Copyright (c) 1998 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

delta = pi / steps;
f0 = 0.5 * delta;
integr = 0;
x2 = 2.*xx.^2;

for ifi = 1:steps
  fi = f0+(ifi-1)*delta;
  sii = sqrt(1-mui.^2);
  sio = sqrt(1-muo.^2);
  cofi = cos(fi);
  si2fi = 1-cofi.^2;
  coste = mui .* muo + sii .* sio .* cofi;
  si2chi = 0.5 .* (1+coste.^2);
  func = si2chi./(1+(1-coste).*x2).^2;
  integr = integr + func .* delta;
end

integr = integr ./ pi;
 