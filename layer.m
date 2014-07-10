function D = layer(ri,si,ti,Ti,Tgnd,Tsky)

%   calculates the upwelling brightness temperatures D (see Note 6
%   or Wiesmann and Matzler 1999)  
%
%   D = layer(ri,si,ti,Ti,Tgnd,Tsky)
%       D:    upwelling brightness temperature
%       ri:   layer reflectivity
%       si:   interface reflectivity
%       ti:   layer transmissivity
%       Ti:   physical temperature [K]
%       Tgnd: brightness temperature of the soil below the snowpack
%       Tsky: brightness temperature of the sky
%
%   Version history:
%      1.0    wi 15.7.95
%      1.1    wi 26.9.97  handles also the case of a single layer
%      1.2    wi 02.03.99 fixed error in 1 layer handling 
%  
%   Uses: -
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

N = length(ri);
ei = 1 - ri - ti;
if N < 1 
   disp('ERROR: No scattering layer')
   return
end
if N == 1 
%   k1 = (ri(1)*(1-si(1))*Tgnd+ti(1)*(1-si(2))*Tsky+ei(1)*Ti(1))/(1-ri(1)*si(1));
%   k2 = ti(1)*si(1)/(1-ri(1)*si(1));
%   k3 = 1-ri(1)*si(2)-ti(1)*si(1)*k2;
%   D = (ti(1)*si(1)*k1+ti(1)*(1-si(1))*Tgnd+ri(1)*(1-si(2))*Tsky+ei(1)*Ti(1))/k3;

    k1 = (1-ri(1)*si(1))*(1-ri(1)*si(2)) - ti(1)*si(1)*ti(1)*si(2);
    D = ti(1)*si(1) * ((1-si(1))*ri(1)*Tgnd + (1-si(2))*Tsky*ti(1) + ei(1)*Ti(1)) / (k1) + ...
        (1-ri(1)*si(1)) * ((1-si(1))*Tgnd*ti(1) + (1-si(2))*Tsky*ri(1) + ei(1)*Ti(1)) / (k1);
   
else
   
   M1 = diag(ri.*si(1:N));
   H = diag(ti(1:N-1).*(1-si(2:N)));
   M1(1:N-1,2:N) = M1(1:N-1,2:N) + H;

   I = eye(size(M1));

   M2 = diag(ti.*si(2:N+1));
   H = diag(ri(2:N).*(1-si(2:N)));
   M2(2:N,1:N-1) = M2(2:N,1:N-1) + H;

   M3 = diag(ti.*si(1:N));
   H = diag(ri(1:N-1).*(1-si(2:N)));
   M3(1:N-1,2:N) = M3(1:N-1,2:N) + H;

   M4 = diag(ri.*si(2:N+1));
   H = diag(ti(2:N).*(1-si(2:N)));
   M4(2:N,1:N-1) = M4(2:N,1:N-1) + H;

   E = ei .* Ti;
   E(1) = E(1) + ri(1) .* (1 - si(1)) .* Tgnd;
   E(N) = E(N) + ti(N) .* (1 - si(N+1)) .* Tsky;
   E;

   F = ei .* Ti;
   F(1) = F(1) + ti(1) .* (1 - si(1)) .* Tgnd;
   F(N) = F(N) + ri(N) .* (1 - si(N+1)) .* Tsky;
   F;

   M5 = M3 * (inv(I - M1) * M2) + M4;

   D = inv(I - M5) * (M3 * inv(I - M1) * E + F);
end

