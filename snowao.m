   function A = snowao(v)

%   computes the depolarization factor of oblate snow grains
%   Note 10, MŠtzler 1997
%
%   A = snowao(v)
%       A:    depolarization factor of oblate snow grains
%       v:    volume fraction of ice
%
%   Version history:
%      1.0    wi 29.5.98
%
%   Uses: -
%       
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

N = length(v);
A = zeros(N,1);
A = A + 0.3;
for i=1:N
   if v(i) < 0.55
      A(i) = 0.476 - 0.64 * v(i);
   end
   if v(i) <= 0.333
      A(i) = 0.1 + 0.5 * v(i);
   end
end

