function A = snowap(v)

%   computes the depolarization factor of prolate snow grains
%   Note 10, Mätzler 1997
%
%   A = snowap(v)
%       A:    depolarization factor of prolate snow grains
%       v:    volume fraction of ice
%
%   Version history:
%      1.0    wi 29.5.98
%
%   Uses: -
%       
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

N = length(v);
A = zeros(N,1);
A = A + 0.496;
for i=1:N
   if v(i) < 0.51
      A(i) = 0.18 - 0.62 * v(i);
   end
   if v(i) <= 0.333
      A(i) = 0.531 - 0.44 * v(i);
   end
   if v(i) < 0.071
      A(i) = 0.5;
   end
end

