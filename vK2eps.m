function result = vK2eps(roi,graintype)
%   Computes volume fraction, K^2 (the squared ratio between internal 
%   and external E-field) and real diel. constant in the
%   Effective-Medium Approximation (Polder and van Santen, 1946) of 
%   dry snow with given depolarization factors, A=A1=A2 and A3=1-2A 
%   all versus snow density, roi (graintype=1),  for spheres
%   (gaintype=2), or spherical shells (graintype=3).
%
%       kq:    squared E field ratio
%       roi:   dry snow density (g/cm3)
%       graintype: 1 (snow), 2 (spheres), 3 (thin sph. shells)
%
%   Version history:
%      1.0    wi 29.5.98, 
%      3.0    lb 2.3.2007, cm 2.4.2007
%
%   Uses: epsr 
%
%   Copyright (c) 2007 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

v = roi ./ 0.917;
N = length(roi);
eps2 = 3.185; eps1=1;

if graintype==1  % snow with the empirical formula
    A = zeros(N,1);
    A = A + 1/3;
    for i=1:N
       if v(i) < 0.55
          A(i) = 0.476 - 0.64 * v(i);
       end
        if v(i) <= 0.333
            A(i) = 0.1 + 0.5 * v(i);
        end
    end
    epseff=epsr(roi);
    A3 = (1 - 2 .* A);
    ea = epseff .* (1-A) + A;
    ea3 = epseff .* (1-A3) + A3;
    K1 = (ea  ./ (ea+A  .* (eps2-1))).^2;
    K3 = (ea3 ./(ea3+A3 .* (eps2-1))).^2;
    kq = (2 .* K1 + K3) ./ 3;
   
elseif graintype==2  %case of small spherical scatterers
    % Solve Polder v. Santen with A=1/3, eps1=1
    aa=2; bb=eps2-2-3*v*(eps2-1); cc=-eps2;
    epseff=(-bb+sqrt(bb.*bb-4*aa*cc))/2/aa; 
    ea=2*epseff/3+1/3;
    kq=(ea./(ea+(eps2-eps1)/3)).^2;
      
elseif graintype==3  %case of thin spherical shells
    kq = (2/3) + (eps1.^2)./(3*(eps2.^2));  
    kq = kq*ones(N,1);
    epseff=eps1+v*(eps2-eps1)*(2+eps1/eps2)./(3-v*(1-eps1/eps2));
end
result=[v, kq, epseff];