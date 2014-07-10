function result = memlsmain(fGHz,tetad,s0h,s0v,ifile,Tsky,Tgnd,sccho)
%   Basic MEMLS program, computes the brightness temperatures Tv and Th
%   of a snowpack at a given frequency and incidence angle.
%     fGHz:    frequency [GHz]
%     tetad: incidence angle [deg]
%     s0h:   snow-ground reflectivity, h-pol 
%     s0v:   snow-ground reflectivity, v-pol
%     ifile: Snowpack-Input File with 
% layer-number, temp [K], volume fraction of liquid water, density [kg/m3],
% thickness [cm], Salinity (0 - 0.1) [ppt], expon.corr.length [mm] 
%     Tsky:  sky brightness temperature [K]
%     Tgnd:  ground temperature [K]
%     sccho: type of scattering coefficient (11 recommended)
%
%   Version history:
%      3.0    ma 02.04.2007
%      3.1    ma 05.04.2007 adapted to saline snow
%   Uses:
%   vK2eps, epswet, epsaliceimag, pfadi, pfadc, polmix, fresnelyc, 
%   slred, sccoeff, rt, layer
%
%   Copyright (c) 2007 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

%  Aux. input parameter for dielectric model of snow:
graintype=1;	% 1, 2 or 3 for different assumptions:
%"1" empirical snow measurements, "2" small spheres, "3" thin shells
teta = (tetad * pi) / 180; % Transformation from degrees to radians
c0=0.299793;               % vac. speed of light in m/ns
y=load(ifile)
di=y(:,5);
y1=find(di>0);
num  = y(y1,1);
Ti   = y(:,2);
Wi   = y(:,3);
roi  = y(:,4);
Sppt = y(:,6);
pci  = y(:,7)
N = length(num);
if N == 0 
   return    % test if there is there a layer at all
end
roi = roi./1000; % transforms density to units of g/cm3
di  = di./100;   % transforms thicknesses to units of m
cc=vK2eps(roi, graintype);     
v =cc(:,1);                  % ice volume fraction
kq=cc(:,2);                  % K^2 ratio
epsid=cc(:,3);               % epsilon of dry snow
nid=sqrt(epsid);             % real refract.index of dry snow
freq=fGHz;
   eii=epsaliceimag(freq,Ti,Sppt); % imag epsilon of saline ice
   epsiid=v.*eii.*kq.*nid;  % imag epsilon of dry snow component
   epsd=epsid+i*epsiid;
   eps  = epswet(freq,Ti,Wi,epsd);
   epsi =real(eps)         % real epsilon of snow, dry or wet
   epsii=imag(eps);         % imag epsilon of snow, dry or wet
   gai = (4*pi*freq).*imag(sqrt(eps))./c0;  %absorption coeff (1/m)
   ns  = sqrt(epsi);        % approx. real refractive index of snow
   tei = [asin(sin(teta)./ns);teta]
   dei = pfadi(tei,di);
   [sih,siv] = fresnelyc(tei,[epsi;1]);
[rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rkq]= slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,kq);
  [gbih,gbiv,gs6,ga2i] = sccoeff(rroi,rTi,rpci,freq,repsi,rgai,sccho,rkq);
  [rdei,rtei,tscat] = pfadc(teta,rdi,repsi,gs6);
  rsih = [s0h;rsih];
  rsiv = [s0v;rsiv];
  [rsih,rsiv] = polmix(tscat,rsih,rsiv);
  rtei.*180./pi;
  [ri,ti]  = rt(ga2i,gbih,rdei); 
  Dh   = layer(ri,rsih,ti,rTi,Tgnd,Tsky);
  N = length(rroi);
  Tbh = (1-rsih(N+1))*Dh(N) + rsih(N+1)*Tsky;
  [ri,ti]  = rt(ga2i,gbiv,rdei);
  Dv   = layer(ri,rsiv,ti,rTi,Tgnd,Tsky);
  Tbv = (1-rsiv(N+1))*Dv(N) + rsiv(N+1)*Tsky;
result=[Tbv,Tbh];
