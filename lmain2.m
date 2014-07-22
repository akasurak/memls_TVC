function result = lmain2(f1,f2,tetad,s0h,s0v,ifile,Tgnd,sccho)
%   MEMLS main program, plots and lists the frequencies(GHz), emissivities 
%   eh(%), ev(%), the emitted brightness temperatures Tebh(K), Tebv(K), 
%   and the effective physical temperatures Teffh(K), Teffv(K) 
%   of the emitting surface (ground and snowpack) of a layered snowpack. 
%     f1:    start frequency [GHz]
%     f2:    stop frequency [GHz]
%     tetad: incidence angle [deg]
%     s0h:   snow-ground reflectivity 
%     s0v:   snow-ground reflectivity
%     ifile: Snowpack-Input File with
%  layer-number, temp [K], vol. fraction of liquid water, density [kg/m3],
%  thickness [cm], Salinity (0 - 0.1) [ppt], expon.corr.length [mm] 
%     Tgnd:  ground temperature [K]
%     sccho: type of scattering coefficient (11 recommended)
%
%   Version history:
%      3.0   ma 04.04.2007 adapted from lmain
%      3.1   ma 05.04.2007 adapted to saline snow
%   Uses:
%     vK2eps, epswet, epsaliceimag, pfadi, pfadc, polmix, fresnelyc, 
%     slred, sccoeff, rt, layer
%
%   Copyright (c) 2007 by the Institute of Applied Physics, 
%   University of Bern, Switzerland
%
%  Aux. input parameter:
graintype=1;	% 1, 2 or 3 for different assumptions:
%"1" empirical snow measurements, "2" spheres, "3" thin shells
teta = (tetad * pi) / 180; % Transformation from deg. to radian
c0=0.299793;               % vac. speed of light in m/ns
y=load(ifile)
di=y(:,5);
y1=find(di>0);
num  = y(y1,1);
Ti   = y(:,2);
Wi   = y(:,3);
roi  = y(:,4);
Sppt = y(:,6);
pci  = y(:,7);
N = length(num);
if N == 0 
   return        % tests the presence of snow layers
end
roi = roi./1000; % transforms density to units of g/cm3
di  = di./100;   % transforms thicknesses to units of m
cc=vK2eps(roi, graintype);     
v =cc(:,1);                  % ice volume fraction
kq=cc(:,2);                  % K^2 ratio
epsid=cc(:,3);               % epsilon of dry snow
nid=sqrt(epsid);             % real refract.index of dry snow
% Plot variables
x  = f1:f2; yh = f1:f2; yv = f1:f2; yeh = f1:f2; yev = f1:f2;
%  loop over frequency
for freq = f1:f2
   eii=epsaliceimag(freq,Ti,Sppt); % imag epsilon of saline ice
   epsiid=v.*eii.*kq.*nid;  % imag epsilon of dry snow component
   epsd=epsid+i*epsiid;
   eps  = epswet(freq,Ti,Wi,epsd);
   epsi =real(eps);         % real epsilon of snow, dry or wet
   epsii=imag(eps);         % imag epsilon of snow, dry or wet
   gai = (4*pi*freq).*imag(sqrt(eps))./c0;  %absorption coeff (1/m)
   ns  = sqrt(epsi);        % approx. real refractive index of snow
   tei = [asin(sin(teta)./ns);teta];
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
% emitted brightness temperature (Tsky=0): 
  Dh   = layer(ri,rsih,ti,rTi,Tgnd,0);
  N = length(rroi);
  Tbeh = (1-rsih(N+1))*Dh(N);
  yeh(freq-f1+1) = Tbeh; 
  [ri,ti]  = rt(ga2i,gbiv,rdei);
  Dv   = layer(ri,rsiv,ti,rTi,Tgnd,0);
  Tbev = (1-rsiv(N+1))*Dv(N);
  yev(freq-f1+1) = Tbev;
% Tsky=100K:
  Dh   = layer(ri,rsih,ti,rTi,Tgnd,100);
  N = length(rroi);
  Tbh = (1-rsih(N+1))*Dh(N) + rsih(N+1)*100;
  yh(freq-f1+1) = Tbh; 
  [ri,ti]  = rt(ga2i,gbiv,rdei);
  Dv   = layer(ri,rsiv,ti,rTi,Tgnd,100);
  Tbv = (1-rsiv(N+1))*Dv(N) + rsiv(N+1)*100;
  yv(freq-f1+1) = Tbv;
end
yv=1-(yv-yev)/100;  % snowpack&ground system emissivity
yh=1-(yh-yeh)/100;

%  graphical output
figure;
h = semilogx(x,yh,'k--');
hold on
k = semilogx(x,yv,'r-');
legend([h,k],'eh','ev')
xlabel('Frequency [GHz]')
ylabel('Emissivity')
title(['snow layer model - inputdata: ',ifile])
disp(['T ground: ',num2str(Tgnd),'K;  snow-ground reflectivity (h/v): ',num2str(s0h),'/',num2str(s0v)])
%gtext(['T ground: ',num2str(Tgnd),'K;  snow-ground reflectivity (h/v): ',num2str(s0h),'/',num2str(s0v)])
text(x(2),yh(length(yh)),['T ground: ',num2str(Tgnd),'K;  snow-ground reflectivity (h/v): ',num2str(s0h),'/',num2str(s0v)])
hold off

figure;
semilogx(x,yeh,'k--',x,yev,'r-',x,yeh./yh,'k-',x,yev./yv,'r.');
hold on
legend('Teh','Tev','Teff,h','Teff,v')
xlabel('Frequency [GHz]')
ylabel('Emitted Brightness Temp,  Effective Temp [K]')
title(['snow layer model - inputdata: ',ifile])
%gtext(['T ground: ',num2str(Tgnd),'K;  snow-ground reflectivity (h/v): ',num2str(s0h),'/',num2str(s0v)])
text(x(2),yeh(length(yeh)),['T ground: ',num2str(Tgnd),'K;  snow-ground reflectivity (h/v): ',num2str(s0h),'/',num2str(s0v)])
hold off
result=[x',yh',yv',yeh',yev',yeh'./yh',yev'./yv'];
% Frequency(GHz), eh(%), ev(%), Tebh(K), Tebv(K), Teffh(K), Teffv(K) 