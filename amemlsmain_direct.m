%function result = amemlsmain(fGHz,tetad,s0h,s0v,ss0h,ss0v,ifile,Tgnd,sccho,m,q)
function result = amemlsmain(fGHz,tetad,s0h,s0v,ss0h,ss0v,y,Tsky,Tgnd,sccho,m,q)
%   Basic MEMLS program, computes the (brightness temperatures Tv and Th),
%   diffuse and specular reflectivities and backscattering coefficients
%   of a snowpack at a given frequency and incidence angle.
% Input parameters:
%     fGHz:    frequency [GHz]
%     tetad: incidence angle [deg]
%     s0h:   snow-ground reflectivity, h-pol 
%     s0v:   snow-ground reflectivity, v-pol
%     ss0h:  specular part of snow-ground reflectivity, h-pol 
%     ss0v:  specular part of snow-ground reflectivity, v-pol
%     ifile: Snowpack-Input File with 
% layer-number, temp [K], volume fraction of liquid water, density [kg/m3],
% thickness [cm], Salinity (0 - 0.1) [ppt], expon.corr.length [mm] 
%     (Tsky: sky brightness temperature [K], not used presently)
%     Tgnd:  ground temperature [K]
%     sccho: type of scattering coefficient (11 recommended)
%     m:     mean slope of surface undulations (typical 0.05 to 0.1)
%     q:     cross pol fraction (typical 0.1 to 0.3)
%   Version history:
%      3.0    ma 02.04.2007
%      3.1    ma 05.04.2007 adapted to saline snow
%      4.0    ma 2013 addition of backscatter coefficient 'MEMLS-Active'      
%   Uses:
%   vK2eps, epswet, epsaliceimag, pfadi, pfadc, polmix, fresnelyc, 
%   slred, sccoeff, rt, layer
%
%  Aux. input parameter for dielectric model of snow:
graintype=1;	% 1, 2 or 3 for different assumptions:
%"1" empirical snow measurements, "2" small spheres, "3" thin shells
forwa=4;  % emprical enhancement factor of extinction by diffraction
% alternative would be the use of the formula of Hallikainen & al 1996
teta = (tetad * pi) / 180; % Transformation from degrees to radians
c0=0.299793;               % vac. speed of light in m/ns

%y=load(ifile);
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
   epsi =real(eps);         % real epsilon of snow, dry or wet
   epsii=imag(eps);         % imag epsilon of snow, dry or wet
   gai = (4*pi*freq).*imag(sqrt(eps))./c0;  %absorption coeff (1/m)
   ns  = sqrt(epsi);        % approx. real refractive index of snow
   tei = [asin(sin(teta)./ns);teta];
   dei = pfadi(tei,di);
   [sih,siv] = fresnelyc(tei,[epsi;1]);
[rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rkq]= slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,kq);
  N = length(rnum);        % number of incoherent layers
  [gbih,gbiv,gs6,ga2i] = sccoeff(rroi,rTi,rpci,freq,repsi,rgai,sccho,rkq);
  [rdei,rtei,tscat] = pfadc(teta,rdi,repsi,gs6);
  
  % Comp of specular reflectivity rsh, rsv, rs0:
  rsih = [s0h;rsih];
  rsiv = [s0v;rsiv];
  rextc= rgai+forwa*gs6;  % extinct coeff for specular/coherent reflections
  Rh(1)=ss0h; Rv(1)=ss0v; % coherent refl at the ground surface
  u=exp(-rextc.*rdei); 
  u2=u.*u;
  for j=2:N+1
      Rh(j)=rsih(j)+Rh(j-1)*((1-rsih(j))*u(j-1))^2/(1-u2(j-1)*rsih(j)*Rh(j-1));
      Rv(j)=rsiv(j)+Rv(j-1)*((1-rsiv(j))*u(j-1))^2/(1-u2(j-1)*rsiv(j)*Rv(j-1));
  end;
  rsh=Rh(N+1);
  rsv=Rv(N+1);
  rs0=0.5*(rsh+rsv);
      
  [rsih,rsiv] = polmix(tscat,rsih,rsiv);
  rtei.*180./pi;
  
% The following 6 lines are to be used if Tbh and Tbv are required:   
  [ri,ti]= rt(ga2i,gbih,rdei); 
  Dh  = layer(ri,rsih,ti,rTi,Tgnd,Tsky);
  Tbh = (1-rsih(N+1))*Dh(N) + rsih(N+1)*Tsky;
  [ri,ti]= rt(ga2i,gbiv,rdei);
  Dv  = layer(ri,rsiv,ti,rTi,Tgnd,Tsky);
  Tbv = (1-rsiv(N+1))*Dv(N) + rsiv(N+1)*Tsky;

% Tb for Tsky = 0 K:
  [ri,ti]= rt(ga2i,gbih,rdei); 
  Dh0    = layer(ri,rsih,ti,rTi,Tgnd,0);
  Tbh0   = (1-rsih(N+1))*Dh0(N);
  [ri,ti]= rt(ga2i,gbiv,rdei);
  Dv0    = layer(ri,rsiv,ti,rTi,Tgnd,0);
  Tbv0   = (1-rsiv(N+1))*Dv0(N);

  % Tb for Tsky = 100 K:  
  [ri,ti]= rt(ga2i,gbih,rdei); 
  Dh100  = layer(ri,rsih,ti,rTi,Tgnd,100);
  Tbh100 = (1-rsih(N+1))*Dh100(N) + rsih(N+1)*100;
  [ri,ti]= rt(ga2i,gbiv,rdei);
  Dv100  = layer(ri,rsiv,ti,rTi,Tgnd,100);
  Tbv100 = (1-rsiv(N+1))*Dv100(N) + rsiv(N+1)*100;
  
  rv     = 0.01*(Tbv100 - Tbv0); % Total reflectivity v pol
  rh     = 0.01*(Tbh100 - Tbh0); % Total reflectivity h pol
  rdv    = rv-rsv;  % diffuse reflectivity v pol
  rdh    = rh-rsh;  % diffuse reflectivity h pol
   
mu=cos(teta);
mu2=mu.*mu;
m2=m*m;
sigma0dv =4*rdv.*mu2;
sigma0dh =4*rdh.*mu2;
sigma0dvv= (1-q)*sigma0dv;
sigma0dhh= (1-q)*sigma0dh;
sigma0hv = q*0.5*(sigma0dv+sigma0dh);
xpon     = -(tan(teta)).^2/(2*m2);
sigma0s  = rs0*exp(xpon)/(2*m2.*mu2.*mu2);
sigma0vv = sigma0dvv + sigma0s;
sigma0hh = sigma0dhh + sigma0s;
result.reflec = [rv,rh,rdv,rdh,rsv,rsh,rs0];
result.sigma0 = [sigma0vv,sigma0hh,sigma0hv];
result.Tb     = [Tbv,Tbh];