 function result = ars0(fGHz,ss0,forwa, ifile,sccho)
%   MEMLS program for specular reflectivity at vertical incidence
%     fGHz:    frequency [GHz]
%     s0:   snow-ground reflectivity
%     ifile: Snowpack-Input File with 
% layer-number, temp [K], volume fraction of liquid water, density [kg/m3],
% thickness [cm], Salinity (0 - 0.1) [ppt], expon.corr.length [mm] 
%     sccho: type of scattering coefficient (11 recommended)
%
%   Version history:
%     ma 2013 for active version      
%   Uses:
%   vK2eps, epswet, epsaliceimag, pfadi, fresnelyc, 
%   slred, sccoeff,
%
%  Aux. input parameter for dielectric model of snow:
graintype=1;	% 1, 2 or 3 for different assumptions:
%"1" empirical snow measurements, "2" small spheres, "3" thin shells
teta = 0; % Incidence angle
c0=0.299793;               % vac. speed of light in m/ns
y=load(ifile);
di=y(:,5);
y1=find(di>0);
num  = y(y1,1);
Ti   = y(:,2);
Wi   = y(:,3);
roi  = y(:,4);
Sppt = y(:,6);
pci  = y(:,7);
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
   tei = [ns*0;0];          % vertical incidence
   dei = di;                % path length = layer thickness
   [sih,siv] = fresnelyc(tei,[epsi;1]);
[rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rkq]= slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,kq);
  N = length(rroi)        % number of incoherent layers
  [gbih,gbiv,gs6,ga2i] = sccoeff(rroi,rTi,rpci,freq,repsi,rgai,sccho,rkq)
  rextc= rgai+forwa*gs6; % extinct coeff for specular/coherent reflections
  % for rextc: alternatively we could use formula of Hallikainen & al 1996
  rsi0 = [ss0;rsih];
  R0(1) = ss0; % Specular refl at ground surface
  u=exp(-rextc.*rdei); 
  u2=u.*u;
  for j=2:N+1
      R0(j)=rsi0(j)+R0(j-1)*((1-rsi0(j))*u(j-1))^2/(1-u2(j-1)*rsi0(j)*R0(j-1));
  end;
  rs0=R0(N+1);
result=rs0;