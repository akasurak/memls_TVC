function  [num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,Wi,gai,kq] = slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,kq) 

%   locates and treats coherent layers in a snowpack
%     see Technote 11
%     with repo=1 the snow layer table is printed after substeps
%
%       num:  index of the layer in the original snowpack
%       roi:  density [g/cm^3]
%       epsi:
%       epsii:
%       tei:  local incidence angle
%       sih:  layer reflectivity at h pol
%       siv:  layer reflectivity at v pol
%       di:   layer thickness
%       dei:  local path length [m]
%       Ti:   physical snow temperature [K]
%       pci:  correlation length [mm]
%       Wi:   wetness  
%       gai:  absorption coefficient
%       freq: frequency
%
%      1.0    wi 21.8.95
%      2.0    wi 13.8.98 completely rewritten
%      3.0    ma 3.04.2007 kq added
%   
%   Uses:
%       fresnelrc
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

%  constants
cc = 0.299793;
FIC = 4 * pi * freq / cc;
fc = 4.712;
repo = 0;
%   is there any layer -> checked in main

N = length(roi);
theta = tei(N+1);
ns = sqrt(epsi);
fi = FIC .* di .* ns .* cos(tei(1:N));


if repo == 1
   disp ([num,di,roi,Ti,pci,ns,gai,kq,epsi,epsii,fi,tei(1:N).*180./pi]);
end
%   find thin layers 
i = find(fi < fc);
A = 0.* roi;
A(i) = A(i) + 1;
%   bottom layer is assumed noncoherent
A(1) = 0;

if length(i) > 0
  disp ([num2str(length(i)),' coherent layers of ',num2str(N),' detected: ',num2str(freq),' GHz']);

  %   identify succeeding coherent layers and mark the packs from 2 ... scmax
  pl = 0;
  sc = 1;
  scmax = 0;
  ml = 0;
  mlo = 0;
  for m = 2:N
    if A(m) == 1 & pl == 1
      if ml == 0
         sc = sc + 1;
         ml = 1;
         A(m-1) = sc;
      end
      % A(m) = 1, pl = 1   
      A(m) = sc;
      scmax = sc;
    else 
      if pl == 1;
      % A(m) = 0, pl = 1 -> non coherent   
        pl = 0;
      else  
        if A(m) == 1 
          % A(m) = 1, pl = 0 -> first coherent
          pl = 1;
          ml = 0;
        end
      end
    end
  end

  %  combine succeeding coherent layers by weighting with the phase
  if scmax > 0
    for m = 2:scmax
      B = find(A == m);
      fitot = sum(fi(B));
      fitv = fi(B) ./ fitot;
      tal = max(B);
      di(tal)   = sum(di(B));
      dei(tal)  = sum(dei(B));
      roi(tal)  = sum(roi(B).* fitv);
      Ti(tal)   = sum(Ti(B).* fitv);
      Wi(tal)   = sum(Wi(B).* fitv);
      pci(tal)  = sum(pci(B).* fitv);
      ns(tal)   = sum(ns(B).* fitv);
      gai(tal)  = sum(gai(B).* fitv);
      kq(tal)   = sum(kq(B).* fitv);
      epsi(tal) = sum(epsi(B).* fitv);
      epsii(tal)= sum(epsii(B).* fitv);
      tei(tal)  = sum(tei(B).* fitv);
      fi(tal) = fitot;
      A(tal) = 1;
    end
  end     
  i = find(A < 2);
  if length(i) > 0
    num=num(i);roi=roi(i);tei=tei([i;N+1]);di=di(i);Ti=Ti(i);pci=pci(i);Wi=Wi(i);gai=gai(i);
    ns=ns(i);epsi=epsi(i);epsii=epsii(i);fi=fi(i);dei=dei(i); kq=kq(i);
  end  
  N = length(roi);
  if (repo == 1)
    disp ([num,di,roi,Ti,pci,ns,gai,kq,epsi,epsii,fi,tei(1:N).*180./pi]);
  end

  i = find(fi < fc);
  A = roi .* 0;
  A(i) = A(i) + 1;
  A(1) = 0;
  sih = zeros(size(roi));
  siv = zeros(size(roi));
  X   = zeros(size(roi));
     
  %   calculate interface reflection coefficients
  [FH,FV] = fresnelrc(tei,[epsi;1]);

  %   reduction on layers of type 0 (coherent layer effects are 
  %     taken into account in the layer reflectivities)
  %     for layers of type 0 shi = FH^2
  i = find(A == 0);
  sih(i) = FH(i).^2;
  siv(i) = FV(i).^2;


  %     for layers of type 1 shi-1 = ...
  i = find(A == 1);
  X(i) =  2.*FH(i).*FH(i-1).*cos(fi(i));
  sih(i-1) = (FH(i).^2+FH(i-1).^2+X(i))./(1+FH(i).^2.*FH(i-1).^2+X(i)); 
  X(i) =  2.*FV(i).*FV(i-1).*cos(fi(i));
  siv(i-1) = (FV(i).^2+FV(i-1).^2+X(i))./(1+FV(i).^2.*FV(i-1).^2+X(i)); 
  N = length(di);
     
  if (repo == 1)
    disp ([num,di,roi,Ti,pci,ns,gai,kq,epsi,epsii,fi,tei(1:N).*180./pi,FH,FV,sih,siv]);
  end
   
  %     remove layers of type 1
  i = find(A == 0);
  ns=ns(i);num=num(i);roi=roi(i);tei=tei(i);di=di(i);Ti=Ti(i);pci=pci(i);Wi=Wi(i);gai=gai(i);kq=kq(i);
  FH=FH(i);FV=FV(i);sih=sih(i);siv=siv(i);fi=fi(i);epsi=epsi(i);epsii=epsii(i);dei=dei(i);
  N = length(di);
  tei = [tei;theta];
  if (repo == 1)
    disp ([num,di,roi,Ti,pci,ns,gai,kq,epsi(1:N),epsii,fi,tei(1:N).*180./pi,FH,FV,sih,siv])
  end
   
else
   disp (['no coherent layers detected: ',num2str(freq),' GHz']);
   if (repo == 1)
      disp ([num,di,roi,Ti,pci,ns,gai,kq,epsi,epsii,fi,tei(1:N).*180./pi,sih,siv]);
   end
end
