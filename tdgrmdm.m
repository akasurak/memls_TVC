function n_s=tdgrmdm(T,rd,mg,freq)
%n_s  is a complex refractive index n_s = n+ik
%n-refractive index, k-attenuation constant
%T-soil temperatur in Celsius
%rd-dry density of soil in g/cm^3
%mg-gravimetric moisture of soil g/g
%freq-frequency of EM wave in GHz


%T=T-6;
Tf=-6; % Soil freezing temperature
e0=8.85; % Vacuum permittivity
Tk=273.15+T;

% If soil is frozen
if (T<Tf)


e0b=8.19; %T<Tf
e0t=24; %T<Tf
e0iw=5.57; %T<Tf

betta0b=-4.8e-3;  %T<Tf
betta0t=-5.62e-3;  %T<Tf
betta0iw=-5.93e-3;  %T<Tf

T0b=-30; %T<Tf
T0t=-20; %T<Tf
T0iw=-20; %T<Tf

eoob=5.92; %T<Tf
eoot=5.05; %T<Tf
eooiw=4.07; %T<Tf

bettaoob=-4.37e-3;  %T<Tf
bettaoot=0;  %T<Tf
bettaooiw=-2.4e-3;  %T<Tf

Toob=-20; %T<Tf
Toot=-20; %T<Tf
Tooiw=-20; %T<Tf

dHb=1540; %T<Tf
dHt=1493; %T<Tf
dHiw=1352; %T<Tf

dSb=0.1; %T<Tf
dSt=0.45; %T<Tf
dSiw=0.39; %T<Tf

sigb=0.07; %T<Tf
sigt=1.6; %T<Tf
sigiw=0.077; %T<Tf

bettasigb=3.37e-3;  %T<Tf
bettasigt=44.7e-3;  %T<Tf
bettasigiw=0.69e-3;  %T<Tf

Tsigb=-20; %T<Tf
Tsigt=-15; %T<Tf
Tsigiw=-20; %T<Tf

tau_b=48*exp(dHb./Tk-dSb)./Tk;
tau_t=48*exp(dHt./Tk-dSt)./Tk;
%tau_t=24.5284-0.77777*T+0.0112*T^2;
tau_iw=48*exp(dHiw./Tk-dSiw)./Tk;
    
sig_b=sigb+bettasigb*(T-Tsigb);
%sig_b=0.14402+0.00382*T;
sig_t=sigt+bettasigt*(T-Tsigt);
sig_iw=sigiw+bettasigiw*(T-Tsigiw);

F0b=log((e0b-1.)/(e0b+2.));
F0t=log((e0t-1.)/(e0t+2.));
F0iw=log((e0iw-1.)/(e0iw+2.));

e0_b=(1+2*exp(F0b-betta0b*(T-T0b)))./(1-exp(F0b-betta0b*(T-T0b)));
%e0_b=13.6695+0.22604*T+0.00166*T^2;
e0_t=(1+2*exp(F0t-betta0t*(T-T0t)))./(1-exp(F0t-betta0t*(T-T0t)));
e0_iw=(1+2*exp(F0iw-betta0iw*(T-T0iw)))./(1-exp(F0iw-betta0iw*(T-T0iw)));

Foob=log((eoob-1.)/(eoob+2.));
Foot=log((eoot-1.)/(eoot+2.));
Fooiw=log((eooiw-1.)/(eooiw+2.));

eoo_b=(1+2*exp(Foob-bettaoob*(T-Toob)))./(1-exp(Foob-bettaoob*(T-Toob)));
eoo_t=(1+2*exp(Foot-bettaoot*(T-Toot)))./(1-exp(Foot-bettaoot*(T-Toot)));
eoo_iw=(1+2*exp(Fooiw-bettaooiw*(T-Tooiw)))./(1-exp(Fooiw-bettaooiw*(T-Tooiw)));

e_b=eoo_b+(e0_b-eoo_b)./(1-2*i*pi*freq*tau_b*1.e-3)+i*sig_b/(2*pi*freq*e0*1.d-3);
e_t=eoo_t+(e0_t-eoo_t)./(1-2*i*pi*freq*tau_t*1.e-3)+i*sig_t/(2*pi*freq*e0*1.d-3);
e_iw=eoo_iw+(e0_iw-eoo_iw)./(1-2*i*pi*freq*tau_iw*1.e-3)+i*sig_iw/(2*pi*freq*e0*1.d-3);

riw=0.917; 

nm_1_rm=0.58806-0.000432*T;
km_rm=0.00752-(0.000038)*T;
n_m=nm_1_rm+i*km_rm;

mg1=0.164+0.0462*exp(T/13.441); %-0.2
mg2=0.35+0.757*exp(T/3.311); %+0.08;

n_b=real(sqrt(e_b))+i*imag(sqrt(e_b));
n_t=real(sqrt(e_t))+i*imag(sqrt(e_t));
n_iw=real(sqrt(e_iw))+i*imag(sqrt(e_iw));

%n_m=0.567+i*0.007;

n_s=1+rd*(n_m+(n_b-1).*(mg+(mg1-mg)*st(mg-mg1))+...
    (n_t-1).*((mg-mg1)*st(mg-mg1)+(mg2-mg)*st(mg-mg2))+ ...
     ((n_iw-1)/riw).*(mg-mg2)*st(mg-mg2) );

else

e0b=13.66; %T>Tf
e0u=86.41; %T>Tf


betta0b=-2.87e-3;  %T>Tf
betta0u=0.12e-3;  %T>Tf


T0b=0; %T>Tf
T0u=20; %T>Tf


eoob=9.29; %T>Tf
eoou=12.13; %T>Tf

bettaoob=-4.0e-3;  %T>Tf
bettaoou=-0.15e-3;  %T>Tf


Toob=20; %T>Tf
Toou=20; %T>Tf


dHb=1540; %T>Tf
dHu=2128; %T>Tf


dSb=0.1; %T>Tf
dSu=2.86; %T>Tf


bettasigb=3.76e-3;  %T>Tf
bettasigu=9.45e-3;  %T>Tf


Tsigb=20; %T>Tf
Tsigu=0; %T>Tf


sigb=0.23; %T>Tf
sigu=1.36; %T>Tf


tau_b=48*exp(dHb./Tk-dSb)./Tk;
tau_u=48*exp(dHu./Tk-dSu)./Tk;



sig_b=sigb+bettasigb*(T-Tsigb);
sig_u=sigu+bettasigu*(T-Tsigu);


F0b=log((e0b-1.)/(e0b+2.));
F0u=log((e0u-1.)/(e0u+2.));


e0_b=(1+2*exp(F0b-betta0b*(T-T0b)))./(1-exp(F0b-betta0b*(T-T0b)));
e0_u=(1+2*exp(F0u-betta0u*(T-T0u)))./(1-exp(F0u-betta0u*(T-T0u)));


Foob=log((eoob-1.)/(eoob+2.));
Foou=log((eoou-1.)/(eoou+2.));


eoo_b=(1+2*exp(Foob-bettaoob*(T-Toob)))./(1-exp(Foob-bettaoob*(T-Toob)));
eoo_u=(1+2*exp(Foou-bettaoou*(T-Toou)))./(1-exp(Foou-bettaoou*(T-Toou)));


e_b=eoo_b+(e0_b-eoo_b)./(1-2*i*pi*freq*tau_b*1.e-3)+i*sig_b/(2*pi*freq*e0*1.d-3);
e_u=eoo_u+(e0_u-eoo_u)./(1-2*i*pi*freq*tau_u*1.e-3)+i*sig_u/(2*pi*freq*e0*1.d-3);


nm_1_rm=0.52181-0.00221*T;
km_rm=0.00684-0.000048*T;

n_m=nm_1_rm+i*km_rm;


mg1=0.324-0.0651*exp(T/33.2);


n_b=real(sqrt(e_b))+i*imag(sqrt(e_b));
n_u=real(sqrt(e_u))+i*imag(sqrt(e_u));


n_s=1+rd*(n_m+(n_b-1).*(mg+(mg1-mg)*st(mg-mg1))+...
        (n_u-1).*(mg-mg1).*st(mg-mg1));

end
end



function y=st(x)

if(x>0)
    y=1.0;
else
    y=0.0;
end 
    

end