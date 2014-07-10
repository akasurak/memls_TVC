%
%  [er,ei]=higheps(theta,f,sand,clay,rhob,rhos)  
%
%  Computes the er and ei given sand and clay fraction and frequency
%  er and ei are the real and imaginary part of complex dielectrics, resp.
%
%    theta              Volumetric water content 
%    f                  frequency (Hz)
%    sand               sand fraction (e.g. 0.50, not 50%)
%    clay               clay fraction
%    rhob               Bulk soil density (typically about 1.6)
%    rhos               Density of soil particles (typically about 2.6)
%
%
%  This function implements the model of Dobson et al. (1985) for the
%  dielectric constant of soils at frequencies of 1.4 Ghz to 18 Ghz.
%  The coefficients from the later paper by Peplinski et al. (with 
%  errata!) are used.  These differ slightly from the original paper.
%
function [er,ei]=higheps(theta,f,sand,clay,rhob,rhos)
TwoPitauw=0.58e-10;
ew0=80.1;
e0=8.85e-12;
ewinf=4.9;
a=0.65;
beta1=1.2748-0.519*sand-0.152*clay;
beta2=1.33797-0.603*sand-0.166*clay;
sigeff=-1.645+1.939*rhob-2.25622*sand+1.594*clay;
efw1=ewinf+(ew0-ewinf)/(1+(TwoPitauw*f)^2);
efw2=TwoPitauw*f*(ew0-ewinf)/(1+(TwoPitauw*f)^2)+sigeff*(rhos-rhob)/(2*pi*e0*f*rhos*theta);
es=((1.01+0.44*rhos)^2)-0.062;
er=(1+rhob*(es^a-1)/rhos+theta^beta1*efw1^a-theta)^(1/a);
ei=(theta^beta2*efw2^a)^(1/a);

